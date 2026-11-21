# -*- coding: utf-8 -*-
"""
03_ortho_clip_mosaic_v2.2.0.py
====================================
ZY1-02E AHSI 步骤3（方案A：先光谱后几何，v2.2.0）：RPC 正射 + 按研究区裁剪 + 云感知镶嵌

本步把步骤2 的**原始格网反射率**变成 5 个研究区（JLHK / LSK / LXJ / ZX / SC）
可以直接进入步骤4（蚀变矿物提取）、步骤5（机器学习）的“分析就绪”数据：

    逐景反射率（原始格网 + 校正后 RPC）
        ├─ RPC 正射到统一网格（EPSG/分辨率从步骤2 georef 读，默认 32650/30 m）
        ├─ 按 研究区 ∩ 单景足迹 裁剪 → 逐景瓦片（tile）
        └─ 按区镶嵌（清晰像元优先、云像元兜底）→ 区域 mosaic

两类成果都保留：
    tiles/   裁剪后的数据（一等成果，可人工重拼）
    mosaic/  拼接后的数据（派生成果，可从 tiles 随时重建）

每区输出（OUT_ROOT_STEP3/<区名>/）：
    tiles/<scene>_clip.tif           逐景裁剪反射率（N 波段=步骤2 剔除后，Float32，带 nodata）
    tiles/<scene>_cloud_mask.tif     同网格云掩膜（0=清晰 1=云 2=云影 255=无效）
    mosaic/mosaic_reflectance.tif    区域镶嵌反射率（N 波段，EPSG 从 georef 读取，默认 32650/30 m）
    mosaic/mosaic_mask.tif           0=无数据 1=清晰 2=云填充 3=云影填充
    mosaic/mosaic_source_scene.tif   每像元来源景编号（0=无）
    mosaic/mosaic_obs.tif            2 波段：总观测数 / 清晰观测数
    mosaic/mosaic_info.json          参与景、日期、优先级、覆盖率、云量、网格
    mosaic/wavelengths_clean.npy     波段中心波长（与步骤2 完全一致）
    preview_quality.png              真彩色/假彩色/NDVI/镶嵌掩膜 质检图

_catalog/ 下输出：
    tiles_catalog.csv    逐瓦片清单（景/区/日期/重叠率/云量/路径）
    area_coverage.csv    逐区汇总（覆盖率/清晰率/云填充率/tiles 数/日期）
    scene_status.csv     全部景中哪些已完成步骤2、哪些 pending
    area_grids.json      5 个研究区的输出网格定义（便于人工重拼/复核）
    step3_batch.log      批处理日志

处理要点
--------
1. 批量：扫描步骤2 输出目录，逐景处理；缺失的景只记为 pending，不报错。
2. 网格：每个研究区固定一个统一网格（EPSG/res 来自 georef，包围盒按 res 取整），
   所有瓦片和 mosaic 都严格对齐到该网格，瓦片可在 QGIS/GDAL 中直接叠置重拼。
3. 正射：使用步骤1 校正后的 RPC（tif 的 RPC / georef.json 的 rpc_corrected），
   不能用原始 .rpb（有约 2.4% 格网偏差）；近天底观测，无需 DEM 做正射。
4. 裁剪：按“研究区（外扩 AREA_BUFFER_M）∩ 单景足迹”的交集裁剪，节省磁盘；
   瓦片包围盒同样对齐到区域网格。
5. 云感知镶嵌：逐像元“先清晰、后云、同等级按优先级”合成：
       - 任何一景的清晰像元优先于有云像元；
       - 全部有云时按优先级用云像元填充（mask 标 2/3，第四步可据此过滤）；
       - 输出 source_scene 索引，第四/五步可按日期/景筛选样本。
   云掩膜来源：auto（AHSI 可见光-NIR 阈值）/ external（外源 MODIS/FY/手工）/ none。
6. 增量：tiles/cloud mask 是缓存；新增景后只补瓦片，`--stage mosaic` 重建 mosaic，
   不需要重跑正射。验证阶段只有 1 景时，允许部分覆盖（默认不要求 100%）。
7. 稳健性：原子写出、跳过已完成、默认出错即停（--continue-on-error 可显式关闭）、
   面积/网格/波段一致性校验、缺失云掩膜时明确降级或报错、输出 mask/来源景便于追溯。
   坐标真源：目标 CRS/分辨率从步骤2 georef 读取（CLI --epsg/--res 可覆盖），
   不再硬编码 32650/30；多景 georef 不一致即停，避免混合网格。

用法
----
    python 03_orthorectify_clip_mosaic_batch.py --dry-run
    python 03_orthorectify_clip_mosaic_batch.py                     # 全流程
    python 03_orthorectify_clip_mosaic_batch.py --stage tiles       # 只补瓦片
    python 03_orthorectify_clip_mosaic_batch.py --stage mosaic      # 只重拼
    python 03_orthorectify_clip_mosaic_batch.py --stage cloud       # 只重算云掩膜
    python 03_orthorectify_clip_mosaic_batch.py --area LSK --scene "*1104951*"

依赖：GDAL >= 3.0（Warp/RPC/OGR）、numpy；matplotlib 可选（质检图）；
      建议用装好 GDAL 的 Anaconda 解释器，例如
      C:\\ProgramData\\anaconda3\\python.exe

后续约定：步骤4 读取 mosaic/mosaic_reflectance.tif + mosaic_mask.tif +
wavelengths_clean.npy，做光谱指数/连续统去除/SAM/MF 等蚀变矿物提取；
步骤5 用同一套 cube + source_scene 索引做机器学习分类。
"""

import argparse
import csv
import datetime
import fnmatch
import glob
import json
import logging
import os
import re
import shutil
import sys
import time

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "当前 Python 解释器缺少 numpy：%s\n"
        "请用装好 GDAL 的解释器运行本脚本，例如：\n"
        "    C:\\ProgramData\\anaconda3\\python.exe 03_orthorectify_clip_mosaic_batch.py"
        % exc)

try:
    from osgeo import gdal, ogr, osr
except ImportError as exc:
    raise SystemExit(
        "当前 Python 解释器缺少 GDAL（osgeo）：%s\n"
        "请用装好 GDAL 的解释器运行本脚本，例如：\n"
        "    C:\\ProgramData\\anaconda3\\python.exe 03_orthorectify_clip_mosaic_batch.py"
        % exc)

gdal.UseExceptions()
ogr.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")

# v2.2.0：路径/网格统一收口到仓库根 config.py
import os as _os, sys as _sys
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _REPO_ROOT not in _sys.path:
    _sys.path.insert(0, _REPO_ROOT)
import config as _C

# ╔══════════════════════════════════════════════════════════════╗
# ║  【配置区】按需修改                                            ║
# ╚══════════════════════════════════════════════════════════════╝

# ────────────【输入/输出】步骤2 成果目录、本步成果根目录 ────────────
IN_ROOT_STEP2 = _C.OUT_ROOT_STEP2
OUT_ROOT_STEP3 = _C.OUT_ROOT_STEP3

# ────────────【输入】研究区矢量目录（每个子目录一个 .shp） ────────────
AREA_ROOT = _C.AREA_ROOT

# 全部原始景目录：仅用于统计“全部景里哪些还没完成步骤2（pending）”
EXPECTED_SCENES_DIR = _C.SRC_DIR

# ────────────【参数】输出网格 ────────────
TARGET_EPSG = _C.TARGET_EPSG  # 与步骤1/2 的 georef target_crs_epsg 保持一致
TARGET_RES = _C.TARGET_RES    # 与步骤1/2 的 georef target_res 保持一致
RESAMPLING = _C.RESAMPLING    # 反射率重采样；NEAREST 只建议用于掩膜
AREA_BUFFER_M = 100.0        # 研究区外扩（吸收 Xian80→WGS84 的基准面差异）
MIN_OVERLAP_PERCENT = 1.0    # 单景与研究区重叠低于该比例则不生成瓦片

# ────────────【参数】无效值 ────────────
NODATA_OUT = -9999.0

# v2.2.0：镶嵌前确定性剔除的“完全无用”波段——仅水汽强吸收带尾段：
#   * 1.4 μm 水汽带尾段 1.42~1.48 μm（1.4253/1.4422/1.4590/1.4759）
#   * 1.8 μm 水汽带前缘 1.79~1.82 μm（1.7954/1.8120）
# 这些波段透过率极低、无矿物信息，且多数景 DN=0（NoData 严重），直接剔除。
# 其余“有用波段”（2.0~2.5 μm 等）保留，其残留 NoData 由光谱方向补值（见 build_area_mosaic）。
DROP_BAD_WL_UM = [1.4253, 1.4422, 1.4590, 1.4759, 1.7954, 1.8120]

# ────────────【参数】S-G 光谱平滑（v2.2.0）────────────
# 沿光谱维对镶嵌反射率做 Savitzky-Golay 平滑（只处理研究区内像元，bbox 留白不动）。
# AHSI SWIR 采样约 16.7 nm，2200 nm Al-OH 吸收峰半高宽仅 1~2 个波段，窗口取 3 最保峰；
# 窗越大越平滑但越易抹平窄吸收峰（5 点轻微、7 点明显）。
SGO_WINDOW = 3                        # 平滑窗口（奇数，3/5/7）
SGO_POLYORDER = 2                     # 多项式阶（2/3）
SGO_PEAK_DEPTH_MAX_DROP = 0.10        # 平滑后诊断峰深度相对下降上限（超限打告警）
SGO_ENABLED = True                    # 是否启用 S-G 平滑

# ────────────【参数】云掩膜 ────────────
# auto     : 用 AHSI 可见光-NIR 阈值自动检测（无需外部数据）
# external : 读取外源云掩膜（MODIS/FY/人工），找不到时按 CLOUD_MASK_FALLBACK_AUTO 降级
# none     : 不做云检测（已知全无云）
CLOUD_MASK_MODE = "auto"
CLOUD_MASK_DIR = ""          # 外源掩膜目录：<CLOUD_MASK_DIR>/<scene>.tif
CLOUD_MASK_FALLBACK_AUTO = True

# 自动云检测的波段目标（μm）与阈值；阈值需在验证景上复核
# 只用可见光-NIR：SWIR 在步骤2 有较多 DN=0 造成的 nodata，会误伤云检测有效性
CLOUD_BAND_TARGETS = {
    "BLUE": 0.47, "GREEN": 0.55, "RED": 0.66, "NIR": 0.86,
}
CLOUD_THRESHOLDS = {
    "blue_min": 0.32,
    "red_min": 0.28,
    "brightness_min": 0.30,
    "ndvi_max": 0.18,
    "whiteness_min": 0.78,
    "saturation_value": 0.99,
    "saturation_ratio_min": 0.50,
}
# 云影（可选）：默认关闭；开启后用“暗像元”粗判，容易把水体误判，需人工复核
CLOUD_SHADOW_ENABLE = False
CLOUD_SHADOW_THRESHOLDS = {
    "blue_max": 0.05,
    "nir_max": 0.08,
    "brightness_max": 0.10,
}

# ────────────【参数】镶嵌优先级（清晰像元优先是固定规则，这里只排同级顺序） ──
# date_desc: 新日期优先 | date_asc: 老日期优先
# overlap_desc: 与研究区重叠大的景优先 | cloud_asc: 云量小的景优先
MOSAIC_PRIORITY = "date_desc"

# ────────────【参数】覆盖要求 ────────────
# 验证阶段只有 1 景，默认允许部分覆盖；批量多景时可开 --require-full-coverage
REQUIRE_FULL_COVERAGE = False
MIN_FULL_COVERAGE = 99.0

# ────────────【参数】批处理行为 ────────────
OVERWRITE = False            # True: 已有瓦片也重做
GDAL_CACHE_MB = 1024

GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]
GEOTIFF_COMPRESS = "DEFLATE"
GEOTIFF_COMPRESS_LEVEL = 6

# ────────────【参数】质检成果图 ────────────
MAKE_QC = True
QC_MAX_SIZE = 1600           # 预览图降采样后长边像素
QC_DPI = 150

# 文件名约定
IN_REFLECTANCE_SUFFIX = "_cal_badband_ortho_atmcorr.tif"
LEGACY_REFLECTANCE_SUFFIXES = (
    "_ortho_atmcorr.tif",
    "reflectance_atm.tif",
)
IN_WAVELENGTHS_NAME = "wavelengths_clean.npy"
IN_GEOREF_NAME = "georef.json"
IN_PARAMS_NAME = "used_parameters.json"

TILE_SUFFIX = "_clip.tif"
CLOUD_SUFFIX = "_cloud_mask.tif"

OUT_MOSAIC_NAME = "mosaic_reflectance.tif"
OUT_MASK_NAME = "mosaic_mask.tif"
OUT_SOURCE_NAME = "mosaic_source_scene.tif"
OUT_OBS_NAME = "mosaic_obs.tif"
OUT_INFO_NAME = "mosaic_info.json"
OUT_QC_NAME = "preview_quality.png"

LOG_NAME = "zy1e.step3"


# ╔══════════════════════════════════════════════════════════════╗
# ║  通用工具                                                     ║
# ╚══════════════════════════════════════════════════════════════╝

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="ZY1E AHSI 步骤3（方案A）：RPC 正射 + 研究区裁剪 + 云感知镶嵌")
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP2,
                        help="步骤2 成果根目录（默认: %(default)s）")
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP3,
                        help="本步成果根目录（默认: %(default)s）")
    parser.add_argument("--areas", dest="area_root", default=AREA_ROOT,
                        help="研究区矢量根目录（默认: %(default)s）")
    parser.add_argument("--area", action="append", default=None,
                        help="只处理匹配的研究区，可多次指定，如 LSK 或 \"2.*\"")
    parser.add_argument("--scene", action="append", default=None,
                        help="只处理匹配的场景名，可多次指定，支持通配符")
    parser.add_argument("--stage", default="all",
                        choices=("all", "tiles", "cloud", "mosaic"),
                        help="all=瓦片+镶嵌；tiles=只补瓦片；cloud=只重算云掩膜；"
                             "mosaic=只重拼（默认: all）")
    parser.add_argument("--priority", default=MOSAIC_PRIORITY,
                        choices=("date_desc", "date_asc", "overlap_desc", "cloud_asc"),
                        help="同级镶嵌优先级（默认: %(default)s）")
    parser.add_argument("--cloud-mode", default=CLOUD_MASK_MODE,
                        choices=("auto", "external", "none"),
                        help="云掩膜来源（默认: %(default)s）")
    parser.add_argument("--cloud-mask-dir", default=CLOUD_MASK_DIR,
                        help="外源云掩膜目录（<dir>/<scene>.tif）")
    parser.add_argument("--buffer-m", type=float, default=AREA_BUFFER_M,
                        help="研究区外扩米数（默认: %(default)s）")
    parser.add_argument("--min-overlap", type=float, default=MIN_OVERLAP_PERCENT,
                        help="景-区最小重叠百分比（默认: %(default)s）")
    parser.add_argument("--epsg", type=int, default=None,
                        help="输出目标 EPSG（默认从步骤2 georef.target_crs_epsg 读取，"
                             "缺失时用 %d）" % TARGET_EPSG)
    parser.add_argument("--res", type=float, default=None,
                        help="输出像元大小 m（默认从步骤2 georef.target_res 读取，"
                             "缺失时用 %g）" % TARGET_RES)
    parser.add_argument("--require-full-coverage", action="store_true",
                        default=REQUIRE_FULL_COVERAGE,
                        help="某区覆盖率低于阈值时视为失败（默认只告警）")
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE,
                        help="已有瓦片也重做（默认跳过）")
    parser.add_argument("--limit", type=int, default=0,
                        help="只处理前 N 景（默认 0 = 全部）")
    parser.add_argument("--dry-run", action="store_true",
                        help="只列出计划处理的景/区/瓦片，不执行")
    parser.add_argument("--no-qc", action="store_true", help="跳过质检图")
    parser.add_argument("--continue-on-error", action="store_true",
                        default=False,
                        help="单景/单区失败时继续（默认关闭=出错即停；"
                             "无视错误会导致数据不可用，仅在明确知情时开启）")
    parser.add_argument("--log", default=None,
                        help="日志文件（默认 <out>/step3_batch.log）")
    parser.add_argument("--verbose", action="store_true", help="输出调试日志")
    return parser.parse_args(argv)


def setup_logging(log_path, verbose=False):
    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(fmt)
    logger.addHandler(console)
    if log_path:
        try:
            parent = os.path.dirname(os.path.abspath(log_path))
            if parent:
                os.makedirs(parent, exist_ok=True)
            file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(fmt)
            logger.addHandler(file_handler)
        except OSError as exc:
            console.handle(logging.LogRecord(
                LOG_NAME, logging.WARNING, __file__, 0,
                "无法写入日志文件 %s: %s", (log_path, exc), None))
    logger.propagate = False
    return logger


def human_size(num_bytes):
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024.0 or unit == "TB":
            return "%.1f %s" % (size, unit)
        size /= 1024.0
    return "%.1f TB" % size


def fmt_elapsed(seconds):
    seconds = max(0.0, float(seconds))
    return "%d:%04.1f" % (int(seconds // 60), seconds % 60.0)


def atomic_write_json(payload, path):
    tmp_path = path + ".tmp"
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    os.replace(tmp_path, path)


def atomic_write_text(text, path):
    tmp_path = path + ".tmp"
    with open(tmp_path, "w", encoding="utf-8") as handle:
        handle.write(text)
    os.replace(tmp_path, path)


def safe_area_name(name):
    return re.sub(r"\s+", "", str(name))


def parse_scene_date(scene_name):
    match = re.search(r"_(\d{8})_", scene_name)
    if match:
        try:
            return datetime.datetime.strptime(match.group(1), "%Y%m%d")
        except ValueError:
            return None
    return None


def nearest_band(wavelengths, target):
    return int(np.argmin(np.abs(np.asarray(wavelengths, dtype=np.float64) - target)))


def build_gtiff_creation_options(compress=True):
    options = list(GTIFF_OPTIONS)
    options.append("TILED=YES")
    options.append("INTERLEAVE=BAND")
    if compress and GEOTIFF_COMPRESS:
        options.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS,
                        "PREDICTOR=3",
                        "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL])
    return options


def build_byte_creation_options():
    return ["BIGTIFF=IF_SAFER", "TILED=YES", "COMPRESS=DEFLATE", "ZLEVEL=6"]


# ╔══════════════════════════════════════════════════════════════╗
# ║  坐标系 / 几何工具                                             ║
# ╚══════════════════════════════════════════════════════════════╝

def make_srs(epsg, traditional=True):
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(int(epsg))
    if traditional:
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def transform_geometry(geom, src_epsg, dst_epsg):
    if geom is None:
        return None
    src = make_srs(src_epsg)
    dst = make_srs(dst_epsg)
    transformer = osr.CoordinateTransformation(src, dst)
    result = geom.Clone()
    result.Transform(transformer)
    return result


def footprint_geometry(footprint):
    """由 georef.json 的四角经纬度构造 WGS84 多边形。"""
    order = [key for key in ("TL", "TR", "BR", "BL") if key in footprint]
    if len(order) < 3:
        order = list(footprint.keys())
    ring = [(float(footprint[key][0]), float(footprint[key][1])) for key in order]
    if not ring:
        return None
    ring.append(ring[0])
    wkt = "POLYGON((%s))" % ", ".join("%.10f %.10f" % point for point in ring)
    return ogr.CreateGeometryFromWkt(wkt)


def write_geojson_cutline(geom_wgs84, path):
    """把 WGS84 几何写成 GeoJSON（GDAL 会按 WGS84 解释并重投影到目标 CRS）。"""
    geojson = {
        "type": "FeatureCollection",
        "features": [{"type": "Feature", "properties": {},
                      "geometry": json.loads(geom_wgs84.ExportToJson())}],
    }
    gdal.FileFromMemBuffer(path, json.dumps(geojson).encode("utf-8"))


def grid_from_bbox(bbox, res):
    xmin, xmax, ymin, ymax = [float(v) for v in bbox]
    x0 = np.floor(xmin / res) * res
    x1 = np.ceil(xmax / res) * res
    y0 = np.ceil(ymax / res) * res
    y1 = np.floor(ymin / res) * res
    width = int(round((x1 - x0) / res))
    height = int(round((y0 - y1) / res))
    if width <= 0 or height <= 0:
        raise ValueError("研究区包围盒无效: %s" % (bbox,))
    return {"x0": float(x0), "y0": float(y0), "res": float(res),
            "width": width, "height": height}


def rasterize_polygon(geom, grid):
    """把研究区多边形栅格化到区域网格，用于计算“面积内覆盖率”。"""
    srs = make_srs(TARGET_EPSG)
    memory_driver = ogr.GetDriverByName("Memory")
    source = memory_driver.CreateDataSource("area_mask")
    layer = source.CreateLayer("area", srs=srs, geom_type=ogr.wkbUnknown)
    feature = ogr.Feature(layer.GetLayerDefn())
    feature.SetGeometry(geom)
    layer.CreateFeature(feature)
    raster = gdal.GetDriverByName("MEM").Create(
        "", grid["width"], grid["height"], 1, gdal.GDT_Byte)
    raster.SetGeoTransform((grid["x0"], grid["res"], 0.0,
                            grid["y0"], 0.0, -grid["res"]))
    raster.SetProjection(srs.ExportToWkt())
    gdal.RasterizeLayer(raster, [1], layer, burn_values=[1])
    mask = raster.GetRasterBand(1).ReadAsArray() == 1
    raster = None
    source = None
    return mask


# ╔══════════════════════════════════════════════════════════════╗
# ║  研究区加载                                                    ║
# ╚══════════════════════════════════════════════════════════════╝

def load_areas(area_root, target_epsg, buffer_m, filters, logger):
    if not os.path.isdir(area_root):
        raise NotADirectoryError("研究区目录不存在: %s" % area_root)
    shapefiles = sorted(glob.glob(os.path.join(area_root, "*", "*.shp")))
    if not shapefiles:
        raise FileNotFoundError("研究区目录下未找到 *.shp: %s" % area_root)

    # v2.2.0：按目录分组，每个研究区目录优先用"与目录名（去掉数字/点前缀）同名"的主 shp，
    # 并 union 该 shp 内全部要素，避免旧版只取第一个要素导致多面研究区漏算。
    by_dir = {}
    for shp in shapefiles:
        by_dir.setdefault(os.path.dirname(shp), []).append(shp)

    areas = []
    for dirpath in sorted(by_dir):
        label = os.path.basename(os.path.normpath(dirpath))
        key = safe_area_name(label)
        letters = re.sub(r"[^A-Za-z]+", "", label)
        shps = by_dir[dirpath]
        main_shp = next(
            (s for s in shps
             if re.sub(r"[^A-Za-z]+", "", os.path.splitext(os.path.basename(s))[0]) == letters),
            shps[0])

        source_ds = ogr.Open(main_shp)
        if source_ds is None:
            logger.warning("无法打开研究区矢量: %s", main_shp)
            continue
        layer = source_ds.GetLayer(0)
        srs = layer.GetSpatialRef()
        if srs is None:
            logger.warning("研究区矢量缺少投影: %s", main_shp)
            source_ds = None
            continue
        code = srs.GetAuthorityCode(None)
        try:
            code = int(code)
        except (TypeError, ValueError):
            logger.warning("研究区矢量无法识别 EPSG: %s", main_shp)
            source_ds = None
            continue

        union_geom = None
        feature = layer.GetNextFeature()
        while feature is not None:
            g = feature.GetGeometryRef()
            if g is not None and not g.IsEmpty():
                union_geom = g.Clone() if union_geom is None else union_geom.Union(g)
            feature = layer.GetNextFeature()
        source_ds = None
        if union_geom is None or union_geom.IsEmpty():
            logger.warning("研究区矢量缺少有效要素: %s", main_shp)
            continue

        geom_target = transform_geometry(union_geom, code, target_epsg)
        if buffer_m and buffer_m > 0:
            geom_target = geom_target.Buffer(float(buffer_m))
        if geom_target is None or geom_target.IsEmpty():
            logger.warning("研究区投影/缓冲后为空: %s", main_shp)
            continue
        geom_wgs84 = transform_geometry(geom_target, target_epsg, 4326)

        grid = grid_from_bbox(geom_target.GetEnvelope(), TARGET_RES)
        area = {
            "key": key,
            "label": label,
            "shp": os.path.abspath(main_shp),
            "source_epsg": code,
            "geom": geom_target,
            "geom_wgs84": geom_wgs84,
            "grid": grid,
            "area_km2": float(geom_target.GetArea()) / 1e6,
        }
        if filters and not _match_any(area["key"], filters) \
                and not _match_any(area["label"], filters):
            continue
        areas.append(area)

    if not areas:
        raise FileNotFoundError("没有匹配的研究区（--area %s）"
                                % (", ".join(filters) if filters else ""))
    logger.info("研究区: %d 个", len(areas))
    for area in areas:
        grid = area["grid"]
        logger.info("  %-8s EPSG:%-5d 面积 %7.1f km2  网格 %d×%d @ %g m",
                    area["key"], area["source_epsg"], area["area_km2"],
                    grid["width"], grid["height"], grid["res"])
    return areas


def _match_any(text, patterns):
    text = str(text).lower()
    for pattern in patterns:
        pattern = str(pattern).lower()
        if fnmatch.fnmatch(text, pattern) or pattern in text:
            return True
    return False


# ╔══════════════════════════════════════════════════════════════╗
# ║  场景加载与景-区重叠                                           ║
# ╚══════════════════════════════════════════════════════════════╝

def discover_scenes(in_root):
    if not os.path.isdir(in_root):
        raise NotADirectoryError("步骤2 成果目录不存在: %s" % in_root)
    if resolve_reflectance_path(in_root):
        return [in_root]
    scene_dirs = sorted(glob.glob(os.path.join(in_root, "*")))
    return [path for path in scene_dirs
            if os.path.isdir(path) and resolve_reflectance_path(path)]


def resolve_reflectance_path(scene_dir):
    renamed = sorted(glob.glob(os.path.join(
        scene_dir, "*" + IN_REFLECTANCE_SUFFIX)))
    if renamed:
        return renamed[0]
    for suffix in LEGACY_REFLECTANCE_SUFFIXES:
        pattern_name = ("*" + suffix if suffix.startswith("_") else suffix)
        legacy = sorted(glob.glob(os.path.join(scene_dir, pattern_name)))
        if legacy:
            return legacy[0]
    return None


def read_scene_info(scene_dir):
    scene = os.path.basename(os.path.normpath(scene_dir))
    reflectance = resolve_reflectance_path(scene_dir)
    if reflectance is None:
        raise FileNotFoundError(
            "找不到步骤2反射率（*_cal_badband_ortho_atmcorr.tif "
            "或兼容旧格式成果）：%s" % scene_dir)
    wavelengths_path = os.path.join(scene_dir, IN_WAVELENGTHS_NAME)
    georef_path = os.path.join(scene_dir, IN_GEOREF_NAME)
    params_path = os.path.join(scene_dir, IN_PARAMS_NAME)

    wavelengths = np.load(wavelengths_path, allow_pickle=False)
    wavelengths = np.asarray(wavelengths, dtype=np.float64).ravel()
    if wavelengths.size == 0 or np.any(wavelengths <= 0):
        raise ValueError("波长数据异常: %s" % wavelengths_path)
    source_dataset = gdal.Open(reflectance)
    if source_dataset is None:
        raise RuntimeError("无法打开步骤2 反射率: %s" % reflectance)
    if source_dataset.RasterCount != wavelengths.size:
        raise ValueError("反射率波段数 %d != 波长数 %d（%s）"
                         % (source_dataset.RasterCount, wavelengths.size,
                            reflectance))
    source_dataset = None

    georef = None
    if os.path.isfile(georef_path):
        with open(georef_path, "r", encoding="utf-8") as handle:
            georef = json.load(handle)
    params = None
    if os.path.isfile(params_path):
        with open(params_path, "r", encoding="utf-8") as handle:
            params = json.load(handle)

    date_text = None
    if params and params.get("date"):
        date_text = str(params["date"])[:10]
    if not date_text:
        scene_date = parse_scene_date(scene)
        date_text = scene_date.strftime("%Y-%m-%d") if scene_date else None
    aot = params.get("aot550") if params else None
    elevation = params.get("elevation_km") if params else None

    footprint = None
    footprint_geom = None
    if georef and georef.get("footprint_wgs84"):
        footprint = georef["footprint_wgs84"]
        footprint_geom = footprint_geometry(footprint)

    rpc = None
    if georef and georef.get("rpc_corrected"):
        rpc = {str(key): str(value) for key, value in georef["rpc_corrected"].items()}

    # 目标坐标真源：优先用步骤2 透传的 georef（单一真源），否则留空由 main 用默认
    target_epsg = None
    target_res = None
    if georef:
        if georef.get("target_crs_epsg") is not None:
            try:
                target_epsg = int(georef["target_crs_epsg"])
            except (TypeError, ValueError):
                target_epsg = None
        if georef.get("target_res") is not None:
            try:
                target_res = float(georef["target_res"])
            except (TypeError, ValueError):
                target_res = None
    if params:
        if target_epsg is None and params.get("target_crs_epsg") is not None:
            try:
                target_epsg = int(params["target_crs_epsg"])
            except (TypeError, ValueError):
                pass
        if target_res is None and params.get("target_res") is not None:
            try:
                target_res = float(params["target_res"])
            except (TypeError, ValueError):
                pass

    return {
        "scene": scene,
        "dir": scene_dir,
        "reflectance": reflectance,
        "wavelengths": wavelengths,
        "n_bands": int(wavelengths.size),
        "footprint": footprint,
        "footprint_geom": footprint_geom,
        "rpc": rpc,
        "date": date_text,
        "aot550": float(aot) if aot is not None else None,
        "elevation_km": float(elevation) if elevation is not None else None,
        "target_epsg": target_epsg,           # v2.1.0：来自 georef 的目标 CRS
        "target_res": target_res,             # v2.1.0：来自 georef 的目标分辨率
    }


def scene_area_overlap(scene, area, min_overlap):
    if scene["footprint_geom"] is None:
        return None
    footprint = transform_geometry(scene["footprint_geom"], 4326, TARGET_EPSG)
    if footprint is None or footprint.IsEmpty():
        return None
    intersection = footprint.Intersection(area["geom"])
    if intersection is None or intersection.IsEmpty():
        return None
    area_value = area["geom"].GetArea()
    if area_value <= 0:
        return None
    percent = 100.0 * intersection.GetArea() / area_value
    if percent < min_overlap:
        return None
    return {"geometry": intersection, "percent": float(percent)}


def tile_window(area, intersection_geom):
    """由交集几何得到对齐到区域网格的瓦片窗口与 outputBounds。"""
    grid = area["grid"]
    res = grid["res"]
    xmin, xmax, ymin, ymax = intersection_geom.GetEnvelope()
    col0 = int(np.floor((xmin - grid["x0"]) / res))
    col1 = int(np.ceil((xmax - grid["x0"]) / res))
    row0 = int(np.floor((grid["y0"] - ymax) / res))
    row1 = int(np.ceil((grid["y0"] - ymin) / res))
    col0 = max(0, min(col0, grid["width"]))
    col1 = max(0, min(col1, grid["width"]))
    row0 = max(0, min(row0, grid["height"]))
    row1 = max(0, min(row1, grid["height"]))
    if col1 <= col0 or row1 <= row0:
        return None
    width = col1 - col0
    height = row1 - row0
    x_left = grid["x0"] + col0 * res
    y_top = grid["y0"] - row0 * res
    bounds = (x_left, y_top - height * res, x_left + width * res, y_top)
    return {"bounds": bounds, "width": width, "height": height,
            "col_off": col0, "row_off": row0}


# ╔══════════════════════════════════════════════════════════════╗
# ║  云掩膜                                                        ║
# ╚══════════════════════════════════════════════════════════════╝

def cloud_mask_path(tile_path):
    if tile_path.endswith(TILE_SUFFIX):
        return tile_path[:-len(TILE_SUFFIX)] + CLOUD_SUFFIX
    return tile_path + CLOUD_SUFFIX


def _write_uint8_mask(mask, reference_path, output_path, description,
                      nodata=255, metadata=None):
    driver = gdal.GetDriverByName("GTiff")
    tmp_path = output_path + ".tmp.tif"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    reference = gdal.Open(reference_path)
    width, height = reference.RasterXSize, reference.RasterYSize
    geo_transform = reference.GetGeoTransform()
    projection = reference.GetProjection()
    reference = None
    dataset = driver.Create(tmp_path, width, height, 1, gdal.GDT_Byte,
                            options=build_byte_creation_options())
    dataset.SetGeoTransform(geo_transform)
    dataset.SetProjection(projection)
    band = dataset.GetRasterBand(1)
    # 先设置 nodata 再写数据：GDAL/GTiff 在“先写后设 nodata”时可能把全零块
    # 当成空块并用 nodata 填充，导致“全清晰”掩膜被误写成全 255
    band.SetNoDataValue(int(nodata))
    band.WriteArray(mask.astype(np.uint8))
    band.SetDescription(description)
    if metadata:
        for key, value in metadata.items():
            dataset.SetMetadataItem(str(key), str(value))
    dataset.FlushCache()
    dataset = None
    os.replace(tmp_path, output_path)


def compute_auto_cloud_mask(tile_path, wavelengths, mask_path, logger):
    """基于 AHSI 可见光-NIR 阈值生成云掩膜（0=清晰 1=云 2=云影 255=无效）。"""
    dataset = gdal.Open(tile_path)
    if dataset is None:
        raise RuntimeError("无法打开瓦片: %s" % tile_path)
    width, height = dataset.RasterXSize, dataset.RasterYSize
    arrays = {}
    for name, target in CLOUD_BAND_TARGETS.items():
        index = nearest_band(wavelengths, target) + 1
        arrays[name] = dataset.GetRasterBand(index).ReadAsArray().astype(np.float64)
    dataset = None

    valid = np.ones((height, width), dtype=bool)
    for values in arrays.values():
        valid &= np.isfinite(values) & (values != NODATA_OUT)
    blue = arrays["BLUE"]
    green = arrays["GREEN"]
    red = arrays["RED"]
    nir = arrays["NIR"]

    denominator = nir + red
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = np.where(denominator > 0, (nir - red) / denominator, 0.0)
    brightness = (blue + green + red) / 3.0
    stack = np.stack([blue, green, red, nir], axis=0)
    stack_max = np.maximum(stack.max(axis=0), 1e-6)
    whiteness = stack.min(axis=0) / stack_max
    saturation_ratio = np.mean(stack >= CLOUD_THRESHOLDS["saturation_value"],
                               axis=0)

    cloud = valid & (
        ((blue >= CLOUD_THRESHOLDS["blue_min"])
         & (red >= CLOUD_THRESHOLDS["red_min"])
         & (brightness >= CLOUD_THRESHOLDS["brightness_min"])
         & (ndvi <= CLOUD_THRESHOLDS["ndvi_max"]))
        | ((whiteness >= CLOUD_THRESHOLDS["whiteness_min"])
           & (brightness >= CLOUD_THRESHOLDS["brightness_min"]))
        | (saturation_ratio >= CLOUD_THRESHOLDS["saturation_ratio_min"])
    )

    shadow = np.zeros_like(cloud)
    if CLOUD_SHADOW_ENABLE:
        shadow = (valid & ~cloud
                  & (blue <= CLOUD_SHADOW_THRESHOLDS["blue_max"])
                  & (nir <= CLOUD_SHADOW_THRESHOLDS["nir_max"])
                  & (brightness <= CLOUD_SHADOW_THRESHOLDS["brightness_max"]))

    mask = np.full((height, width), 255, dtype=np.uint8)
    mask[valid] = 0
    mask[shadow] = 2
    mask[cloud] = 1
    _write_uint8_mask(
        mask, tile_path, mask_path,
        "cloud mask (0=clear 1=cloud 2=shadow 255=nodata)",
        metadata={"CLOUD_MODE": "auto",
                  "THRESHOLDS": json.dumps(CLOUD_THRESHOLDS)},
    )
    stats = mask_statistics(mask)
    logger.info("  云掩膜: 清晰 %.1f%%, 云 %.1f%%, 云影 %.1f%%",
                stats["clear_percent"], stats["cloud_percent"],
                stats["shadow_percent"])
    return stats


def mask_statistics(mask):
    total = int(mask.size)
    valid = mask != 255
    valid_count = int(np.sum(valid))
    clear = int(np.sum(mask == 0))
    cloud = int(np.sum(mask == 1))
    shadow = int(np.sum(mask == 2))
    return {
        "valid_pixels": valid_count,
        "clear_pixels": clear,
        "cloud_pixels": cloud,
        "shadow_pixels": shadow,
        "total_pixels": total,
        "clear_percent": 100.0 * clear / valid_count if valid_count else 0.0,
        "cloud_percent": 100.0 * cloud / valid_count if valid_count else 0.0,
        "shadow_percent": 100.0 * shadow / valid_count if valid_count else 0.0,
        "nodata_percent": 100.0 * (total - valid_count) / total if total else 0.0,
    }


def find_external_cloud_mask(scene, cloud_dir):
    candidates = [os.path.join(scene["dir"], "cloud_mask.tif")]
    if cloud_dir:
        candidates.append(os.path.join(cloud_dir, scene["scene"] + ".tif"))
        candidates.append(os.path.join(cloud_dir, scene["scene"] + ".tiff"))
    for path in candidates:
        if path and os.path.isfile(path):
            return os.path.abspath(path)
    return None


def warp_external_cloud_mask(scene, external_path, area, window, cutline_path,
                             mask_path, logger):
    """把外源云掩膜重采样到瓦片网格（nearest）。"""
    source = gdal.Open(external_path)
    if source is None:
        raise RuntimeError("无法打开外源云掩膜: %s" % external_path)
    rpc = source.GetMetadata("RPC")
    geo_transform = source.GetGeoTransform()
    projection = source.GetProjection()
    source = None

    vrt = None
    use_rpc = False
    if rpc:
        use_rpc = True
        source_path = external_path
    elif projection and not _is_identity_geotransform(geo_transform):
        use_rpc = False
        source_path = external_path
    elif scene.get("rpc"):
        vrt = gdal.Translate("/vsimem/step3_ext_cloud.vrt", external_path,
                             format="VRT")
        vrt.SetMetadata(scene["rpc"], "RPC")
        source_path = "/vsimem/step3_ext_cloud.vrt"
        use_rpc = True
    else:
        source_path = external_path
        use_rpc = False

    try:
        options = gdal.WarpOptions(
            format="GTiff",
            dstSRS="EPSG:%d" % TARGET_EPSG,
            xRes=TARGET_RES, yRes=TARGET_RES,
            outputBounds=window["bounds"],
            targetAlignedPixels=True,
            resampleAlg="nearest",
            srcNodata=255, dstNodata=255,
            cutlineDSName=cutline_path, cropToCutline=True,
            multithread=True,
            warpOptions=["NUM_THREADS=ALL_CPUS"],
            creationOptions=build_byte_creation_options(),
            rpc=use_rpc,
            outputType=gdal.GDT_Byte,
        )
        tmp_path = mask_path + ".tmp.tif"
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        result = gdal.Warp(tmp_path, source_path, options=options)
        if result is None:
            raise RuntimeError("外源云掩膜 warp 失败: %s" % mask_path)
        result = None
        dataset = gdal.Open(tmp_path, gdal.GA_Update)
        band = dataset.GetRasterBand(1)
        band.SetDescription("cloud mask (0=clear 1=cloud 2=shadow 255=nodata)")
        dataset.SetMetadataItem("CLOUD_MODE", "external")
        dataset.SetMetadataItem("CLOUD_MASK_SOURCE", external_path)
        dataset = None
        os.replace(tmp_path, mask_path)
    finally:
        vrt = None
        if use_rpc and source_path.startswith("/vsimem/"):
            try:
                gdal.Unlink(source_path)
            except Exception:
                pass


def _is_identity_geotransform(geo_transform):
    if geo_transform is None:
        return True
    return (abs(geo_transform[0]) < 1e-9 and abs(geo_transform[1] - 1.0) < 1e-9
            and abs(geo_transform[2]) < 1e-9 and abs(geo_transform[3]) < 1e-9
            and abs(geo_transform[4]) < 1e-9 and abs(geo_transform[5] - 1.0) < 1e-9)


def build_tile_cloud_mask(scene, area, window, cutline_path, tile_path,
                          mask_path, cloud_mode, cloud_mask_dir, logger):
    external = None
    if cloud_mode == "external":
        external = find_external_cloud_mask(scene, cloud_mask_dir)
        if external is None:
            if CLOUD_MASK_FALLBACK_AUTO:
                logger.warning("场景 %s 未找到外源云掩膜，降级为 auto", scene["scene"])
                cloud_mode = "auto"
            else:
                raise FileNotFoundError(
                    "场景 %s 未找到外源云掩膜（%s）"
                    % (scene["scene"], cloud_mask_dir))
    if cloud_mode == "external":
        warp_external_cloud_mask(scene, external, area, window, cutline_path,
                                 mask_path, logger)
        dataset = gdal.Open(mask_path)
        mask = dataset.GetRasterBand(1).ReadAsArray()
        dataset = None
        return mask_statistics(mask)
    if cloud_mode == "none":
        dataset = gdal.Open(tile_path)
        band1 = dataset.GetRasterBand(1).ReadAsArray()
        dataset = None
        valid = np.isfinite(band1) & (band1 != NODATA_OUT)
        mask = np.full(band1.shape, 255, dtype=np.uint8)
        mask[valid] = 0
        _write_uint8_mask(mask, tile_path, mask_path,
                          "cloud mask (0=clear 255=nodata)", nodata=255,
                          metadata={"CLOUD_MODE": "none"})
        return mask_statistics(mask)
    return compute_auto_cloud_mask(tile_path, scene["wavelengths"],
                                   mask_path, logger)


# ╔══════════════════════════════════════════════════════════════╗
# ║  逐景瓦片（正射 + 裁剪）                                       ║
# ╚══════════════════════════════════════════════════════════════╝

def warp_reflectance_tile(scene, area, overlap, window, cutline_path,
                          tile_path, logger):
    source_path = scene["reflectance"]
    source = gdal.Open(source_path)
    if source is None:
        raise RuntimeError("无法打开步骤2 反射率: %s" % source_path)
    rpc = source.GetMetadata("RPC")
    source = None
    vrt = None
    use_rpc = True
    if not rpc:
        if not scene.get("rpc"):
            raise ValueError("步骤2 反射率缺少 RPC: %s" % source_path)
        vrt = gdal.Translate("/vsimem/step3_src_%s.vrt" % scene["scene"],
                             source_path, format="VRT")
        vrt.SetMetadata(scene["rpc"], "RPC")
        source_path = "/vsimem/step3_src_%s.vrt" % scene["scene"]
    try:
        options = gdal.WarpOptions(
            format="GTiff",
            dstSRS="EPSG:%d" % TARGET_EPSG,
            xRes=TARGET_RES, yRes=TARGET_RES,
            outputBounds=window["bounds"],
            targetAlignedPixels=True,
            resampleAlg=RESAMPLING,
            dstNodata=NODATA_OUT,
            cutlineDSName=cutline_path, cropToCutline=True,
            multithread=True,
            warpOptions=["NUM_THREADS=ALL_CPUS"],
            creationOptions=build_gtiff_creation_options(),
            rpc=use_rpc,
            outputType=gdal.GDT_Float32,
        )
        tmp_path = tile_path + ".tmp.tif"
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        result = gdal.Warp(tmp_path, source_path, options=options)
        if result is None:
            raise RuntimeError("RPC 正射失败: %s" % tile_path)
        result = None

        dataset = gdal.Open(tmp_path, gdal.GA_Update)
        for index, wavelength in enumerate(scene["wavelengths"], start=1):
            band = dataset.GetRasterBand(index)
            band.SetDescription("%.4f um" % float(wavelength))
        dataset.SetMetadataItem("SCENE", scene["scene"])
        dataset.SetMetadataItem("DATE", scene["date"] or "unknown")
        dataset.SetMetadataItem("AREA", area["key"])
        dataset.SetMetadataItem("OVERLAP_PERCENT", "%.3f" % overlap["percent"])
        if scene.get("aot550") is not None:
            dataset.SetMetadataItem("AOT550", "%.6f" % scene["aot550"])
        if scene.get("elevation_km") is not None:
            dataset.SetMetadataItem("ELEVATION_KM", "%.6f" % scene["elevation_km"])
        dataset.SetMetadataItem("TARGET_CRS", "EPSG:%d" % TARGET_EPSG)
        dataset.SetMetadataItem("TARGET_RES", "%g" % TARGET_RES)
        dataset.SetMetadataItem(
            "GEOMETRY", "RPC orthorectified (step1 corrected RPC)")
        dataset.SetMetadataItem("SOURCE_REFLECTANCE",
                                os.path.basename(scene["reflectance"]))
        dataset = None
        # 不在此处提交；build_tile 负责瓦片和掩膜成对提交
    finally:
        vrt = None
        if source_path.startswith("/vsimem/"):
            try:
                gdal.Unlink(source_path)
            except Exception:
                pass


def build_tile(scene, area, overlap, args, logger):
    window = tile_window(area, overlap["geometry"])
    if window is None:
        logger.warning("  %s 与 %s 交集窗口为空，跳过",
                       scene["scene"], area["key"])
        return "skip"
    area_dir = os.path.join(args.out_root, area["key"])
    tile_dir = os.path.join(area_dir, "tiles")
    os.makedirs(tile_dir, exist_ok=True)
    tile_path = os.path.join(tile_dir, scene["scene"] + TILE_SUFFIX)
    mask_path = cloud_mask_path(tile_path)

    need_tile = args.overwrite or not os.path.isfile(tile_path)
    need_mask = args.overwrite or not os.path.isfile(mask_path)
    # v2.2.0：瓦片波段数与当前景不一致、或上游反射率比瓦片新，均需重建瓦片
    if not need_tile:
        _stale_bands = False
        try:
            _tds = gdal.Open(tile_path)
            _stale_bands = _tds is None or _tds.RasterCount != scene["n_bands"]
            _tds = None
        except Exception:
            _stale_bands = True
        _stale_src = False
        try:
            _stale_src = os.path.getmtime(scene["reflectance"]) > os.path.getmtime(tile_path)
        except OSError:
            _stale_src = False
        if _stale_bands or _stale_src:
            need_tile = True
            need_mask = True
    if not need_tile and not need_mask:
        return "skip"

    tile_tmp = tile_path + ".tmp.tif"
    tile_backup = tile_path + ".bak"
    mask_backup = mask_path + ".bak"
    tile_had_backup = need_tile and os.path.isfile(tile_path)
    mask_had_backup = need_mask and os.path.isfile(mask_path)

    if os.path.exists(tile_tmp):
        os.remove(tile_tmp)
    if (tile_had_backup and os.path.exists(tile_backup)) \
            or (mask_had_backup and os.path.exists(mask_backup)):
        raise FileExistsError("检测到上次未完成的瓦片备份，请先人工检查 .bak 文件")

    cutline_path = "/vsimem/step3_cutline_%s_%s.geojson" % (
        area["key"].replace(".", "_"), scene["scene"][-12:])
    intersection_wgs84 = transform_geometry(overlap["geometry"], TARGET_EPSG, 4326)
    write_geojson_cutline(intersection_wgs84, cutline_path)
    try:
        if tile_had_backup:
            os.replace(tile_path, tile_backup)
        if mask_had_backup:
            os.replace(mask_path, mask_backup)
        if need_tile:
            logger.info("  正射裁剪: %s  [%s × %s] 重叠 %.1f%%",
                        scene["scene"][-16:], area["key"], scene["date"],
                        overlap["percent"])
            warp_reflectance_tile(scene, area, overlap, window, cutline_path,
                                  tile_path, logger)
            if not os.path.isfile(tile_tmp):
                raise RuntimeError("瓦片临时成果不存在: %s" % tile_tmp)
            os.replace(tile_tmp, tile_path)
        if need_mask:
            stats = build_tile_cloud_mask(
                scene, area, window, cutline_path, tile_path, mask_path,
                args.cloud_mode, args.cloud_mask_dir, logger)
            logger.info("  瓦片云量: 清晰 %.1f%%, 云 %.1f%%",
                        stats["clear_percent"], stats["cloud_percent"])
        if tile_had_backup and os.path.exists(tile_backup):
            os.remove(tile_backup)
        if mask_had_backup and os.path.exists(mask_backup):
            os.remove(mask_backup)
    except Exception:
        if os.path.exists(tile_tmp):
            os.remove(tile_tmp)
        if need_tile and os.path.exists(tile_path):
            os.remove(tile_path)
        if tile_had_backup and os.path.exists(tile_backup):
            os.replace(tile_backup, tile_path)
        if need_mask and os.path.exists(mask_path):
            os.remove(mask_path)
        if mask_had_backup and os.path.exists(mask_backup):
            os.replace(mask_backup, mask_path)
        raise
    finally:
        try:
            gdal.Unlink(cutline_path)
        except Exception:
            pass

    return "ok" if need_tile else "mask-only"


# ╔══════════════════════════════════════════════════════════════╗
# ║  瓦片扫描 / 镶嵌                                               ║
# ╚══════════════════════════════════════════════════════════════╝

def scan_tiles(area):
    area_dir = area["out_dir"]
    tile_dir = os.path.join(area_dir, "tiles")
    if not os.path.isdir(tile_dir):
        return []
    tiles = []
    for path in sorted(glob.glob(os.path.join(tile_dir, "*" + TILE_SUFFIX))):
        try:
            dataset = gdal.Open(path)
            if dataset is None:
                continue
            metadata = dataset.GetMetadata()
            geo_transform = dataset.GetGeoTransform()
            width, height = dataset.RasterXSize, dataset.RasterYSize
            raster_count = dataset.RasterCount
            dataset = None
            mask_path = cloud_mask_path(path)
            if not os.path.isfile(mask_path):
                continue
            grid = area["grid"]
            res = grid["res"]
            col_off = int(round((geo_transform[0] - grid["x0"]) / res))
            row_off = int(round((grid["y0"] - geo_transform[3]) / res))
            expected_col = (geo_transform[0] - grid["x0"]) / res
            expected_row = (grid["y0"] - geo_transform[3]) / res
            if abs(expected_col - col_off) > 1e-3 or abs(expected_row - row_off) > 1e-3:
                continue
            mask_stats = read_mask_stats(mask_path)
            tiles.append({
                "path": path,
                "mask_path": mask_path,
                "scene": metadata.get("SCENE", os.path.basename(path)),
                "date": metadata.get("DATE", ""),
                "overlap": float(metadata.get("OVERLAP_PERCENT", 0.0) or 0.0),
                "aot550": _to_float(metadata.get("AOT550")),
                "width": width,
                "height": height,
                "raster_count": raster_count,
                "col_off": col_off,
                "row_off": row_off,
                "cloud_percent": mask_stats["cloud_percent"],
                "clear_percent": mask_stats["clear_percent"],
                "valid_pixels": mask_stats["valid_pixels"],
                "total_pixels": mask_stats["total_pixels"],
            })
        except Exception:
            continue
    return tiles


def _to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def read_mask_stats(mask_path):
    dataset = gdal.Open(mask_path)
    mask = dataset.GetRasterBand(1).ReadAsArray()
    dataset = None
    return mask_statistics(mask)


def sort_tiles(tiles, priority):
    if priority == "date_asc":
        return sorted(tiles, key=lambda item: (item["date"] or "", item["scene"]))
    if priority == "overlap_desc":
        return sorted(tiles, key=lambda item: -item["overlap"])
    if priority == "cloud_asc":
        return sorted(tiles, key=lambda item: item["cloud_percent"])
    return sorted(tiles, key=lambda item: (item["date"] or "", item["scene"]),
                  reverse=True)


def mosaic_paths(area):
    mosaic_dir = os.path.join(area["out_dir"], "mosaic")
    return {
        "dir": mosaic_dir,
        "reflectance": os.path.join(mosaic_dir, OUT_MOSAIC_NAME),
        "mask": os.path.join(mosaic_dir, OUT_MASK_NAME),
        "source": os.path.join(mosaic_dir, OUT_SOURCE_NAME),
        "obs": os.path.join(mosaic_dir, OUT_OBS_NAME),
        "info": os.path.join(mosaic_dir, OUT_INFO_NAME),
        "qc": os.path.join(mosaic_dir, OUT_QC_NAME),
        "wavelengths": os.path.join(mosaic_dir, IN_WAVELENGTHS_NAME),
    }


def mosaic_is_current(area, tiles, wavelengths=None):
    paths = mosaic_paths(area)
    needed = [paths["reflectance"], paths["mask"], paths["source"],
              paths["obs"], paths["info"]]
    if not all(os.path.isfile(path) for path in needed):
        return False
    # v2.2.0：瓦片波段数与波长不一致 → 判为过期，强制重建，避免波段错配
    if wavelengths is not None:
        for tile in tiles:
            try:
                ds = gdal.Open(tile["path"])
                if ds is None or ds.RasterCount != int(len(wavelengths)):
                    ds = None
                    return False
                ds = None
            except Exception:
                return False
    try:
        mosaic_mtime = min(os.path.getmtime(path) for path in needed)
    except OSError:
        return False
    for tile in tiles:
        try:
            if os.path.getmtime(tile["path"]) > mosaic_mtime:
                return False
            if os.path.getmtime(tile["mask_path"]) > mosaic_mtime:
                return False
        except OSError:
            return False
    return True


def _scene_mean_spectrum(tif_path, area_mask):
    """逐波段计算研究区内有效像元的平均反射率，返回 (n_bands,) 数组。"""
    ds = gdal.Open(tif_path)
    n = ds.RasterCount
    means = np.full(n, np.nan, dtype=np.float64)
    for b in range(1, n + 1):
        arr = ds.GetRasterBand(b).ReadAsArray().astype(np.float64)
        m = (arr != NODATA_OUT) & area_mask
        if np.any(m):
            means[b - 1] = float(arr[m].mean())
    ds = None
    return means


def _diagnostic_peaks(wavelengths, spectrum):
    """诊断吸收峰（2200/2340/1000-1100nm）的峰位与深度，返回 {名称: (峰位nm, 深度)}。"""
    metrics = {}
    for lo, hi, name in [(2.15, 2.25, "AlOH_2200"), (2.28, 2.40, "CO3_2340"),
                         (1.00, 1.12, "Fe_1100")]:
        idx = [i for i, w in enumerate(wavelengths) if lo <= w <= hi]
        if len(idx) < 2:
            continue
        sub = spectrum[idx]
        if not np.all(np.isfinite(sub)):
            continue
        baseline = 0.5 * (float(sub[0]) + float(sub[-1]))
        if baseline <= 0:
            continue
        depth = 1.0 - float(sub.min()) / baseline
        peak_nm = float(wavelengths[idx[int(np.argmin(sub))]])
        metrics[name] = (peak_nm, depth)
    return metrics


def spectral_savgol_smooth(tif_path, wavelengths, area_mask, logger):
    """v2.2.0：沿光谱维对镶嵌反射率做 S-G 平滑（滑窗，只处理研究区内像元）。

    返回 (平滑前诊断峰, 平滑后诊断峰)，供 QC/留痕。
    """
    from scipy.signal import savgol_filter
    ds = gdal.Open(tif_path, gdal.GA_Update)
    n = ds.RasterCount
    before = _scene_mean_spectrum(tif_path, area_mask)
    peaks_before = _diagnostic_peaks(wavelengths, before)
    half = SGO_WINDOW // 2
    order = min(SGO_POLYORDER, SGO_WINDOW - 1)
    for b in range(1, n + 1):
        b0 = max(1, b - half)
        b1 = min(n, b + half)
        stack = np.stack([ds.GetRasterBand(bb).ReadAsArray().astype(np.float64)
                          for bb in range(b0, b1 + 1)])      # (窗口宽, rows, cols)
        wl = stack.shape[0]
        if wl >= 3 and order < wl:
            smoothed = savgol_filter(stack, wl, order, axis=0)
        else:
            smoothed = stack
        mid = b - b0
        arr = ds.GetRasterBand(b).ReadAsArray().astype(np.float64)
        valid = area_mask & (arr != NODATA_OUT)
        arr[valid] = smoothed[mid][valid]
        ds.GetRasterBand(b).WriteArray(arr.astype(np.float32))
        if b == 1 or b % 20 == 0 or b == n:
            logger.info("    S-G 平滑 %3d/%d 波段", b, n)
    ds.FlushCache()
    ds = None
    after = _scene_mean_spectrum(tif_path, area_mask)
    peaks_after = _diagnostic_peaks(wavelengths, after)
    for name in peaks_before:
        if name in peaks_after:
            _, d0 = peaks_before[name]
            _, d1 = peaks_after[name]
            drop = (d0 - d1) / max(abs(d0), 1e-9)
            logger.info("S-G 诊断峰 %s：深度 %.3f -> %.3f（相对变化 %+.1f%%）",
                        name, d0, d1, 100.0 * drop)
            if drop > SGO_PEAK_DEPTH_MAX_DROP:
                logger.warning("S-G 后 %s 峰深下降 %.1f%% 超过阈值 %.0f%%，建议减小 SGO_WINDOW",
                               name, 100.0 * drop, 100.0 * SGO_PEAK_DEPTH_MAX_DROP)
    return peaks_before, peaks_after


def build_area_mosaic(area, tiles, args, logger, wavelengths):
    grid = area["grid"]
    width, height = grid["width"], grid["height"]
    res = grid["res"]
    paths = mosaic_paths(area)
    os.makedirs(paths["dir"], exist_ok=True)

    priority_tiles = sort_tiles(tiles, args.priority)
    logger.info("  镶嵌: %d 个瓦片，优先级 %s", len(priority_tiles), args.priority)
    # v2.2.0：镶嵌前校验所有瓦片波段数与波长一致，禁止波段错配
    for tile in priority_tiles:
        if tile.get("raster_count") is not None and tile.get("raster_count") != len(wavelengths):
            raise ValueError("瓦片 %s 波段数 %s != 波长数 %d，禁止镶嵌（波段错配）"
                             % (tile.get("scene"), tile.get("raster_count"), len(wavelengths)))
    for index, tile in enumerate(priority_tiles, start=1):
        logger.info("    [%d] %s  %s  重叠 %.1f%%  云 %.1f%%",
                    index, tile["scene"][-16:], tile["date"],
                    tile["overlap"], tile["cloud_percent"])

    # v2.2.0：镶嵌前剔除跨景 NoData 严重的波段（按波长匹配），避免光谱掉点
    keep_idx = [i for i, w in enumerate(wavelengths)
                if not any(abs(float(w) - d) < 0.0015 for d in DROP_BAD_WL_UM)]
    if not keep_idx:
        raise ValueError("波段剔除后为空，无法镶嵌")
    good_wavelengths = [float(wavelengths[i]) for i in keep_idx]
    if len(keep_idx) != len(wavelengths):
        logger.info("  镶嵌前剔除 %d 个 NoData 严重波段，保留 %d 波段",
                    len(wavelengths) - len(keep_idx), len(keep_idx))

    aux = []
    for tile in priority_tiles:
        dataset = gdal.Open(tile["path"])
        band1 = dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
        dataset = None
        mask_dataset = gdal.Open(tile["mask_path"])
        mask = mask_dataset.GetRasterBand(1).ReadAsArray()
        mask_dataset = None
        if mask.shape != band1.shape:
            raise ValueError("瓦片与云掩膜尺寸不一致: %s" % tile["path"])
        valid = np.isfinite(band1) & (band1 != NODATA_OUT) & (mask != 255)
        aux.append({"tile": tile, "band1": band1, "mask": mask,
                    "valid": valid, "col_off": tile["col_off"],
                    "row_off": tile["row_off"],
                    "width": tile["width"], "height": tile["height"]})

    chosen = np.full((height, width), -1, dtype=np.int16)
    state = np.zeros((height, width), dtype=np.uint8)
    observations = np.zeros((height, width), dtype=np.uint16)
    clear_observations = np.zeros((height, width), dtype=np.uint16)

    # 统计每个像元的观测数 / 清晰观测数
    for item in aux:
        region = np.s_[item["row_off"]:item["row_off"] + item["height"],
                       item["col_off"]:item["col_off"] + item["width"]]
        observations[region] += item["valid"].astype(np.uint16)
        clear_observations[region] += (item["valid"] & (item["mask"] == 0)).astype(np.uint16)

    # 第一轮：清晰像元优先（不受场景优先级影响，任何清晰观测 > 有云观测）
    for index, item in enumerate(aux):
        region = np.s_[item["row_off"]:item["row_off"] + item["height"],
                       item["col_off"]:item["col_off"] + item["width"]]
        clear = item["valid"] & (item["mask"] == 0)
        chosen_region = chosen[region]
        new_pixels = clear & (chosen_region == -1)
        chosen_region[new_pixels] = index
        state_region = state[region]
        state_region[new_pixels] = 1

    # 第二轮：没有清晰观测时，用有云像元按优先级兜底
    for index, item in enumerate(aux):
        region = np.s_[item["row_off"]:item["row_off"] + item["height"],
                       item["col_off"]:item["col_off"] + item["width"]]
        cloudy = item["valid"] & (item["mask"] != 0)
        chosen_region = chosen[region]
        new_pixels = cloudy & (chosen_region == -1)
        if not np.any(new_pixels):
            continue
        chosen_region[new_pixels] = index
        state_region = state[region]
        mask_sub = item["mask"]
        state_region[new_pixels] = np.where(mask_sub[new_pixels] == 1, 2, 3)

    # 研究区多边形掩膜：覆盖率以“研究区面积内像元”为分母
    area_mask = rasterize_polygon(area["geom"], grid)
    coverage = mosaic_coverage(state, area_mask)

    # 逐波段组装（波段数 = 步骤2 剔除后的 N 波段，来自 wavelengths）
    driver = gdal.GetDriverByName("GTiff")
    tmp_reflectance = paths["reflectance"] + ".tmp.tif"
    if os.path.exists(tmp_reflectance):
        os.remove(tmp_reflectance)
    output = driver.Create(tmp_reflectance, width, height, len(keep_idx),
                           gdal.GDT_Float32,
                           options=build_gtiff_creation_options())
    output.SetGeoTransform((grid["x0"], res, 0.0, grid["y0"], 0.0, -res))
    target_srs = make_srs(TARGET_EPSG)
    output.SetProjection(target_srs.ExportToWkt())
    output.SetMetadataItem("AREA", area["key"])
    output.SetMetadataItem("AREA_LABEL", area["label"])
    output.SetMetadataItem("TILE_COUNT", str(len(priority_tiles)))
    output.SetMetadataItem("PRIORITY", args.priority)
    output.SetMetadataItem("CLOUD_MODE", args.cloud_mode)
    output.SetMetadataItem("TARGET_CRS", "EPSG:%d" % TARGET_EPSG)
    output.SetMetadataItem("TARGET_RES", "%g" % TARGET_RES)
    output.SetMetadataItem("COVERAGE_PERCENT", "%.3f" % coverage["coverage_percent"])
    output.SetMetadataItem("CLEAR_PERCENT", "%.3f" % coverage["clear_percent"])
    output.SetMetadataItem("CLOUD_FILLED_PERCENT",
                           "%.3f" % coverage["cloud_filled_percent"])
    output.SetMetadataItem("SOURCE_STEP2", args.in_root)
    output.SetMetadataItem("WAVELENGTHS_FILE", IN_WAVELENGTHS_NAME)

    sources = [gdal.Open(item["tile"]["path"]) for item in aux]
    total_bands = int(len(keep_idx))

    def _assemble_band(in_b):
        """v2.2.0：主景优先 + 逐波段兜底——先填主景（band1 选中景）清晰值保证光谱连贯，再逐波段补。"""
        clear_result = np.full((height, width), NODATA_OUT, dtype=np.float32)
        cloudy_result = np.full((height, width), NODATA_OUT, dtype=np.float32)
        # 主景优先：主景在该波段的清晰值先写入（即使优先级更低也不被覆盖）
        for index, item in enumerate(aux):
            region = np.s_[item["row_off"]:item["row_off"] + item["height"],
                           item["col_off"]:item["col_off"] + item["width"]]
            data = sources[index].GetRasterBand(in_b).ReadAsArray().astype(np.float32)
            valid_b = np.isfinite(data) & (data != NODATA_OUT)
            primary = (chosen[region] == index) & valid_b & (item["mask"] == 0)
            tgt = clear_result[region]
            tgt[primary] = data[primary]
        # 逐波段兜底：清晰优先、云兜底（只填仍 NoData 的）
        for index, item in enumerate(aux):
            region = np.s_[item["row_off"]:item["row_off"] + item["height"],
                           item["col_off"]:item["col_off"] + item["width"]]
            data = sources[index].GetRasterBand(in_b).ReadAsArray().astype(np.float32)
            valid_b = np.isfinite(data) & (data != NODATA_OUT)
            mask_sub = item["mask"]
            clear = valid_b & (mask_sub == 0)
            tgt_clear = clear_result[region]
            new_clear = clear & (tgt_clear == NODATA_OUT)
            tgt_clear[new_clear] = data[new_clear]
            cloudy = valid_b & (mask_sub != 0)
            tgt_cloudy = cloudy_result[region]
            new_cloudy = cloudy & (tgt_cloudy == NODATA_OUT)
            tgt_cloudy[new_cloudy] = data[new_cloudy]
        return np.where(clear_result != NODATA_OUT,
                        clear_result, cloudy_result).astype(np.float32)

    def _write_band(out_b, arr):
        band = output.GetRasterBand(out_b)
        band.SetNoDataValue(NODATA_OUT)
        band.WriteArray(arr)
        band.SetDescription("%.4f um" % good_wavelengths[out_b - 1])

    last_valid = None       # (out_b, arr) 最近一个无缺失的有效波段
    pending = []            # [(out_b, arr), ...] 待补值的连续缺失波段（缓冲<=3）
    gap_filled_px = 0
    interp_px = 0
    copy_px = 0

    def _flush_pending(next_idx, next_arr):
        """对缓冲的缺失波段做光谱补值并写出。

        next_idx 非空 = 找到了后续有效波段：缺口长度 <=3 用上下邻带线性插值，
        否则（连续缺口>3 / 末尾无后续）用最近有效波段拷贝兜底。
        """
        nonlocal last_valid, gap_filled_px, interp_px, copy_px
        cur = list(pending)      # 必须复制：pending.clear() 会清空原列表
        pending.clear()
        lv_idx, lv_arr = (last_valid if last_valid is not None else (None, None))
        for pb, p_arr in cur:
            miss = (p_arr == NODATA_OUT) & area_mask
            if np.any(miss):
                if lv_arr is not None and next_arr is not None:
                    gap_len = next_idx - lv_idx - 1
                    both = miss & (lv_arr != NODATA_OUT) & (next_arr != NODATA_OUT)
                    if gap_len <= 3 and next_idx > lv_idx:
                        frac = (pb - lv_idx) / float(next_idx - lv_idx)
                        p_arr[both] = lv_arr[both] + (next_arr[both] - lv_arr[both]) * frac
                        interp_px += int(both.sum())
                    else:
                        p_arr[both] = lv_arr[both]
                        copy_px += int(both.sum())
                    only_l = miss & ~both & (lv_arr != NODATA_OUT)
                    p_arr[only_l] = lv_arr[only_l]
                    copy_px += int(only_l.sum())
                    only_u = miss & (p_arr == NODATA_OUT) & (next_arr != NODATA_OUT)
                    p_arr[only_u] = next_arr[only_u]
                    copy_px += int(only_u.sum())
                elif lv_arr is not None:
                    ok = miss & (lv_arr != NODATA_OUT)
                    p_arr[ok] = lv_arr[ok]
                    copy_px += int(ok.sum())
                elif next_arr is not None:
                    ok = miss & (next_arr != NODATA_OUT)
                    p_arr[ok] = next_arr[ok]
                    copy_px += int(ok.sum())
                gap_filled_px += int(miss.sum())
            _write_band(pb, p_arr)
        if next_idx is not None:
            last_valid = (next_idx, next_arr)
        elif cur:
            last_valid = (cur[-1][0], cur[-1][1])

    # v2.2.0：逐波段独立选景 + 有用波段残留 NoData 光谱补值
    # （线性插值，连续缺口>3 用拷贝兜底；bbox 留白保持 NoData）
    area_px = max(int(area_mask.sum()), 1)
    for out_b, in_b0 in enumerate(keep_idx, start=1):
        arr = _assemble_band(int(in_b0) + 1)
        miss = (arr == NODATA_OUT) & area_mask
        miss_count = int(miss.sum())
        # 阈值过滤“多边形边界栅格化噪声”：只有缺口占比 >0.1% 才视为真实波段缺失并补值
        if miss_count > 0 and miss_count / float(area_px) > 0.001:
            pending.append((out_b, arr))
            if len(pending) > 3:
                _flush_pending(None, None)      # 连续缺口>3：拷贝兜底
        else:
            if pending:
                _flush_pending(out_b, arr)      # 用本带做上邻带插值
            _write_band(out_b, arr)
            last_valid = (out_b, arr)
        if out_b == 1 or out_b % 20 == 0 or out_b == total_bands:
            logger.info("    写出 %3d/%d 波段", out_b, total_bands)
    if pending:
        _flush_pending(None, None)              # 末尾剩余缺口：拷贝兜底
    for _i in range(len(sources)):
        sources[_i] = None
    output.FlushCache()
    output = None
    os.replace(tmp_reflectance, paths["reflectance"])

    # v2.2.0：S-G 光谱平滑（沿光谱维，保诊断吸收峰）+ 诊断峰核验
    sgo_peaks_before = None
    sgo_peaks_after = None
    if SGO_ENABLED:
        sgo_peaks_before, sgo_peaks_after = spectral_savgol_smooth(
            paths["reflectance"], good_wavelengths, area_mask, logger)

    # mask：0=无数据 1=清晰 2=云填充 3=云影填充
    _write_grid_mask(state, paths["mask"],
                     "mosaic mask (0=nodata 1=clear 2=cloud-filled 3=shadow-filled)",
                     grid, NODATA_OUT)
    # 来源景（1 基；0=无）
    source_index = np.where(chosen >= 0, chosen + 1, 0).astype(np.uint16)
    _write_grid_raster(source_index, paths["source"], gdal.GDT_UInt16,
                       grid, 0, "source tile index (1-based; 0=none)")
    observations_byte = np.minimum(observations, 255).astype(np.uint8)
    clear_byte = np.minimum(clear_observations, 255).astype(np.uint8)
    _write_grid_2band(observations_byte, clear_byte, paths["obs"], grid)

    info = {
        "area": area["key"],
        "area_label": area["label"],
        "source_shp": area["shp"],
        "source_epsg": area["source_epsg"],
        "target_crs_epsg": TARGET_EPSG,
        "target_res": TARGET_RES,
        "grid": grid,
        "band_count": total_bands,
        "dropped_bands_um": [float(w) for w in wavelengths if w not in good_wavelengths],
        "gap_filled_pixels": int(gap_filled_px),
        "gap_fill_interp_pixels": int(interp_px),
        "gap_fill_copy_pixels": int(copy_px),
        "savgol": {
            "enabled": bool(SGO_ENABLED),
            "window": int(SGO_WINDOW),
            "polyorder": int(SGO_POLYORDER),
            "peaks_before": {k: {"peak_nm": v[0], "depth": v[1]}
                             for k, v in (sgo_peaks_before or {}).items()},
            "peaks_after": {k: {"peak_nm": v[0], "depth": v[1]}
                            for k, v in (sgo_peaks_after or {}).items()},
        },
        "priority": args.priority,
        "cloud_mode": args.cloud_mode,
        "cloud_thresholds": dict(CLOUD_THRESHOLDS),
        "cloud_shadow_enabled": bool(CLOUD_SHADOW_ENABLE),
        "tile_count": len(priority_tiles),
        "date_range": sorted({item["tile"]["date"] for item in aux
                              if item["tile"].get("date")}),
        "tiles": [
            {
                "index": index + 1,
                "scene": item["tile"]["scene"],
                "date": item["tile"]["date"],
                "overlap_percent": item["tile"]["overlap"],
                "cloud_percent": item["tile"]["cloud_percent"],
                "clear_percent": item["tile"]["clear_percent"],
                "tile": item["tile"]["path"],
                "cloud_mask": item["tile"]["mask_path"],
            }
            for index, item in enumerate(aux)
        ],
        "coverage": coverage,
        "mosaic": {
            "reflectance": paths["reflectance"],
            "mask": paths["mask"],
            "source_scene": paths["source"],
            "observations": paths["obs"],
            "wavelengths": paths["wavelengths"],
        },
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    atomic_write_json(info, paths["info"])
    wavelengths_path = paths["wavelengths"]
    tmp_wavelengths = wavelengths_path + ".tmp.npy"
    np.save(tmp_wavelengths, np.asarray(good_wavelengths, dtype=np.float64))
    os.replace(tmp_wavelengths, wavelengths_path)
    logger.info("  mosaic 覆盖: 清晰 %.1f%%, 云填充 %.1f%%, 无数据 %.1f%%",
                coverage["clear_percent"], coverage["cloud_filled_percent"],
                coverage["nodata_percent"])
    return info


def mosaic_coverage(state, area_mask=None):
    """覆盖率统计：以研究区多边形内像元为分母（area_mask）。"""
    if area_mask is None:
        area_mask = np.ones_like(state, dtype=bool)
    area_pixels = int(np.sum(area_mask))
    covered = (state > 0) & area_mask
    clear = (state == 1) & area_mask
    cloud = (state == 2) & area_mask
    shadow = (state == 3) & area_mask
    covered_count = int(np.sum(covered))
    return {
        "total_pixels": area_pixels,
        "bbox_pixels": int(state.size),
        "nodata_pixels": area_pixels - covered_count,
        "clear_pixels": int(np.sum(clear)),
        "cloud_filled_pixels": int(np.sum(cloud)),
        "shadow_filled_pixels": int(np.sum(shadow)),
        "coverage_percent": 100.0 * covered_count / area_pixels if area_pixels else 0.0,
        "clear_percent": 100.0 * int(np.sum(clear)) / area_pixels if area_pixels else 0.0,
        "cloud_filled_percent": 100.0 * int(np.sum(cloud)) / area_pixels if area_pixels else 0.0,
        "shadow_filled_percent": 100.0 * int(np.sum(shadow)) / area_pixels if area_pixels else 0.0,
        "nodata_percent": 100.0 * (area_pixels - covered_count) / area_pixels if area_pixels else 0.0,
    }


def _write_grid_raster(array, path, dtype, grid, nodata, description):
    driver = gdal.GetDriverByName("GTiff")
    tmp_path = path + ".tmp.tif"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    dataset = driver.Create(tmp_path, grid["width"], grid["height"], 1,
                            dtype, options=build_byte_creation_options())
    dataset.SetGeoTransform((grid["x0"], grid["res"], 0.0,
                             grid["y0"], 0.0, -grid["res"]))
    dataset.SetProjection(make_srs(TARGET_EPSG).ExportToWkt())
    band = dataset.GetRasterBand(1)
    band.SetNoDataValue(nodata)
    band.WriteArray(array)
    band.SetDescription(description)
    dataset.FlushCache()
    dataset = None
    os.replace(tmp_path, path)


def _write_grid_mask(mask, path, description, grid, nodata):
    _write_grid_raster(mask.astype(np.uint8), path, gdal.GDT_Byte, grid,
                       0, description)


def _write_grid_2band(band1, band2, path, grid):
    driver = gdal.GetDriverByName("GTiff")
    tmp_path = path + ".tmp.tif"
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    dataset = driver.Create(tmp_path, grid["width"], grid["height"], 2,
                            gdal.GDT_Byte, options=build_byte_creation_options())
    dataset.SetGeoTransform((grid["x0"], grid["res"], 0.0,
                             grid["y0"], 0.0, -grid["res"]))
    dataset.SetProjection(make_srs(TARGET_EPSG).ExportToWkt())
    for index, data, description in ((1, band1, "observation count"),
                                     (2, band2, "clear observation count")):
        band = dataset.GetRasterBand(index)
        band.SetNoDataValue(255)
        band.WriteArray(data)
        band.SetDescription(description)
    dataset.FlushCache()
    dataset = None
    os.replace(tmp_path, path)


# ╔══════════════════════════════════════════════════════════════╗
# ║  质检成果图                                                    ║
# ╚══════════════════════════════════════════════════════════════╝

def _stretch(array, valid):
    values = array[valid]
    if values.size == 0:
        return np.zeros_like(array)
    low, high = np.percentile(values, [2, 98])
    if high <= low:
        high = low + 1e-6
    return np.clip((array - low) / (high - low), 0.0, 1.0)


def write_area_qc(area, info, logger):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import ListedColormap
    except Exception as exc:
        logger.warning("matplotlib 不可用，跳过质检图: %s", exc)
        return False

    paths = mosaic_paths(area)
    if not os.path.isfile(paths["reflectance"]):
        return False
    dataset = gdal.Open(paths["reflectance"])
    width, height = dataset.RasterXSize, dataset.RasterYSize
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))
    buf_w, buf_h = max(1, width // scale), max(1, height // scale)
    wavelengths = np.load(paths["wavelengths"])

    band_map = {
        "R": nearest_band(wavelengths, 0.66),
        "G": nearest_band(wavelengths, 0.55),
        "B": nearest_band(wavelengths, 0.47),
        "NIR": nearest_band(wavelengths, 0.86),
    }
    arrays = {}
    for name, index in band_map.items():
        arrays[name] = dataset.GetRasterBand(index + 1).ReadAsArray(
            buf_xsize=buf_w, buf_ysize=buf_h).astype(np.float64)
    geo_transform = dataset.GetGeoTransform()
    dataset = None

    mask_dataset = gdal.Open(paths["mask"])
    mask = mask_dataset.GetRasterBand(1).ReadAsArray(
        buf_xsize=buf_w, buf_ysize=buf_h)
    mask_dataset = None

    valid = np.ones((buf_h, buf_w), dtype=bool)
    for values in arrays.values():
        valid &= np.isfinite(values) & (values != NODATA_OUT)

    def combine(names):
        channels = [arrays[name] for name in names]
        rgb = np.ones(channels[0].shape + (3,), dtype=np.float64)
        for index, channel in enumerate(channels):
            rgb[..., index] = _stretch(channel, valid)
        rgb[~valid] = 1.0
        return np.clip(rgb, 0.0, 1.0)

    true_rgb = combine(("R", "G", "B"))
    false_rgb = combine(("NIR", "R", "G"))
    nir, red = arrays["NIR"], arrays["R"]
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = np.where((nir + red) > 0, (nir - red) / (nir + red), np.nan)
    ndvi = np.ma.masked_invalid(ndvi)

    extent = [geo_transform[0] / 1000.0,
              (geo_transform[0] + geo_transform[1] * width) / 1000.0,
              (geo_transform[3] + geo_transform[5] * height) / 1000.0,
              geo_transform[3] / 1000.0]

    figure, axes = plt.subplots(2, 2, figsize=(13, 10))
    ax_true, ax_false, ax_ndvi, ax_mask = axes.ravel()
    ax_true.imshow(true_rgb, extent=extent, origin="upper", interpolation="nearest")
    ax_true.set_title("True color reflectance", fontsize=10)
    ax_false.imshow(false_rgb, extent=extent, origin="upper", interpolation="nearest")
    ax_false.set_title("False color (NIR/R/G)", fontsize=10)
    ax_ndvi.imshow(ndvi, extent=extent, origin="upper", cmap="RdYlGn",
                   vmin=-0.2, vmax=0.9, interpolation="nearest")
    ax_ndvi.set_title("NDVI", fontsize=10)
    colour_map = ListedColormap(["#ffffff", "#2e7d32", "#e53935", "#fb8c00"])
    ax_mask.imshow(mask, extent=extent, origin="upper", cmap=colour_map,
                   vmin=0, vmax=3, interpolation="nearest")
    ax_mask.set_title("Mosaic mask: 0 nodata / 1 clear / 2 cloud / 3 shadow",
                      fontsize=9)

    # 研究区边界
    points = None
    try:
        ring = area["geom"].GetGeometryRef(0)
        points = ring.GetPoints() if ring is not None else None
    except Exception:
        points = None
    if points:
        xs = [point[0] / 1000.0 for point in points]
        ys = [point[1] / 1000.0 for point in points]
        for axis in (ax_true, ax_false, ax_ndvi, ax_mask):
            axis.plot(xs, ys, color="black", linewidth=1.0)

    for axis in (ax_true, ax_false, ax_ndvi, ax_mask):
        axis.set_xlabel("Easting (km)", fontsize=9)
        axis.set_ylabel("Northing (km)", fontsize=9)
        axis.tick_params(labelsize=8)

    coverage = info.get("coverage", {})
    dates = sorted({item["date"] for item in info.get("tiles", []) if item.get("date")})
    date_text = "%s~%s" % (dates[0], dates[-1]) if dates else "unknown"
    figure.suptitle(
        "%s | %d tiles | %s | coverage %.1f%% / clear %.1f%% / cloud-filled %.1f%%"
        % (area["label"], info.get("tile_count", 0), date_text,
           coverage.get("coverage_percent", 0.0),
           coverage.get("clear_percent", 0.0),
           coverage.get("cloud_filled_percent", 0.0)),
        fontsize=10)
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    figure.savefig(paths["qc"], dpi=QC_DPI, bbox_inches="tight")
    plt.close(figure)
    logger.info("  质检成果图: %s", paths["qc"])
    return True


# ╔══════════════════════════════════════════════════════════════╗
# ║  编目 / 状态表                                                 ║
# ╚══════════════════════════════════════════════════════════════╝

def expected_scene_names(directory):
    if not directory or not os.path.isdir(directory):
        return []
    names = set()
    for path in glob.glob(os.path.join(directory, "*.tar*")):
        base = os.path.basename(path)
        lower = base.lower()
        for ext in (".tar.gz", ".tgz", ".tar"):
            if lower.endswith(ext):
                base = base[:-len(ext)]
                break
        names.add(base)
    return sorted(names)


def _write_csv(path, header, rows):
    tmp_path = path + ".tmp"
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)
    os.replace(tmp_path, path)


def write_catalogs(out_root, areas, scenes, available_names, logger):
    catalog_dir = os.path.join(out_root, "_catalog")
    os.makedirs(catalog_dir, exist_ok=True)

    grid_payload = {}
    for area in areas:
        grid_payload[area["key"]] = {
            "label": area["label"],
            "shp": area["shp"],
            "source_epsg": area["source_epsg"],
            "target_crs_epsg": TARGET_EPSG,
            "target_res": TARGET_RES,
            "area_km2": area["area_km2"],
            "grid": area["grid"],
        }
    atomic_write_json(grid_payload, os.path.join(catalog_dir, "area_grids.json"))

    tile_rows = []
    coverage_rows = []
    for area in areas:
        tiles = scan_tiles(area)
        for tile in tiles:
            tile_rows.append([
                area["key"], tile["scene"], tile["date"],
                "%.3f" % tile["overlap"],
                tile["valid_pixels"], tile["total_pixels"],
                "%.2f" % tile["clear_percent"],
                "%.2f" % tile["cloud_percent"],
                "%.4f" % tile["aot550"] if tile["aot550"] is not None else "",
                tile["path"], tile["mask_path"],
            ])
        info = None
        paths = mosaic_paths(area)
        if os.path.isfile(paths["info"]):
            try:
                with open(paths["info"], "r", encoding="utf-8") as handle:
                    info = json.load(handle)
            except Exception:
                info = None
        coverage = (info or {}).get("coverage", {})
        dates = sorted({tile["date"] for tile in tiles if tile["date"]})
        coverage_rows.append([
            area["key"], area["label"], "%.1f" % area["area_km2"],
            len(tiles), len({tile["scene"] for tile in tiles}),
            "~".join([dates[0], dates[-1]]) if dates else "",
            "%.2f" % coverage.get("coverage_percent", 0.0),
            "%.2f" % coverage.get("clear_percent", 0.0),
            "%.2f" % coverage.get("cloud_filled_percent", 0.0),
            "%.2f" % coverage.get("nodata_percent", 0.0),
            paths["reflectance"] if os.path.isfile(paths["reflectance"]) else "",
        ])

    _write_csv(os.path.join(catalog_dir, "tiles_catalog.csv"),
               ["area", "scene", "date", "overlap_percent", "valid_pixels",
                "total_pixels", "clear_percent", "cloud_percent", "aot550",
                "tile", "cloud_mask"], tile_rows)
    _write_csv(os.path.join(catalog_dir, "area_coverage.csv"),
               ["area", "label", "area_km2", "tiles", "scenes", "dates",
                "coverage_percent", "clear_percent", "cloud_filled_percent",
                "nodata_percent", "mosaic"], coverage_rows)

    expected = expected_scene_names(EXPECTED_SCENES_DIR)
    processed = {scene["scene"]: scene for scene in scenes}
    status_rows = []
    for name in sorted(set(expected) | set(available_names)):
        scene = processed.get(name)
        date = scene["date"] if scene and scene.get("date") else ""
        if not date:
            parsed = parse_scene_date(name)
            date = parsed.strftime("%Y-%m-%d") if parsed else ""
        ready = name in available_names
        status_rows.append([name, "Y" if ready else "N", date,
                            "ready" if ready else "pending (步骤2 未完成)"])
    _write_csv(os.path.join(catalog_dir, "scene_status.csv"),
               ["scene", "step2_ready", "date", "note"], status_rows)
    logger.info("编目文件: %s", catalog_dir)


# ╔══════════════════════════════════════════════════════════════╗
# ║  分阶段处理                                                    ║
# ╚══════════════════════════════════════════════════════════════╝

def stage_tiles(areas, scenes, args, logger):
    summary = {}
    for area in areas:
        planned = done = skipped = failed = 0
        for scene in scenes:
            overlap = scene_area_overlap(scene, area, args.min_overlap)
            if overlap is None:
                continue
            planned += 1
            try:
                status = build_tile(scene, area, overlap, args, logger)
                if status == "skip":
                    skipped += 1
                else:
                    done += 1
            except Exception as exc:
                failed += 1
                logger.exception("瓦片失败: %s / %s", scene["scene"], area["key"])
                if not args.continue_on_error:
                    raise RuntimeError("瓦片失败即停（--continue-on-error 可关闭）: %s / %s"
                                       % (scene["scene"], area["key"])) from exc
        summary[area["key"]] = {"planned": planned, "done": done,
                                "skipped": skipped, "failed": failed}
        logger.info("区 %s 瓦片: 计划 %d，完成 %d，跳过 %d，失败 %d",
                    area["key"], planned, done, skipped, failed)
    return summary


def stage_cloud(areas, scenes, args, logger):
    summary = {}
    for area in areas:
        done = failed = 0
        for scene in scenes:
            overlap = scene_area_overlap(scene, area, args.min_overlap)
            if overlap is None:
                continue
            window = tile_window(area, overlap["geometry"])
            if window is None:
                continue
            tile_path = os.path.join(area["out_dir"], "tiles",
                                     scene["scene"] + TILE_SUFFIX)
            if not os.path.isfile(tile_path):
                continue
            mask_path = cloud_mask_path(tile_path)
            cutline_path = "/vsimem/step3_cloud_%s_%s.geojson" % (
                area["key"].replace(".", "_"), scene["scene"][-12:])
            intersection_wgs84 = transform_geometry(
                overlap["geometry"], TARGET_EPSG, 4326)
            write_geojson_cutline(intersection_wgs84, cutline_path)
            try:
                build_tile_cloud_mask(
                    scene, area, window, cutline_path, tile_path, mask_path,
                    args.cloud_mode, args.cloud_mask_dir, logger)
                done += 1
            except Exception as exc:
                failed += 1
                logger.exception("云掩膜重算失败: %s / %s",
                                 scene["scene"], area["key"])
                if not args.continue_on_error:
                    raise RuntimeError("云掩膜失败即停: %s / %s"
                                       % (scene["scene"], area["key"])) from exc
            finally:
                try:
                    gdal.Unlink(cutline_path)
                except Exception:
                    pass
        summary[area["key"]] = {"done": done, "failed": failed}
        logger.info("区 %s 云掩膜: 完成 %d，失败 %d", area["key"], done, failed)
    return summary


def resolve_mosaic_wavelengths(tiles, scenes, in_root, logger):
    lookup = {scene["scene"]: scene for scene in scenes}
    for tile in tiles:
        scene = lookup.get(tile["scene"])
        if scene is not None:
            return scene["wavelengths"]
        path = os.path.join(in_root, tile["scene"], IN_WAVELENGTHS_NAME)
        if os.path.isfile(path):
            return np.load(path, allow_pickle=False)
    raise FileNotFoundError("无法找到任何瓦片对应的 wavelengths_clean.npy")


def stage_mosaic(areas, scenes, args, logger):
    summary = {}
    for area in areas:
        tiles = scan_tiles(area)
        if not tiles:
            logger.warning("区 %s 没有可用瓦片，跳过镶嵌", area["key"])
            summary[area["key"]] = {"tiles": 0, "coverage": 0.0,
                                     "failed": False}
            continue
        wavelengths = resolve_mosaic_wavelengths(tiles, scenes,
                                                 args.in_root, logger)
        if args.stage == "all" and not args.overwrite \
                and mosaic_is_current(area, tiles, wavelengths):
            logger.info("区 %s mosaic 已是最新，跳过（--stage mosaic 可强制重拼）",
                        area["key"])
            summary[area["key"]] = {"tiles": len(tiles), "coverage": None,
                                     "failed": False}
            continue
        info = build_area_mosaic(area, tiles, args, logger, wavelengths)
        if MAKE_QC and not args.no_qc:
            qc_ok = write_area_qc(area, info, logger)
            if not qc_ok:
                raise RuntimeError("质检图生成失败（用户必需成果，不允许静默跳过）: %s" % area["key"])
        coverage = info["coverage"]["coverage_percent"]
        failed = bool(args.require_full_coverage
                      and coverage < MIN_FULL_COVERAGE)
        if failed:
            logger.error("区 %s 覆盖率 %.1f%% 低于要求 %.1f%%",
                         area["key"], coverage, MIN_FULL_COVERAGE)
        summary[area["key"]] = {"tiles": len(tiles), "coverage": coverage,
                                 "failed": failed}
    return summary


# ╔══════════════════════════════════════════════════════════════╗
# ║  批处理主流程                                                  ║
# ╚══════════════════════════════════════════════════════════════╝

def main(argv=None):
    global TARGET_EPSG, TARGET_RES
    args = parse_args(argv)
    out_root = os.path.abspath(args.out_root)
    args.out_root = out_root
    os.makedirs(out_root, exist_ok=True)
    log_path = args.log or os.path.join(out_root, "step3_batch.log")
    logger = setup_logging(log_path, args.verbose)
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)

    in_root = os.path.abspath(args.in_root)
    args.in_root = in_root
    try:
        scene_dirs = discover_scenes(in_root)
    except Exception as exc:
        logger.error("扫描步骤2 成果失败: %s", exc)
        return 2
    available_names = {os.path.basename(os.path.normpath(path))
                       for path in scene_dirs}
    if args.scene:
        scene_dirs = [path for path in scene_dirs
                      if _match_any(os.path.basename(path), args.scene)]
    if args.limit > 0:
        scene_dirs = scene_dirs[: args.limit]

    scenes = []
    for scene_dir in scene_dirs:
        try:
            scenes.append(read_scene_info(scene_dir))
        except Exception as exc:
            logger.exception("读取场景失败: %s", scene_dir)
            if not args.continue_on_error:
                raise RuntimeError("读取场景失败即停: %s" % scene_dir) from exc
    if not scenes:
        logger.error("在 %s 下没有可用的步骤2 成果", in_root)
        return 2

    # v2.1.0 坐标单一真源：CLI > 步骤2 georef > 代码默认。多景不一致即停，避免混合网格。
    georef_epsgs = {s["target_epsg"] for s in scenes if s.get("target_epsg") is not None}
    georef_reses = {s["target_res"] for s in scenes if s.get("target_res") is not None}
    if len(georef_epsgs) > 1:
        raise ValueError("步骤2 各景 georef.target_crs_epsg 不一致: %s，"
                         "混合网格不可镶嵌，请统一后重跑步骤1/2" % sorted(georef_epsgs))
    if len(georef_reses) > 1:
        raise ValueError("步骤2 各景 georef.target_res 不一致: %s，"
                         "混合分辨率不可镶嵌" % sorted(georef_reses))
    if args.epsg is not None:
        TARGET_EPSG = int(args.epsg)                     # CLI 覆盖优先
    elif georef_epsgs:
        TARGET_EPSG = next(iter(georef_epsgs))           # 从 georef 单一真源读取
    if args.res is not None:
        TARGET_RES = float(args.res)                     # CLI 覆盖优先
    elif georef_reses:
        TARGET_RES = next(iter(georef_reses))
    logger.info("目标网格真源: EPSG:%d / %g m（CLI epsg=%s res=%s，"
                "georef epsg=%s res=%s）", TARGET_EPSG, TARGET_RES,
                args.epsg, args.res,
                sorted(georef_epsgs) if georef_epsgs else None,
                sorted(georef_reses) if georef_reses else None)

    try:
        areas = load_areas(args.area_root, TARGET_EPSG, args.buffer_m,
                           args.area, logger)
    except Exception as exc:
        logger.error("加载研究区失败: %s", exc)
        return 2
    for area in areas:
        area["out_dir"] = os.path.join(out_root, area["key"])

    expected = expected_scene_names(EXPECTED_SCENES_DIR)
    pending = sorted(set(expected) - available_names)
    logger.info("=" * 72)
    logger.info("ZY1E AHSI 步骤3（方案A）：RPC 正射 + 研究区裁剪 + 云感知镶嵌")
    logger.info("步骤2 成果: %s", in_root)
    logger.info("输出根目录: %s", out_root)
    logger.info("研究区    : %d 个（%s）", len(areas),
                ", ".join(area["key"] for area in areas))
    logger.info("可用景    : %d 景", len(scenes))
    if expected:
        logger.info("待补景    : %d / %d（步骤2 未完成）", len(pending), len(expected))
    logger.info("目标网格  : EPSG:%d / %g m", TARGET_EPSG, TARGET_RES)
    logger.info("阶段      : %s | 优先级: %s | 云模式: %s",
                args.stage, args.priority, args.cloud_mode)
    logger.info("日志文件  : %s", log_path)
    logger.info("=" * 72)

    if args.dry_run:
        total_pairs = 0
        for area in areas:
            pairs = []
            for scene in scenes:
                overlap = scene_area_overlap(scene, area, args.min_overlap)
                if overlap is None:
                    continue
                window = tile_window(area, overlap["geometry"])
                if window is None:
                    continue
                pairs.append((scene, overlap, window))
            total_pairs += len(pairs)
            logger.info("[%s] 瓦片 %d 个", area["key"], len(pairs))
            for scene, overlap, window in pairs:
                logger.info("    %s  %s  重叠 %.1f%%  瓦片 %d×%d",
                            scene["scene"][-16:], scene["date"],
                            overlap["percent"], window["width"], window["height"])
        logger.info("dry-run: 共计划 %d 个瓦片；未执行任何处理", total_pairs)
        return 0

    started = time.time()
    failures = []
    tile_summary = {}
    cloud_summary = {}
    if args.stage in ("all", "tiles"):
        tile_summary = stage_tiles(areas, scenes, args, logger)
        for area_key, item in tile_summary.items():
            if item.get("failed"):
                failures.append(area_key)
    if args.stage == "cloud":
        cloud_summary = stage_cloud(areas, scenes, args, logger)
        for area_key, item in cloud_summary.items():
            if item.get("failed"):
                failures.append(area_key)
    if args.stage in ("all", "mosaic"):
        mosaic_summary = stage_mosaic(areas, scenes, args, logger)
        for area_key, item in mosaic_summary.items():
            if item.get("failed"):
                failures.append(area_key)

    write_catalogs(out_root, areas, scenes, available_names, logger)

    logger.info("=" * 72)
    logger.info("批处理结束，总用时 %s", fmt_elapsed(time.time() - started))
    for area in areas:
        paths = mosaic_paths(area)
        if os.path.isfile(paths["reflectance"]):
            info = None
            try:
                with open(paths["info"], "r", encoding="utf-8") as handle:
                    info = json.load(handle)
            except Exception:
                info = None
            coverage = (info or {}).get("coverage", {})
            logger.info("  [%s] 覆盖 %.1f%% | 清晰 %.1f%% | 云填充 %.1f%% | %s",
                        area["key"], coverage.get("coverage_percent", 0.0),
                        coverage.get("clear_percent", 0.0),
                        coverage.get("cloud_filled_percent", 0.0),
                        human_size(os.path.getsize(paths["reflectance"])))
        else:
            logger.info("  [%s] 尚未生成 mosaic", area["key"])
    if failures:
        logger.error("有 %d 个区未达到覆盖要求: %s",
                     len(failures), ", ".join(failures))
        return 1
    logger.info("全部完成，成果位于: %s", out_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
