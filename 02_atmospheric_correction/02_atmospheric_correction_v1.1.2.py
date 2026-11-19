# -*- coding: utf-8 -*-
"""
02_atmospheric_correction.py
ZY1E AHSI 步骤2：在步骤1 正射网格上使用 6S 将辐亮度换算为地表反射率（v1.1.2）。
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
import tarfile
import time
import zipfile

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit(
        "当前 Python 解释器缺少 numpy：%s\n"
        "请用装好 GDAL + Py6S 的解释器运行本脚本，例如：\n"
        "    E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe 02_atmospheric_correction_v1.1.2.py"
        % exc)

try:
    from osgeo import gdal, osr
except ImportError as exc:
    raise SystemExit(
        "当前 Python 解释器缺少 GDAL（osgeo）：%s\n"
        "请用装好 GDAL + Py6S 的解释器运行本脚本，例如：\n"
        "    E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe 02_atmospheric_correction_v1.1.2.py"
        % exc)

gdal.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")

# =================================================================
# 配置区（按需修改）
# =================================================================

# 输入步骤1 正射成果目录；输出大气校正成果根目录
IN_ROOT_STEP1 = r"E:\ZhouYang\test\ZY1E_output_step1"
OUT_ROOT_STEP2 = r"E:\ZhouYang\test\ZY1E_output_step2"

# 原始 TAR 目录；读取观测几何和日期，优先使用 georef.json 的 source_tar
TAR_DIR = r"E:\ZhouYang\test\ZY1E_AHSI"

# MODIS MOD04/MYD04 AOD 产品目录（按日期及场景足迹选有效像元）
MODIS_DIR = r"E:\ZhouYang\test\MOD04_L2"
# ASTGTM DEM 目录；可放置 GeoTIFF/IMG 或包含栅格的 ZIP
DEM_DIR = r"E:\ZhouYang\test\DEM"

# 6S 可执行程序路径；留空时依次查环境变量、脚本 _bin 目录和 PATH
SIXS_EXE = ""

# 自动按日期和纬度选择大气模式；关闭时使用手工模式
AUTO_DETECT_ATMOS = True
ATMOS_PROFILE_MANUAL = "MidlatitudeWinter"  # Tropical/MidlatitudeSummer/MidlatitudeWinter/SubarcticSummer/SubarcticWinter/USStandard

# 气溶胶类型（无土地覆盖数据时默认 Continental）
AERO_PROFILE_MANUAL = "Continental"  # Continental/Urban/Maritime/Desert
AUTO_DETECT_AERO = False             # 启用时按经纬度粗略判断

# AOT 来源：auto=本地 MODIS→AERONET→手工；modis/aeronet/manual 可单独指定
AOT_SOURCE = "auto"
AOT550_MANUAL = 0.2                  # 自动来源无有效值时使用的 AOT550
AOT_MIN_VALID = 0.0                  # MODIS 有效 AOD 开区间下限
AOT_MAX_VALID = 5.0                  # MODIS 有效 AOD 闭区间上限
AERONET_ENABLED = True               # 是否允许联网查询 AERONET
AERONET_MAX_DISTANCE_KM = 300        # AERONET 站点最远距离（km）
AERONET_SEARCH_DAYS = 3              # AERONET 日期搜索范围（天）

# 自动从足迹内 DEM 计算平均海拔；失败时使用固定兜底值
AUTO_DETECT_DEM = True
DEM_RES = 90.0                       # DEM 计算网格分辨率（m）
ELEVATION_KM_FALLBACK = 0.2          # DEM 不可用时的海拔（km）

# 水汽为空时使用 6S 默认廓线；指定水汽后可同时指定臭氧
WATER_VAPOR_CM = None               # 水汽含量（g/cm²），None 表示使用大气模式廓线
OZONE_CM = None                     # 臭氧含量（cm-atm），None 时使用下方默认值
OZONE_CM_DEFAULT = 0.3              # 手工水汽廓线下的默认臭氧量（cm-atm）

# 反射率裁剪范围；None 表示不裁剪
REFLECTANCE_CLIP = (0.0, 1.0)
MASK_INVALID_RADIANCE = True         # 是否掩膜低于阈值的输入辐亮度
RADIANCE_MIN_VALID = 0.0             # 输入辐亮度有效下限
NODATA_OUT = -9999.0                 # 输出无效像元值

# XML 缺少观测几何时使用的手工兜底值（角度）
SOLAR_ZENITH_MANUAL = None
SOLAR_AZIMUTH_MANUAL = None
VIEW_ZENITH_MANUAL = None
VIEW_AZIMUTH_MANUAL = None

# 批处理和 GDAL 参数
OVERWRITE = False                    # 是否重做已有成果
GDAL_CACHE_MB = 1024                 # GDAL 缓存大小（MB）

# GeoTIFF 写出设置
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]
GEOTIFF_COMPRESS = "DEFLATE"         # DEFLATE/LZW/None
GEOTIFF_COMPRESS_LEVEL = 6           # 压缩等级

# 质检图设置
MAKE_QC = True                       # 是否生成 preview_quality.png（需要 matplotlib）
QC_GEOREF_PREVIEW = True             # 是否生成带地理参考的降采样预览
QC_PREVIEW_SCALE = 2.0               # 预览像元大小相对输入的倍数
QC_MAX_SIZE = 1400                   # 普通预览图长边最大像素数
QC_DPI = 150                         # 质检图分辨率

# 步骤1 输入文件名：新命名优先，同时兼容旧版固定名称
IN_RADIANCE_SUFFIX = "_cal_badband_ortho.tif"
LEGACY_RADIANCE_NAME = "radiance_ortho.tif"
IN_WAVELENGTHS_NAME = "wavelengths_clean.npy"
IN_GEOREF_NAME = "georef.json"

# 本步输出文件名：保留场景名，并以缩写标注正射与大气校正处理
OUT_REFLECTANCE_SUFFIX = "_cal_badband_ortho_atmcorr.tif"
LEGACY_REFLECTANCE_NAMES = (
    "_ortho_atmcorr.tif",
    "reflectance_atm.tif",
)
OUT_WAVELENGTHS_NAME = "wavelengths_clean.npy"  # 随输出保存波长，供步骤3 使用
OUT_GEOREF_NAME = "georef.json"                # 复制步骤1 元数据用于追溯
OUT_COEF_NAME = "py6s_coefficients.npy"         # 每个波段的 [xa, xb, xc] 系数
OUT_STATS_NAME = "band_stats_reflectance.csv"   # 每波段反射率统计
OUT_PARAMS_JSON = "used_parameters.json"        # 本次参数和数据来源
OUT_PARAMS_TXT = "used_parameters.txt"          # 参数摘要文本
OUT_PREVIEW_GEO_NAME = "preview_rgb_geo.tif"    # 带坐标的降采样预览
OUT_QC_NAME = "preview_quality.png"             # 质检图

LOG_NAME = "zy1e.step2"                         # 日志记录器名称

# 用户可选的 Py6S 大气和气溶胶模式名称
AERO_PROFILE_CHOICES = ("Continental", "Urban", "Maritime", "Desert")
ATMOS_PROFILE_CHOICES = ("Tropical", "MidlatitudeSummer", "MidlatitudeWinter",
                         "SubarcticSummer", "SubarcticWinter", "USStandard")


# =================================================================
# 通用工具
# =================================================================

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="ZY1E AHSI 步骤2：6S 大气校正（正射格网，直接承接步骤1 正射辐亮度）")
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP1,
                        help="步骤1 成果根目录（默认: %(default)s）")
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP2,
                        help="本步成果根目录，每景生成同名子文件夹（默认: %(default)s）")
    parser.add_argument("--scene", action="append", default=None,
                        help="只处理匹配的场景名，可多次指定，支持通配符（如 \"*E117.20*\"）")
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE,
                        help="已有成果时覆盖重做（默认跳过）")
    parser.add_argument("--limit", type=int, default=0,
                        help="只处理前 N 景，便于调试（默认 0 = 全部）")
    parser.add_argument("--dry-run", action="store_true",
                        help="只列出待处理场景与目标目录，不执行")
    parser.add_argument("--no-qc", action="store_true",
                        help="跳过质检图生成（保留统计 CSV 与预览 tif 开关不受影响）")
    parser.add_argument("--log", default=None,
                        help="日志文件路径（默认 <out>/step2_batch.log）")
    parser.add_argument("--verbose", action="store_true", help="输出调试日志")

    # 手工覆盖 6S 关键参数（便于快速试验，不写死文件）
    parser.add_argument("--aot", type=float, default=None,
                        help="直接指定 AOT550（覆盖自动 AOD 检索）")
    parser.add_argument("--aot-source", choices=("auto", "modis", "aeronet", "manual"),
                        default=None,
                        help="覆盖配置区 AOT_SOURCE（默认 auto=MODIS→AERONET→手工）")
    parser.add_argument("--elev", type=float, default=None,
                        help="直接指定目标海拔（km，覆盖 DEM 自动计算）")
    parser.add_argument("--atmos", default=None, choices=ATMOS_PROFILE_CHOICES,
                        help="直接指定 6S 大气模式")
    parser.add_argument("--aero", default=None, choices=AERO_PROFILE_CHOICES,
                        help="直接指定 6S 气溶胶类型")
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


def discover_scene_dirs(in_root):
    """返回步骤1 成果子目录列表；若传入的就是单景目录则直接返回它。"""
    if not os.path.isdir(in_root):
        raise NotADirectoryError("步骤1 成果目录不存在: %s" % in_root)
    if resolve_radiance_path(in_root):
        return [in_root]
    scene_dirs = sorted(glob.glob(os.path.join(in_root, "*")))
    return [path for path in scene_dirs
            if os.path.isdir(path) and resolve_radiance_path(path)]


def resolve_radiance_path(scene_dir):
    renamed = sorted(glob.glob(os.path.join(
        scene_dir, "*" + IN_RADIANCE_SUFFIX)))
    if renamed:
        return renamed[0]
    legacy = os.path.join(scene_dir, LEGACY_RADIANCE_NAME)
    return legacy if os.path.isfile(legacy) else None


def fmt_elapsed(seconds):
    seconds = max(0.0, float(seconds))
    return "%d:%04.1f" % (int(seconds // 60), seconds % 60.0)


def log_stage_start(logger, number, total, description):
    logger.info("[步骤 %d/%d] 开始：%s", number, total, description)
    return time.time()


def log_stage_done(logger, number, total, description, started):
    logger.info("[步骤 %d/%d] 完成：%s，用时 %s",
                number, total, description, fmt_elapsed(time.time() - started))


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


def copy_file_atomic(src, dst):
    tmp_path = dst + ".tmp"
    shutil.copyfile(src, tmp_path)
    os.replace(tmp_path, dst)


def load_georef_json(path):
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


# =================================================================
# 产品元数据解析
# =================================================================

def locate_tar(scene, georef, tar_dir, logger):
    """定位原始 tar 压缩包：先按 georef.json 记录，再按场景名匹配。"""
    candidates = []
    if georef and georef.get("source_tar"):
        candidates.append(os.path.join(tar_dir, georef["source_tar"]))
        candidates.append(georef["source_tar"])
    for ext in (".tar", ".tar.gz", ".tgz"):
        candidates.append(os.path.join(tar_dir, scene + ext))
    for path in candidates:
        if path and os.path.isfile(path):
            return os.path.abspath(path)

    # 最后兜底：按场景名在 TAR_DIR 里模糊匹配
    hits = []
    for ext in (".tar", ".tar.gz", ".tgz"):
        hits.extend(glob.glob(os.path.join(tar_dir, "*" + scene + "*" + ext)))
    hits = sorted(set(os.path.abspath(path) for path in hits))
    if hits:
        if len(hits) > 1:
            logger.warning("场景 %s 匹配到 %d 个压缩包，使用第一个: %s",
                           scene, len(hits), os.path.basename(hits[0]))
        return hits[0]
    return None


def decode_text(raw):
    text = raw.decode("utf-8", errors="ignore")
    return text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")


def read_product_xml(tar_path, logger):
    """从 tar 中读取产品元数据 XML 文本。"""
    with tarfile.open(tar_path, "r:*") as tar:
        hits = [name for name in tar.getnames()
                if name.lower().endswith(".xml")
                and not name.endswith("/")]
        if not hits:
            raise FileNotFoundError("压缩包中没有 XML 元数据: %s" % tar_path)
        # 优先包含 ProductMetaData 的产品 XML
        preferred = None
        for name in sorted(hits):
            handle = tar.extractfile(name)
            if handle is None:
                continue
            with handle:
                text = decode_text(handle.read())
            if "ProductMetaData" in text or "SolarZenith" in text:
                preferred = (name, text)
                break
        if preferred is None:
            name = sorted(hits)[0]
            handle = tar.extractfile(name)
            if handle is None:
                raise FileNotFoundError("无法读取 XML: %s" % name)
            with handle:
                preferred = (name, decode_text(handle.read()))
        if len(hits) > 1:
            logger.debug("压缩包内 XML %d 个，使用: %s", len(hits), preferred[0])
        return preferred[1]


def _find_float(text, tag):
    match = re.search(r"<%s>\s*([-+0-9.eE]+)\s*</%s>" % (tag, tag),
                      text, re.IGNORECASE)
    if not match:
        return None
    try:
        return float(match.group(1))
    except ValueError:
        return None


def parse_product_metadata(xml_text, scene, logger):
    """解析产品 XML：观测几何、日期、中心位置、四角足迹。"""
    meta = {
        "satellite": (_find_text(xml_text, "SatelliteID") or "ZY1E"),
        "sensor": (_find_text(xml_text, "SensorID") or "AHSI"),
        "solar_z": _find_float(xml_text, "SolarZenith"),
        "solar_a": _find_float(xml_text, "SolarAzimuth"),
        "view_z": _find_float(xml_text, "SatelliteZenith"),
        "view_a": _find_float(xml_text, "SatelliteAzimuth"),
        "center_lat": _find_float(xml_text, "CenterLatitude"),
        "center_lon": _find_float(xml_text, "CenterLongitude"),
        "date": None,
        "footprint": {},
    }

    match = re.search(r"<StartTime>\s*(\d{4})-(\d{2})-(\d{2})"
                      r"(?:[ T](\d{2}):(\d{2}):(\d{2}))?", xml_text)
    if match:
        parts = [int(value) for value in match.groups()[:6] if value is not None]
        if len(parts) == 3:
            parts += [0, 0, 0]
        try:
            meta["date"] = datetime.datetime(*parts)
        except ValueError:
            logger.warning("StartTime 解析异常: %s", match.group(0))

    corner_map = {"TL": "TopLeft", "TR": "TopRight",
                  "BR": "BottomRight", "BL": "BottomLeft"}
    for abbr, corner in corner_map.items():
        lat = _find_float(xml_text, corner + "Latitude")
        lon = _find_float(xml_text, corner + "Longitude")
        if lat is not None and lon is not None:
            meta["footprint"][abbr] = (lon, lat)

    missing = [key for key in ("solar_z", "solar_a", "view_z", "view_a")
               if meta[key] is None]
    if missing:
        logger.warning("产品 XML 缺少观测几何: %s", ", ".join(missing))
    return meta


def _find_text(text, tag):
    match = re.search(r"<%s>\s*([^<]+?)\s*</%s>" % (tag, tag), text, re.IGNORECASE)
    return match.group(1).strip() if match else None


def parse_scene_name(scene):
    """从场景名兜底解析日期/中心经纬度。"""
    info = {"date": None, "center_lat": None, "center_lon": None}
    match = re.search(r"_(\d{8})_", scene)
    if match:
        try:
            info["date"] = datetime.datetime.strptime(match.group(1), "%Y%m%d")
        except ValueError:
            pass
    match = re.search(r"_([EW])(\d{1,3}\.\d+)_([NS])(\d{1,3}\.\d+)", scene)
    if match:
        lon = float(match.group(2)) * (1.0 if match.group(1) == "E" else -1.0)
        lat = float(match.group(4)) * (1.0 if match.group(3) == "N" else -1.0)
        info["center_lon"] = lon
        info["center_lat"] = lat
    return info


def resolve_metadata(scene, scene_dir, georef, tar_dir, logger):
    """汇总本景 6S 计算所需的元数据：观测几何、日期、中心位置、足迹、tar 路径。"""
    name_info = parse_scene_name(scene)
    meta = {
        "scene": scene,
        "tar": None,
        "solar_z": None, "solar_a": None, "view_z": None, "view_a": None,
        "date": name_info["date"],
        "center_lat": name_info["center_lat"],
        "center_lon": name_info["center_lon"],
        "footprint": {},
        "satellite": "ZY1E", "sensor": "AHSI",
    }

    tar_path = locate_tar(scene, georef, tar_dir, logger)
    if tar_path:
        meta["tar"] = tar_path
        try:
            xml_text = read_product_xml(tar_path, logger)
            xml_meta = parse_product_metadata(xml_text, scene, logger)
            for key in ("solar_z", "solar_a", "view_z", "view_a",
                        "satellite", "sensor"):
                if xml_meta.get(key) is not None:
                    meta[key] = xml_meta[key]
            if xml_meta.get("date") is not None:
                meta["date"] = xml_meta["date"]
            for key in ("center_lat", "center_lon"):
                if xml_meta.get(key) is not None:
                    meta[key] = xml_meta[key]
            if xml_meta.get("footprint"):
                meta["footprint"] = dict(xml_meta["footprint"])
        except Exception as exc:
            logger.warning("读取产品 XML 失败，改用 georef.json/场景名兜底: %s", exc)
    else:
        logger.warning("未找到原始压缩包，观测几何将使用手工兜底值")

    # 足迹优先用步骤1 的 georef.json 的 footprint_wgs84（WGS84，与步骤1/步骤3 一致）
    if georef and georef.get("footprint_wgs84"):
        footprint = {}
        for key, value in georef["footprint_wgs84"].items():
            footprint[key] = (float(value[0]), float(value[1]))
        meta["footprint"] = footprint

    # 手工几何兜底
    manual = {"solar_z": SOLAR_ZENITH_MANUAL, "solar_a": SOLAR_AZIMUTH_MANUAL,
              "view_z": VIEW_ZENITH_MANUAL, "view_a": VIEW_AZIMUTH_MANUAL}
    for key, value in manual.items():
        if meta[key] is None and value is not None:
            meta[key] = float(value)
            logger.warning("观测几何 %s 使用手工兜底值 %s", key, value)

    missing = [key for key in ("solar_z", "solar_a", "view_z", "view_a")
               if meta[key] is None]
    if missing:
        raise ValueError(
            "缺少观测几何参数 %s；请在脚本配置区填写 *_MANUAL 兜底值，"
            "或确认 TAR_DIR 指向原始压缩包" % ", ".join(missing))

    if not meta["footprint"]:
        raise ValueError(
            "缺少四角足迹（georef.json 与产品 XML 均未提供）；"
            "请确认步骤1 正常完成（georef.json）")

    if meta["center_lat"] is None or meta["center_lon"] is None:
        lons = [value[0] for value in meta["footprint"].values()]
        lats = [value[1] for value in meta["footprint"].values()]
        meta["center_lon"] = float(np.mean(lons))
        meta["center_lat"] = float(np.mean(lats))
        logger.warning("中心经纬度缺失，用足迹均值兜底: (%.4f, %.4f)",
                       meta["center_lat"], meta["center_lon"])

    if meta["date"] is None:
        meta["date"] = datetime.datetime(2000, 1, 1)
        logger.warning("成像日期缺失，6S 月份/日期使用 2000-01-01 兜底")

    return meta


def footprint_ring(footprint):
    """返回闭合的 WGS84 足迹环 [[lon, lat], ...]。"""
    order = [key for key in ("TL", "TR", "BR", "BL") if key in footprint]
    if len(order) < 3:
        order = list(footprint.keys())
    ring = [[float(footprint[key][0]), float(footprint[key][1])] for key in order]
    if ring:
        ring.append(list(ring[0]))
    return ring


def points_in_polygon(x, y, ring):
    """射线法判断点是否在闭合多边形内（ring 首尾可重复）。"""
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    inside = np.zeros(x.shape, dtype=bool)
    points = np.asarray(ring, dtype=np.float64)
    if points.shape[0] > 1 and np.allclose(points[0], points[-1]):
        points = points[:-1]
    count = points.shape[0]
    for idx in range(count):
        x1, y1 = points[idx]
        x2, y2 = points[(idx + 1) % count]
        cond = (y1 > y) != (y2 > y)
        if not np.any(cond):
            continue
        with np.errstate(divide="ignore", invalid="ignore"):
            xin = (x2 - x1) * (y - y1) / (y2 - y1) + x1
        inside ^= cond & (x < xin)
    return inside


# =================================================================
# AOT550：MODIS MOD04_L2 / AERONET / 手工
# =================================================================

def find_mod04_files(modis_dir, date):
    """按成像日期（UTC 年积日）递归查找本地 MOD04/MYD04 AOD 产品。"""
    if not modis_dir or not os.path.isdir(modis_dir):
        return []
    doy = date.strftime("%Y%j")
    patterns = ["MOD04_L2.A%s*" % doy, "MYD04_L2.A%s*" % doy,
                "MOD04_3K.A%s*" % doy, "MYD04_3K.A%s*" % doy]
    hits = set()
    for pattern in patterns:
        hits.update(glob.glob(os.path.join(modis_dir, "**", pattern),
                              recursive=True))
    return sorted(path for path in hits
                  if os.path.isfile(path) and path.lower().endswith((".hdf", ".hdfeos")))


def read_mod04_aod(path, ring, logger):
    """读取一个 MOD04_L2 文件，返回足迹多边形内有效的 AOD550 列表。"""
    logger.info("  检查 MODIS 文件：%s", os.path.basename(path))
    try:
        from pyhdf.SD import SD, SDC
    except ImportError:
        logger.warning("未安装 pyhdf，无法读取 MODIS MOD04 AOD: %s", path)
        return None

    sds_priority = ("AOD_550_Dark_Target_Deep_Blue_Combined",
                    "Deep_Blue_Aerosol_Optical_Depth_550_Land",
                    "Optical_Depth_Land_And_Ocean")
    handle = None
    try:
        handle = SD(path, SDC.READ)
        lat = handle.select("Latitude")[:].astype(np.float64)
        lon = handle.select("Longitude")[:].astype(np.float64)
        lons = np.asarray([point[0] for point in ring])
        lats = np.asarray([point[1] for point in ring])
        bbox_mask = ((lat >= lats.min() - 0.05) & (lat <= lats.max() + 0.05)
                     & (lon >= lons.min() - 0.05) & (lon <= lons.max() + 0.05))
        if not np.any(bbox_mask):
            logger.warning("  MOD04 %s 的经纬度范围不覆盖本景足迹",
                           os.path.basename(path))
            return None
        inside = bbox_mask & points_in_polygon(lon, lat, ring)
        if not np.any(inside):
            logger.warning("  MOD04 %s 与本景外接范围有交集，但没有像元中心落入本景足迹",
                           os.path.basename(path))
            return None

        # 按产品优先级读取 AOD；不同 MOD04 版本可能只包含其中部分数据集。
        for sds_name in sds_priority:
            try:
                sds = handle.select(sds_name)
            except Exception:
                continue
            attrs = sds.attributes()
            raw = sds[:].astype(np.float64)
            scale = float(attrs.get("scale_factor", 1.0))
            offset = float(attrs.get("add_offset", 0.0))
            fill = attrs.get("_FillValue", None)
            # 先应用比例因子和偏移量，再按填充值、有效范围及 AOT 范围筛选。
            data = raw * scale + offset
            mask = inside & np.isfinite(data)
            if fill is not None:
                mask &= raw != float(fill)
            valid_range = attrs.get("valid_range", None)
            if valid_range is not None:
                mask &= (data >= float(valid_range[0])) & (data <= float(valid_range[1]))
            mask &= (data > AOT_MIN_VALID) & (data <= AOT_MAX_VALID)
            values = data[mask]
            if values.size:
                logger.info("  MOD04 %s: %s, 足迹内 %d 个有效像元, "
                            "AOD550 中位数 %.4f",
                            os.path.basename(path), sds_name, values.size,
                            float(np.median(values)))
                return {
                    "file": os.path.abspath(path),
                    "sds": sds_name,
                    "values": values.astype(np.float64),
                }
        logger.warning("  MOD04 %s 有足迹内像元，但没有有效 AOD 像元",
                       os.path.basename(path))
        return None
    except Exception as exc:
        logger.warning("读取 MOD04 失败: %s (%s)", path, exc)
        return None
    finally:
        if handle is not None:
            try:
                handle.end()
            except Exception:
                pass


def read_modis_aot(modis_dir, date, ring, logger):
    """汇总当天所有 MOD04 文件中足迹内的 AOD550，返回中位数与来源信息。"""
    search_started = time.time()
    logger.info("搜索本地 MODIS AOD：目录=%s，日期=%s",
                modis_dir, date.strftime("%Y-%m-%d"))
    files = find_mod04_files(modis_dir, date)
    if not files:
        logger.warning("MODIS 未覆盖/无可用文件：目录 %s 中没有 %s 当天的 "
                       "MOD04/MYD04 产品",
                       modis_dir, date.strftime("%Y-%m-%d"))
        logger.info("MODIS AOD 搜索结束，用时 %s",
                    fmt_elapsed(time.time() - search_started))
        return None
    logger.info("找到 %d 个 MOD04_L2 文件（%s）", len(files),
                date.strftime("%Y-%m-%d"))
    values = []
    used = []
    for path in files:
        result = read_mod04_aod(path, ring, logger)
        if result is None:
            continue
        values.append(result["values"])
        used.append({"file": result["file"], "sds": result["sds"],
                     "cells": int(result["values"].size)})
    if not values:
        logger.warning("MODIS 当天找到 %d 个产品，但没有产品覆盖本景足迹并提供有效 AOD550；"
                       "将尝试后续 AOT 来源", len(files))
        logger.info("MODIS AOD 搜索结束，用时 %s",
                    fmt_elapsed(time.time() - search_started))
        return None
    merged = np.concatenate(values)
    aot = float(np.median(merged))
    detail = {
        "source": "MODIS MOD04_L2",
        "aot550": aot,
        "cells": int(merged.size),
        "min": float(merged.min()),
        "max": float(merged.max()),
        "std": float(merged.std()),
        "files": used,
    }
    logger.info("AOT550 = %.4f（MODIS %d 个像元，范围 %.3f~%.3f）",
                aot, merged.size, merged.min(), merged.max())
    logger.info("MODIS AOD 搜索完成，用时 %s",
                fmt_elapsed(time.time() - search_started))
    return detail


_AERONET_SITE_CACHE = None
_AERONET_DISABLED = False   # 网络不可用/连续失败后，本次批处理不再尝试 AERONET


def _load_aeronet_sites(logger=None):
    """下载并缓存 AERONET 站点列表 [(name, lat, lon), ...]。

    官方文件列顺序为：Site_Name, Longitude, Latitude, Elevation。
    """
    global _AERONET_SITE_CACHE, _AERONET_DISABLED
    if _AERONET_DISABLED:
        return None
    if _AERONET_SITE_CACHE is not None:
        return _AERONET_SITE_CACHE
    try:
        import requests
    except ImportError:
        if logger:
            logger.info("未安装 requests，跳过 AERONET")
        return None
    url = "https://aeronet.gsfc.nasa.gov/aeronet_locations_v3.txt"
    try:
        if logger:
            logger.info("正在连接 AERONET 下载站点列表（最多等待 30 秒）")
        response = requests.get(url, timeout=30)
        response.raise_for_status()
    except Exception as exc:
        _AERONET_DISABLED = True
        if logger:
            logger.info("AERONET 站点列表下载失败（本次批处理不再尝试）: %s", exc)
        return None
    sites = []
    for line in response.text.strip().splitlines():
        parts = line.split(",")
        if len(parts) < 3:
            continue
        try:
            site_name = parts[0].strip()
            lon = float(parts[1])
            lat = float(parts[2])
        except ValueError:
            continue
        if -90.0 <= lat <= 90.0 and -180.0 <= lon <= 180.0:
            sites.append((site_name, lat, lon))
    _AERONET_SITE_CACHE = sites
    if logger:
        logger.info("AERONET 站点列表读取完成：%d 个站点", len(sites))
    return sites


def find_nearest_aeronet_site(center_lat, center_lon, max_distance_km=None,
                              logger=None):
    """从 AERONET 官方站点列表找最近站点（需联网，可选回退）。

    返回 (site_name, distance_km)；未找到时 (None, None)。
    """
    sites = _load_aeronet_sites(logger)
    if not sites:
        return None, None
    best_site, best_dist = None, 1e9
    for site_name, lat, lon in sites:
        dlat = np.radians(lat - center_lat)
        dlon = np.radians(lon - center_lon)
        a = (np.sin(dlat / 2.0) ** 2
             + np.cos(np.radians(center_lat)) * np.cos(np.radians(lat))
             * np.sin(dlon / 2.0) ** 2)
        dist = 6371.0 * 2.0 * np.arcsin(np.sqrt(a))
        if dist < best_dist:
            best_site, best_dist = site_name, dist
    if best_site is None:
        return None, None
    return best_site, best_dist


def _fetch_aeronet_aod_level(site_name, target_date, search_days, level,
                             logger=None):
    """按指定 AERONET level 下载 AOD_500nm 中位数，返回 (中位数, 记录数) 或 None。"""
    global _AERONET_DISABLED
    try:
        import requests
    except ImportError:
        return None
    # 注意：AERONET CGI 的参数名区分大小写，必须用 AOD20 / AOD15（小写会返回空结果）
    aod_key = "AOD20" if level == "2.0" else "AOD15"
    url = "https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3"
    for delta in range(search_days + 1):
        day1 = target_date - datetime.timedelta(days=delta)
        day2 = target_date + datetime.timedelta(days=delta)
        params = {
            "site": site_name,
            "year": day1.year, "month": day1.month, "day": day1.day,
            "year2": day2.year, "month2": day2.month, "day2": day2.day,
            aod_key: 1, "AVG": 10,
        }
        try:
            if logger:
                logger.info("请求 AERONET %s：Level %s，日期窗口 %s 至 %s "
                            "（最多等待 60 秒）",
                            site_name, level, day1.isoformat(), day2.isoformat())
            response = requests.get(url, params=params, timeout=60)
            response.raise_for_status()
        except Exception as exc:
            # 网络不可用时直接停用，避免批处理每景都等超时
            _AERONET_DISABLED = True
            if logger:
                logger.warning("AERONET 请求失败，将停止后续联网查询：%s", exc)
            break
        lines = response.text.strip().splitlines()
        header_idx = None
        for idx, line in enumerate(lines):
            if "AOD_500nm" in line:
                header_idx = idx
                break
        if header_idx is None:
            continue
        header = [item.strip() for item in lines[header_idx].split(",")]
        try:
            col_idx = header.index("AOD_500nm")
        except ValueError:
            continue
        values = []
        for line in lines[header_idx + 1:]:
            parts = line.split(",")
            if len(parts) <= col_idx:
                continue
            try:
                value = float(parts[col_idx])
            except ValueError:
                continue
            if -999 < value < 5:
                values.append(value)
        if len(values) >= 3:
            if logger:
                logger.info("AERONET Level %s 获取到 %d 条有效记录",
                            level, len(values))
            return float(np.median(values)), len(values)
        if logger:
            logger.info("AERONET 当前日期窗口没有足够有效记录")
    return None


def fetch_aeronet_aod(site_name, target_date, search_days=3, prefer_level="2.0",
                      logger=None):
    """从 AERONET 下载站点附近日期的 AOD_500nm（6S 用 AOT550 近似）。

    返回 (AOD, level, records)；无数据时返回 (None, None, None)。
    优先用 Level 2.0，没有数据时自动退回 Level 1.5。
    """
    levels = ("2.0", "1.5") if prefer_level == "2.0" else ("1.5",)
    for level in levels:
        result = _fetch_aeronet_aod_level(site_name, target_date,
                                          search_days, level, logger=logger)
        if result is not None:
            return result[0], level, result[1]
    return None, None, None


def resolve_aot550(meta, ring, args, logger):
    """确定 AOT550：命令行 > 本地 MODIS > AERONET > 手工兜底。"""
    if args.aot is not None:
        if not (AOT_MIN_VALID < args.aot <= AOT_MAX_VALID):
            raise ValueError("--aot 超出有效范围 (%.2f, %.2f]: %s"
                             % (AOT_MIN_VALID, AOT_MAX_VALID, args.aot))
        logger.info("AOT550 = %.4f（命令行指定）", args.aot)
        return float(args.aot), {"source": "命令行 --aot", "aot550": float(args.aot)}

    mode = (args.aot_source or AOT_SOURCE or "auto").lower()
    if mode not in ("auto", "modis", "aeronet", "manual"):
        raise ValueError("AOT_SOURCE 只能是 auto/modis/aeronet/manual: %s" % AOT_SOURCE)

    # auto 模式先使用本地 MODIS；明确选择 modis 时不回退到其他来源。
    if mode in ("auto", "modis"):
        detail = read_modis_aot(MODIS_DIR, meta["date"], ring, logger)
        if detail is not None:
            return detail["aot550"], detail

    # auto 模式下 MODIS 无有效覆盖时，按配置尝试联网 AERONET。
    if mode in ("auto", "aeronet") and AERONET_ENABLED \
            and not _AERONET_DISABLED and meta.get("center_lat") is not None:
        logger.info("尝试 AERONET 回退 ...")
        try:
            site, dist = find_nearest_aeronet_site(
                meta["center_lat"], meta["center_lon"], logger=logger)
            if site is None or dist is None:
                logger.info("AERONET 未找到可用站点")
            elif dist > AERONET_MAX_DISTANCE_KM:
                logger.info("AERONET 最近站点 %s 距离 %.0f km，超过上限 %d km，跳过",
                            site, dist, AERONET_MAX_DISTANCE_KM)
            else:
                aod, level, records = fetch_aeronet_aod(
                    site, meta["date"].date(), AERONET_SEARCH_DAYS,
                    logger=logger)
                if aod is not None:
                    logger.info("AOT550 = %.4f（AERONET %s，距离 %.0f km，"
                                "Level %s，%d 条记录，AOD500 近似 550 nm）",
                                aod, site, dist, level, records)
                    return float(aod), {
                        "source": "AERONET", "site": site,
                        "distance_km": float(dist), "aod500": float(aod),
                        "level": level, "records": int(records),
                        "note": "AERONET AOD_500nm 直接近似为 AOT550",
                    }
                logger.info("AERONET %s 站附近 %d 天无可用 AOD",
                            site, AERONET_SEARCH_DAYS)
        except Exception as exc:
            logger.warning("AERONET 回退失败: %s", exc)

    logger.warning("AOT550 自动获取失败，使用手工兜底值 %.4f", AOT550_MANUAL)
    return float(AOT550_MANUAL), {"source": "AOT550_MANUAL",
                                  "aot550": float(AOT550_MANUAL)}


# =================================================================
# 海拔：本地 ASTGTM DEM
# =================================================================

def find_dem_datasets(dem_dir, logger):
    """列出 DEM 目录下的栅格数据集路径（zip 内 .img 通过 /vsizip/ 读取）。"""
    if not dem_dir or not os.path.isdir(dem_dir):
        return []
    paths = []
    for path in sorted(glob.glob(os.path.join(dem_dir, "*"))):
        lower = path.lower()
        if lower.endswith(".zip"):
            try:
                with zipfile.ZipFile(path) as archive:
                    members = [name for name in archive.namelist()
                               if name.lower().endswith((".img", ".tif", ".tiff"))
                               and not name.endswith("/")]
                if not members:
                    logger.debug("DEM 压缩包内没有栅格: %s", path)
                    continue
                members.sort()
                paths.append("/vsizip/" + os.path.abspath(path).replace("\\", "/")
                             + "/" + members[0])
            except Exception as exc:
                logger.warning("无法读取 DEM 压缩包: %s (%s)", path, exc)
        elif lower.endswith((".tif", ".tiff", ".img")):
            paths.append(path)
    return paths


def compute_mean_elevation(ring, dem_dir, logger):
    """仅统计整景足迹多边形内的有效 DEM 像元并返回高程统计。"""
    started = time.time()
    logger.info("开始从 DEM 计算整景足迹平均海拔：%s", dem_dir)
    paths = find_dem_datasets(dem_dir, logger)
    if not paths:
        logger.warning("未找到 DEM 数据: %s", dem_dir)
        return None
    logger.info("DEM 数据: %d 个瓦片", len(paths))

    vrt_path = "/vsimem/zy1e_step2_dem.vrt"
    vrt = None
    try:
        vrt = gdal.BuildVRT(vrt_path, paths)
        if vrt is None:
            raise RuntimeError("DEM VRT 构建失败")
        first = gdal.Open(paths[0])
        if first is None:
            raise RuntimeError("无法打开 DEM 瓦片: %s" % paths[0])
        srs_wkt = first.GetProjection()
        geo_transform = first.GetGeoTransform()
        first = None
        if not srs_wkt:
            logger.warning("DEM 缺少投影信息，跳过自动海拔")
            return None

        srs_dem = osr.SpatialReference()
        srs_dem.ImportFromWkt(srs_wkt)
        srs_wgs84 = osr.SpatialReference()
        srs_wgs84.ImportFromEPSG(4326)
        srs_dem.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        srs_wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transformer = osr.CoordinateTransformation(srs_wgs84, srs_dem)
        points = [transformer.TransformPoint(lon, lat) for lon, lat in ring]
        min_x = min(point[0] for point in points)
        max_x = max(point[0] for point in points)
        min_y = min(point[1] for point in points)
        max_y = max(point[1] for point in points)

        dataset = gdal.Warp(
            "", vrt, format="MEM",
            dstSRS=srs_dem.ExportToWkt(),
            outputBounds=[min_x, min_y, max_x, max_y],
            xRes=DEM_RES, yRes=DEM_RES,
            resampleAlg="bilinear",
            dstNodata=32767,
            multithread=True,
            warpOptions=["NUM_THREADS=ALL_CPUS"])
        if dataset is None:
            raise RuntimeError("DEM 范围裁剪失败")
        elevation = dataset.GetRasterBand(1).ReadAsArray().astype(np.float64)
        output_gt = dataset.GetGeoTransform()
        logger.info("DEM 足迹矩形网格读取完成：%d×%d @ %g m",
                    dataset.RasterXSize, dataset.RasterYSize, DEM_RES)
        dataset = None
    finally:
        vrt = None
        try:
            gdal.Unlink(vrt_path)
        except Exception:
            pass

    # 将 DEM 像元中心投影到同一坐标系，并排除足迹外接矩形中的足迹外像元。
    ring_dem = [[point[0], point[1]]
                for point in (transformer.TransformPoint(lon, lat)
                              for lon, lat in ring)]
    cols = np.arange(elevation.shape[1], dtype=np.float64) + 0.5
    rows = np.arange(elevation.shape[0], dtype=np.float64) + 0.5
    x_grid = (output_gt[0] + cols[None, :] * output_gt[1]
              + rows[:, None] * output_gt[2])
    y_grid = (output_gt[3] + cols[None, :] * output_gt[4]
              + rows[:, None] * output_gt[5])
    footprint_mask = points_in_polygon(x_grid, y_grid, ring_dem)

    # 排除常见 DEM 填充值及超出合理地形高程范围的值。
    valid = footprint_mask & np.isfinite(elevation)
    valid &= (elevation != 32767) & (elevation != -32768)
    valid &= (elevation > -1000) & (elevation < 9000)
    total = int(footprint_mask.sum())
    count = int(valid.sum())
    if count == 0:
        logger.warning("DEM 整景足迹多边形内没有有效像元")
        return None
    mean_m = float(elevation[valid].mean())
    coverage = 100.0 * count / max(total, 1)
    detail = {
        "source": "ASTGTM DEM (%d tiles)" % len(paths),
        "tiles": paths,
        "elevation_km": mean_m / 1000.0,
        "mean_m": mean_m,
        "std_m": float(elevation[valid].std()),
        "min_m": float(elevation[valid].min()),
        "max_m": float(elevation[valid].max()),
        "valid_pixels": count,
        "coverage_percent": coverage,
        "res_m": DEM_RES,
        "dem_gt": [float(value) for value in geo_transform[:6]],
    }
    logger.info("DEM 整景足迹平均海拔 = %.1f m（范围 %.1f~%.1f m，有效覆盖 %.1f%%，"
                "%d 像元 @ %g m）", mean_m, detail["min_m"], detail["max_m"],
                coverage, count, DEM_RES)
    logger.info("DEM 海拔统计完成，用时 %s",
                fmt_elapsed(time.time() - started))
    return detail


def resolve_elevation(meta, ring, args, logger):
    """确定目标海拔（km）：命令行 > DEM > 手工兜底。"""
    if args.elev is not None:
        logger.info("海拔 = %.3f km（命令行指定）", args.elev)
        return float(args.elev), {"source": "命令行 --elev",
                                  "elevation_km": float(args.elev)}
    if AUTO_DETECT_DEM:
        detail = compute_mean_elevation(ring, DEM_DIR, logger)
        if detail is not None:
            return detail["elevation_km"], detail
    logger.warning("DEM 自动计算失败，使用兜底海拔 %.3f km", ELEVATION_KM_FALLBACK)
    return float(ELEVATION_KM_FALLBACK), {"source": "ELEVATION_KM_FALLBACK",
                                          "elevation_km": float(ELEVATION_KM_FALLBACK)}


# =================================================================
# 6S 大气模式、气溶胶及系数计算
# =================================================================

def auto_select_atmos_profile(date, center_lat, logger):
    """按成像日期和中心纬度选择 6S 大气模式（南半球月份翻转）。"""
    if center_lat is None:
        logger.warning("无中心纬度，使用 %s", ATMOS_PROFILE_MANUAL)
        return ATMOS_PROFILE_MANUAL
    lat_abs = abs(center_lat)
    is_summer = 4 <= date.month <= 9
    if center_lat < 0:
        is_summer = not is_summer
    if lat_abs < 23.0:
        profile = "Tropical"
    elif lat_abs <= 55.0:
        profile = "MidlatitudeSummer" if is_summer else "MidlatitudeWinter"
    else:
        profile = "SubarcticSummer" if is_summer else "SubarcticWinter"
    logger.info("大气模式: 纬度 %.2f°, 月份 %d → %s（自动）",
                lat_abs, date.month, profile)
    return profile


def resolve_atmos_profile(meta, args, logger):
    if args.atmos is not None:
        return args.atmos, "命令行 --atmos"
    if WATER_VAPOR_CM is not None:
        name = "UserWaterAndOzone(water=%.2f, ozone=%.2f)" % (
            float(WATER_VAPOR_CM),
            float(OZONE_CM if OZONE_CM is not None else OZONE_CM_DEFAULT))
        return name, "配置区 WATER_VAPOR_CM"
    if AUTO_DETECT_ATMOS and meta.get("date") is not None:
        return auto_select_atmos_profile(meta["date"], meta.get("center_lat"),
                                         logger), "自动（纬度+月份）"
    return ATMOS_PROFILE_MANUAL, "配置区 ATMOS_PROFILE_MANUAL"


def guess_aero_profile(meta, logger):
    """无离线土地覆盖数据时的经纬度粗判（默认关闭，仅预留）。"""
    lat = meta.get("center_lat")
    lon = meta.get("center_lon")
    if lat is None or lon is None:
        return AERO_PROFILE_MANUAL
    lat_abs = abs(lat)
    desert = (75 < lon < 110) and (35 < lat_abs < 45)
    big_city = ((115 < lon < 120 and 38 < lat_abs < 41)
                or (118 < lon < 123 and 29 < lat_abs < 33)
                or (112 < lon < 115 and 21 < lat_abs < 24))
    if desert:
        profile = "Desert"
    elif big_city:
        profile = "Urban"
    else:
        profile = AERO_PROFILE_MANUAL
    logger.info("气溶胶类型: 经纬度粗判 → %s", profile)
    return profile


def resolve_aero_profile(meta, args, logger):
    if args.aero is not None:
        return args.aero, "命令行 --aero"
    if AUTO_DETECT_AERO:
        return guess_aero_profile(meta, logger), "自动（经纬度粗判）"
    return AERO_PROFILE_MANUAL, "配置区 AERO_PROFILE_MANUAL"


def resolve_sixs_exe(logger):
    """查找 6S 可执行文件。"""
    candidates = []
    if SIXS_EXE:
        candidates.append(SIXS_EXE)
    env_path = os.environ.get("SIXS_EXE")
    if env_path:
        candidates.append(env_path)
    # Conda on Windows usually installs executables in the active environment's Library/bin.
    conda_bin = os.path.join(sys.prefix, "Library", "bin")
    for name in ("sixs.exe", "sixs", "sixsV1.1.exe", "sixsV1.1"):
        candidates.append(os.path.join(conda_bin, name))
    script_dir = os.path.dirname(os.path.abspath(__file__))
    for name in ("sixs.exe", "sixs", "sixsV1.1.exe", "sixsV1.1"):
        candidates.append(os.path.join(script_dir, "_bin", name))
    for name in ("sixs.exe", "sixs", "sixsV1.1.exe", "sixsV1.1"):
        found = shutil.which(name)
        if found:
            candidates.append(found)

    for candidate in candidates:
        if candidate and os.path.isfile(candidate):
            logger.info("6S 可执行文件: %s", candidate)
            return os.path.abspath(candidate)
    raise FileNotFoundError(
        "未找到 6S 可执行文件 sixs.exe。\n"
        "请把 sixs.exe（及 libgfortran/libgcc/libquadmath/libwinpthread DLL）放到 %s，\n"
        "或设置环境变量 SIXS_EXE 指向可执行文件。"
        % os.path.join(script_dir, "_bin"))


def build_sixs(sixs_exe, meta, atmos_name, aero_name, aot550,
               elevation_km, logger):
    """构造 Py6S 对象（复用同一个对象逐波段 run，避免重复查找可执行文件）。"""
    try:
        from Py6S import (SixS, Geometry, AtmosProfile, AeroProfile,
                          Wavelength, AtmosCorr, GroundReflectance)
    except ImportError as exc:
        raise ImportError(
            "未安装 Py6S，请用装有 GDAL+Py6S 的 Python 运行本脚本"
            "（例如 E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe）: %s" % exc)

    atmos_map = {
        "Tropical": AtmosProfile.Tropical,
        "MidlatitudeSummer": AtmosProfile.MidlatitudeSummer,
        "MidlatitudeWinter": AtmosProfile.MidlatitudeWinter,
        "SubarcticSummer": AtmosProfile.SubarcticSummer,
        "SubarcticWinter": AtmosProfile.SubarcticWinter,
        "USStandard": AtmosProfile.USStandard1962,
    }
    aero_map = {
        "Continental": AeroProfile.Continental,
        "Urban": AeroProfile.Urban,
        "Maritime": AeroProfile.Maritime,
        "Desert": AeroProfile.Desert,
    }

    sixs = SixS(sixs_exe)
    sixs.geometry = Geometry.User()
    sixs.geometry.solar_z = float(meta["solar_z"])
    sixs.geometry.solar_a = float(meta["solar_a"])
    sixs.geometry.view_z = float(meta["view_z"])
    sixs.geometry.view_a = float(meta["view_a"])
    if meta.get("date") is not None:
        sixs.geometry.day = int(meta["date"].day)
        sixs.geometry.month = int(meta["date"].month)

    if isinstance(atmos_name, str) and atmos_name.startswith("UserWaterAndOzone"):
        sixs.atmos_profile = AtmosProfile.UserWaterAndOzone(
            float(WATER_VAPOR_CM),
            float(OZONE_CM if OZONE_CM is not None else OZONE_CM_DEFAULT))
    else:
        if atmos_name not in atmos_map:
            raise ValueError("不支持的大气模式: %s（可选: %s）"
                             % (atmos_name, ", ".join(sorted(atmos_map))))
        sixs.atmos_profile = AtmosProfile.PredefinedType(atmos_map[atmos_name])

    if aero_name not in aero_map:
        raise ValueError("不支持的气溶胶类型: %s（可选: %s）"
                         % (aero_name, ", ".join(sorted(aero_map))))
    sixs.aero_profile = AeroProfile.PredefinedType(aero_map[aero_name])
    sixs.aot550 = float(aot550)
    sixs.altitudes.set_sensor_satellite_level()
    sixs.altitudes.set_target_custom_altitude(float(elevation_km))
    sixs.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(0)
    # 地面反射率对 Lambertian 大气校正系数没有影响，显式设为 0 便于阅读 6S 输入文件
    sixs.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.0)
    return sixs


# 6S（1.1.x 的某些构建）在个别波长会输出 NaN 系数（本数据实测 0.5506/0.6708/
# 1.2408 μm 附近约 1 nm 宽度内）。这些波段的系数与相邻波长几乎一致，因此按
# 微小平移重试，失败时再用相邻有效波段线性插值补齐。
SIXS_WL_OFFSETS = (0.001, -0.001, 0.002, -0.002, 0.003, -0.003,
                   0.005, -0.005, -0.008, 0.008)


def _sixs_coefficients_once(sixs, wavelength):
    """在指定波长跑一次 6S，返回 (xa, xb, xc)。"""
    from Py6S import Wavelength
    sixs.wavelength = Wavelength(float(wavelength))
    sixs.run()
    return np.asarray([float(sixs.outputs.coef_xa),
                       float(sixs.outputs.coef_xb),
                       float(sixs.outputs.coef_xc)], dtype=np.float64)


def compute_sixs_coefficients(sixs, wavelengths, logger):
    """逐波段运行 6S，返回系数矩阵 (n, 3)、6S 版本、耗时与重试波段信息。"""
    n_bands = int(len(wavelengths))
    logger.info("开始计算 6S 系数：共 %d 个波段", n_bands)
    coefficients = np.full((n_bands, 3), np.nan, dtype=np.float64)
    version = None
    perturbed = {}
    started = time.time()

    for index, wavelength in enumerate(wavelengths, start=1):
        values = None
        error = None
        try:
            values = _sixs_coefficients_once(sixs, wavelength)
        except Exception as exc:
            error = exc
        if values is not None and np.all(np.isfinite(values)) and values[0] > 0:
            coefficients[index - 1] = values
        else:
            # 6S 在个别波长输出 NaN/异常：按微小平移重试（系数随波长变化平缓）
            for offset in SIXS_WL_OFFSETS:
                candidate = float(wavelength) + offset
                if candidate <= 0.0:
                    continue
                try:
                    values = _sixs_coefficients_once(sixs, candidate)
                except Exception:
                    continue
                if np.all(np.isfinite(values)) and values[0] > 0:
                    coefficients[index - 1] = values
                    perturbed[index] = (float(candidate), offset)
                    logger.warning(
                        "6S 在 %.6f μm 无有效输出，改用 %.6f μm"
                        "（偏移 %+.1f nm）", float(wavelength), candidate,
                        offset * 1000.0)
                    break
            if np.any(~np.isfinite(coefficients[index - 1])):
                logger.error("6S 在 %.6f μm 多次重试仍失败%s，稍后用相邻波段插值",
                             float(wavelength),
                             "（%s）" % error if error else "")

        if version is None and getattr(sixs.outputs, "version", None):
            version = sixs.outputs.version
        if index == 1 or index % 10 == 0 or index == n_bands:
            elapsed = time.time() - started
            eta = elapsed / index * (n_bands - index)
            logger.info("  6S 计算进度 %3d/%d 波段（%.1f%%），已用 %s，预计剩余 %s",
                        index, n_bands, 100.0 * index / n_bands,
                        fmt_elapsed(elapsed), fmt_elapsed(eta))

    # 仍失败的波段：用相邻有效波段按波长线性插值补齐
    invalid = np.where(~np.isfinite(coefficients).all(axis=1))[0]
    if invalid.size:
        valid = np.where(np.isfinite(coefficients).all(axis=1))[0]
        if valid.size < 2:
            raise ValueError("6S 系数几乎全部无效（有效波段 %d 个）" % valid.size)
        order = np.argsort(np.asarray(wavelengths)[valid])
        sorted_wl = np.asarray(wavelengths)[valid][order]
        for column in range(3):
            sorted_values = coefficients[valid, column][order]
            coefficients[invalid, column] = np.interp(
                np.asarray(wavelengths)[invalid], sorted_wl, sorted_values)
        logger.warning("6S 系数用相邻有效波段线性插值补齐 %d 个波段: %s",
                       invalid.size,
                       ", ".join(str(int(pos) + 1) for pos in invalid))

    if not np.all(np.isfinite(coefficients)):
        raise ValueError("6S 系数仍包含 NaN/Inf，结果不可用")
    if np.any(coefficients[:, 0] <= 0):
        bad = np.where(coefficients[:, 0] <= 0)[0]
        raise ValueError("6S 系数 xa <= 0（波段 %s），结果不可用"
                         % ", ".join(str(int(pos) + 1) for pos in bad[:5]))

    if perturbed:
        logger.warning("共 %d 个波段使用了微小平移重试: %s", len(perturbed),
                       ", ".join("波段%d" % band for band in sorted(perturbed)))
    elapsed = time.time() - started
    logger.info("6S 系数计算完成: %d 波段，用时 %s", n_bands, fmt_elapsed(elapsed))
    return coefficients, version, elapsed, perturbed


# =================================================================
# 反射率写出与统计
# =================================================================

def build_gtiff_creation_options():
    creation = list(GTIFF_OPTIONS)
    creation.append("INTERLEAVE=BAND")   # 高光谱立方体用 BSQ：逐波段写入效率高
    if GEOTIFF_COMPRESS:
        creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS,
                         "PREDICTOR=3",
                         "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL])
    return creation


def write_reflectance(src_ds, coefficients, wavelengths, tmp_path, logger):
    """逐波段把正射辐亮度换算成反射率并写出同网格 GeoTIFF，返回逐波段统计。

    输出直接沿用输入正射辐亮度的仿射变换与投影（CRS/尺寸完全一致），
    不做重采样、不写 RPC；步骤1 已完成正射。
    """
    n_bands, n_rows, n_cols = (src_ds.RasterCount,
                               src_ds.RasterYSize, src_ds.RasterXSize)
    if n_bands != int(len(wavelengths)):
        raise ValueError("辐亮度波段数 %d != 波长数 %d"
                         % (n_bands, len(wavelengths)))
    if n_bands != coefficients.shape[0]:
        raise ValueError("辐亮度波段数 %d != 6S 系数行数 %d"
                         % (n_bands, coefficients.shape[0]))

    logger.info("开始反射率逐波段处理：%d 波段，影像 %d×%d，临时输出 %s",
                n_bands, n_cols, n_rows, tmp_path)
    driver = gdal.GetDriverByName("GTiff")
    if os.path.exists(tmp_path):
        os.remove(tmp_path)
    out_ds = driver.Create(tmp_path, n_cols, n_rows, n_bands,
                           gdal.GDT_Float32,
                           options=build_gtiff_creation_options())
    if out_ds is None:
        raise RuntimeError("无法创建 GeoTIFF: %s" % tmp_path)
    # 与输入正射辐亮度严格同网格：直接复制仿射变换与投影（不重采样、不写 RPC）
    out_ds.SetGeoTransform(src_ds.GetGeoTransform())
    out_ds.SetProjection(src_ds.GetProjection())

    clip_low, clip_high = (REFLECTANCE_CLIP if REFLECTANCE_CLIP else (None, None))
    stats = []
    masked_pixels = 0
    total_pixels = int(n_rows) * int(n_cols)

    for index, (wavelength, coef) in enumerate(
            zip(wavelengths, coefficients), start=1):
        source_band = src_ds.GetRasterBand(index)
        radiance = source_band.ReadAsArray().astype(np.float64)
        invalid = ~np.isfinite(radiance)
        if MASK_INVALID_RADIANCE:
            invalid |= radiance <= float(RADIANCE_MIN_VALID)

        y_value = coef[0] * radiance - coef[1]
        denominator = 1.0 + coef[2] * y_value
        with np.errstate(divide="ignore", invalid="ignore"):
            reflectance = y_value / denominator
        invalid |= ~np.isfinite(reflectance) | (denominator <= 0)
        if clip_low is not None:
            np.clip(reflectance, clip_low, clip_high, out=reflectance)

        count = int(invalid.sum())
        if count:
            reflectance[invalid] = NODATA_OUT
            masked_pixels += count

        out_band = out_ds.GetRasterBand(index)
        out_band.WriteArray(reflectance.astype(np.float32))
        out_band.SetNoDataValue(NODATA_OUT)
        out_band.SetDescription("%.4f um" % float(wavelength))

        valid_values = reflectance[~invalid]
        if valid_values.size:
            sample = valid_values[::4]
            p1, p99 = np.percentile(sample, [1, 99])
            stats.append({
                "band": index,
                "wavelength": float(wavelength),
                "min": float(valid_values.min()),
                "max": float(valid_values.max()),
                "mean": float(valid_values.mean()),
                "std": float(valid_values.std()),
                "p1": float(p1),
                "p99": float(p99),
                "valid": int(valid_values.size),
            })
        else:
            logger.warning("  波段 %d 没有有效像元", index)
            stats.append({"band": index, "wavelength": float(wavelength),
                          "min": 0.0, "max": 0.0, "mean": 0.0, "std": 0.0,
                          "p1": 0.0, "p99": 0.0, "valid": 0})

        if index == 1 or index % 10 == 0 or index == n_bands:
            logger.info("  反射率处理进度 %3d/%d 波段（%.1f%%），当前波段有效像元 %.1f%%",
                        index, n_bands, 100.0 * index / n_bands,
                        100.0 * stats[-1]["valid"] / max(total_pixels, 1))

    out_ds.FlushCache()
    out_ds = None
    if MASK_INVALID_RADIANCE and masked_pixels:
        logger.info("  辐亮度无效像元累计: %d 个", masked_pixels)
    return stats


def write_band_stats_csv(stats, output_path):
    tmp_path = output_path + ".tmp"
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "wavelength_um", "min", "max", "mean",
                         "std", "p1", "p99", "valid_pixels"])
        for item in stats:
            writer.writerow([
                item["band"],
                "%.6f" % item["wavelength"],
                "%.6f" % item["min"],
                "%.6f" % item["max"],
                "%.6f" % item["mean"],
                "%.6f" % item["std"],
                "%.6f" % item["p1"],
                "%.6f" % item["p99"],
                item["valid"],
            ])
    os.replace(tmp_path, output_path)


def read_band_stats_csv(path):
    stats = []
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        for row in csv.DictReader(handle):
            try:
                stats.append({
                    "band": int(row["band"]),
                    "wavelength": float(row["wavelength_um"]),
                    "min": float(row["min"]),
                    "max": float(row["max"]),
                    "mean": float(row["mean"]),
                    "std": float(row["std"]),
                    "p1": float(row["p1"]),
                    "p99": float(row["p99"]),
                    "valid": int(row["valid_pixels"]),
                })
            except (KeyError, ValueError):
                continue
    return stats


# =================================================================
# 质检：带坐标预览与成果图
# =================================================================

def select_qc_bands(wavelengths):
    """选择质检用的 5 个波段（真彩色 / 假彩色 / 短波红外）。"""
    band_map = {
        "R": _nearest_band(wavelengths, 0.66),
        "G": _nearest_band(wavelengths, 0.55),
        "B": _nearest_band(wavelengths, 0.47),
        "NIR": _nearest_band(wavelengths, 0.86),
        "SWIR": _nearest_band(wavelengths, 2.20),
    }
    unique = sorted(set(band_map.values()))
    return band_map, unique


def _nearest_band(wavelengths, target):
    return int(np.argmin(np.abs(wavelengths - target)))


def _stretch(array, valid):
    values = array[valid]
    if values.size == 0:
        return np.zeros_like(array)
    low, high = np.percentile(values, [2, 98])
    if high <= low:
        high = low + 1e-6
    return np.clip((array - low) / (high - low), 0.0, 1.0)


def prepare_qc_preview(reflectance_path, wavelengths, preview_tif, logger):
    """对正射反射率的 5 个预览波段直接降采样出带坐标预览，不做 RPC/warp。

    反射率已与步骤1 正射辐亮度同网格（EPSG:32650/30 m），这里只用
    gdal.Translate 按 QC_PREVIEW_SCALE 降采样（默认 60 m），保留原 CRS。
    """
    band_map, unique = select_qc_bands(wavelengths)
    source = gdal.Open(reflectance_path)
    if source is None:
        raise RuntimeError("无法打开反射率: %s" % reflectance_path)
    try:
        geo_transform = source.GetGeoTransform()
        if not geo_transform or abs(geo_transform[1]) <= 0.0:
            raise RuntimeError("反射率缺少有效仿射变换，无法生成带坐标预览")
        preview_res = abs(geo_transform[1]) * QC_PREVIEW_SCALE
        options = gdal.TranslateOptions(
            format="GTiff",
            bandList=[band + 1 for band in unique],
            xRes=preview_res,
            yRes=preview_res,
            resampleAlg="bilinear",
            noData=NODATA_OUT,
            outputType=gdal.GDT_Float32,
            creationOptions=build_gtiff_creation_options(),
        )
        tmp_path = preview_tif + ".tmp.tif"
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        result = gdal.Translate(tmp_path, source, options=options)
        if result is None:
            raise RuntimeError("预览降采样失败: %s" % preview_tif)
        result = None
        update = gdal.Open(tmp_path, gdal.GA_Update)
        for position, band_index in enumerate(unique, start=1):
            band = update.GetRasterBand(position)
            band.SetNoDataValue(NODATA_OUT)
            band.SetDescription("%.4f um" % wavelengths[band_index])
        update = None
        os.replace(tmp_path, preview_tif)
        logger.info("  带坐标预览: %s (%d 波段 @ %g m，降采样自正射反射率)",
                    preview_tif, len(unique), preview_res)
        return band_map, unique
    finally:
        source = None


def _prepare_figure_data(tif_path, preview_tif, band_map, unique,
                         wavelengths, logger):
    """返回 (arrays, extent, xlabel, ylabel, crs_text, size_text, georeferenced)。"""
    if preview_tif and os.path.isfile(preview_tif):
        dataset = gdal.Open(preview_tif)
        if dataset is not None:
            arrays = {}
            for name, band_index in band_map.items():
                position = unique.index(band_index) + 1
                arrays[name] = dataset.GetRasterBand(position).ReadAsArray(
                ).astype(np.float64)
            width, height = dataset.RasterXSize, dataset.RasterYSize
            geo_transform = dataset.GetGeoTransform()
            srs = osr.SpatialReference()
            srs.ImportFromWkt(dataset.GetProjection())
            epsg = srs.GetAuthorityCode(None) or "?"
            extent = [geo_transform[0] / 1000.0,
                      (geo_transform[0] + geo_transform[1] * width) / 1000.0,
                      (geo_transform[3] + geo_transform[5] * height) / 1000.0,
                      geo_transform[3] / 1000.0]
            dataset = None
            crs_text = "EPSG:%s @ %g m preview" % (epsg, abs(geo_transform[1]))
            return (arrays, extent, "Easting (km)", "Northing (km)",
                    crs_text, "%d×%d" % (width, height), True)

    # 回退：无带坐标预览时直接读正射反射率（降采样），坐标轴只画行列号
    dataset = gdal.Open(tif_path)
    if dataset is None:
        return None
    width, height = dataset.RasterXSize, dataset.RasterYSize
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))
    buf_w, buf_h = max(1, width // scale), max(1, height // scale)
    arrays = {}
    for name, band_index in band_map.items():
        arrays[name] = dataset.GetRasterBand(band_index + 1).ReadAsArray(
            buf_xsize=buf_w, buf_ysize=buf_h).astype(np.float64)
    dataset = None
    return (arrays, [0, width, height, 0], "Column", "Row",
            "no georeference", "%d×%d" % (width, height), False)


def write_quality_figure(tif_path, wavelengths, output_png, scene,
                         params_summary, logger, band_stats=None,
                         preview_tif=None):
    """生成质检成果图：真彩色 / 假彩色 / NDVI / 均值光谱。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        logger.warning("matplotlib 不可用，跳过质检图: %s", exc)
        return False

    wavelengths = np.asarray(wavelengths, dtype=np.float64)
    band_map, unique = select_qc_bands(wavelengths)
    prepared = _prepare_figure_data(tif_path, preview_tif, band_map, unique,
                                    wavelengths, logger)
    if prepared is None:
        logger.warning("无法读取成果影像，跳过质检图: %s", tif_path)
        return False
    arrays, extent, xlabel, ylabel, crs_text, size_text, georeferenced = prepared

    def combine(names):
        channels = [arrays[name] for name in names]
        valid = np.ones_like(channels[0], dtype=bool)
        for channel in channels:
            valid &= np.isfinite(channel) & (channel != NODATA_OUT)
        rgb = np.ones(channels[0].shape + (3,), dtype=np.float64)
        for index, channel in enumerate(channels):
            rgb[..., index] = _stretch(channel, valid)
        rgb[~valid] = 1.0
        return np.clip(rgb, 0.0, 1.0), valid

    true_rgb, valid = combine(("R", "G", "B"))
    false_rgb, _ = combine(("NIR", "R", "G"))

    nir = arrays["NIR"]
    red = arrays["R"]
    ndvi_valid = (np.isfinite(nir) & np.isfinite(red)
                  & (nir != NODATA_OUT) & (red != NODATA_OUT)
                  & ((nir + red) > 0))
    with np.errstate(divide="ignore", invalid="ignore"):
        ndvi = (nir - red) / (nir + red)
    ndvi_map = np.ma.masked_invalid(np.where(ndvi_valid, ndvi, np.nan))

    if band_stats:
        figure, axes = plt.subplots(2, 2, figsize=(13, 10))
        ax_true, ax_false, ax_ndvi, ax_spec = axes.ravel()
    else:
        figure, axes = plt.subplots(1, 3, figsize=(15, 5.5))
        ax_true, ax_false, ax_ndvi = axes.ravel()
        ax_spec = None

    wl_r = wavelengths[band_map["R"]]
    wl_g = wavelengths[band_map["G"]]
    wl_b = wavelengths[band_map["B"]]
    wl_nir = wavelengths[band_map["NIR"]]

    ax_true.imshow(true_rgb, extent=extent, origin="upper",
                   interpolation="nearest")
    ax_true.set_title("True color  (R %.2f / G %.2f / B %.2f μm)"
                      % (wl_r, wl_g, wl_b), fontsize=10)
    ax_false.imshow(false_rgb, extent=extent, origin="upper",
                    interpolation="nearest")
    ax_false.set_title("False color  (NIR %.2f / R %.2f / G %.2f μm)"
                       % (wl_nir, wl_r, wl_g), fontsize=10)
    ax_ndvi.imshow(ndvi_map, extent=extent, origin="upper", cmap="RdYlGn",
                   vmin=-0.2, vmax=0.9, interpolation="nearest")
    ax_ndvi.set_title("NDVI  (NIR %.2f, R %.2f μm)" % (wl_nir, wl_r),
                      fontsize=10)
    for axis in (ax_true, ax_false, ax_ndvi):
        axis.set_xlabel(xlabel, fontsize=9)
        axis.set_ylabel(ylabel, fontsize=9)
        axis.tick_params(labelsize=8)

    if ax_spec is not None:
        spectrum = sorted(band_stats, key=lambda item: item["wavelength"])
        spectrum_wl = np.array([item["wavelength"] for item in spectrum])
        spectrum_mean = np.array([item["mean"] for item in spectrum])
        spectrum_std = np.array([item["std"] for item in spectrum])
        ax_spec.plot(spectrum_wl, spectrum_mean, color="black", linewidth=1.0)
        ax_spec.fill_between(spectrum_wl,
                             spectrum_mean - spectrum_std,
                             spectrum_mean + spectrum_std,
                             color="0.75", alpha=0.4, linewidth=0)
        ax_spec.set_xlabel("Wavelength (μm)", fontsize=9)
        ax_spec.set_ylabel("Mean reflectance", fontsize=9)
        ax_spec.set_ylim(bottom=0.0)
        ax_spec.grid(alpha=0.3)
        ax_spec.set_title("Mean reflectance spectrum (full-res stats)",
                          fontsize=10)

    valid_ratio = float(valid.mean()) if valid.size else 0.0
    figure.suptitle("%s | %d bands | %s | %s\n%s | valid %.1f%%"
                    % (scene, len(wavelengths), size_text, crs_text,
                       params_summary, 100.0 * valid_ratio), fontsize=10)
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    figure.savefig(output_png, dpi=QC_DPI, bbox_inches="tight")
    plt.close(figure)
    logger.info("质检成果图: %s", output_png)
    return True


def params_summary_text(payload):
    return ("%s / %s / AOT550=%.3f / elev=%.3f km"
            % (payload.get("atmos_profile", "?"),
               payload.get("aero_profile", "?"),
               float(payload.get("aot550", 0.0)),
               float(payload.get("elevation_km", 0.0))))


def write_parameters(payload, json_path, txt_path, logger):
    atomic_write_json(payload, json_path)
    lines = [
        "# 步骤2 大气校正（6S）参数记录",
        "场景: %s" % payload.get("scene", ""),
        "正射辐亮度输入: %s" % payload.get("source_radiance", ""),
        "原始压缩包: %s" % (payload.get("source_tar") or "缺失"),
        "波段数: %d" % payload.get("band_count", 0),
        "输出网格: EPSG:%s @ %s m（与输入正射辐亮度一致，无重采样、无 RPC）"
        % (payload.get("output_grid_crs_epsg", "?"),
           payload.get("output_grid_res", "?")),
        "成像日期: %s" % (payload.get("date") or "未知"),
        "传感器: %s / %s" % (payload.get("satellite", "?"),
                             payload.get("sensor", "?")),
        "太阳天顶角: %s" % payload.get("solar_zenith"),
        "太阳方位角: %s" % payload.get("solar_azimuth"),
        "卫星天顶角: %s" % payload.get("view_zenith"),
        "卫星方位角: %s" % payload.get("view_azimuth"),
        "中心纬度: %s" % payload.get("center_lat"),
        "中心经度: %s" % payload.get("center_lon"),
        "大气模式: %s（%s）" % (payload.get("atmos_profile"),
                              payload.get("atmos_source")),
        "气溶胶类型: %s（%s）" % (payload.get("aero_profile"),
                               payload.get("aero_source")),
        "AOT550: %.4f（%s）" % (float(payload.get("aot550", 0.0)),
                              payload.get("aot_source")),
        "海拔(km): %.4f（%s）" % (float(payload.get("elevation_km", 0.0)),
                                payload.get("elevation_source")),
        "水汽/臭氧: %s" % (payload.get("water_vapor") or "6S 大气模式自带廓线"),
        "6S 可执行文件: %s" % payload.get("sixs_path"),
        "6S 版本: %s（Py6S %s）" % (payload.get("sixs_version"),
                                  payload.get("py6s_version")),
        "6S 波长重试: %s" % (
            ", ".join("波段%d(%+.1f nm)" % (item["band"], item["offset_nm"])
                      for item in payload.get("sixs_perturbed_bands") or [])
            or "无"),
        "校正公式: %s" % payload.get("formula"),
        "反射率裁剪: %s" % (payload.get("reflectance_clip") or "不裁剪"),
        "生成时间: %s" % payload.get("created"),
    ]
    atomic_write_text("\n".join(lines) + "\n", txt_path, )
    logger.info("参数记录: %s / %s", json_path, txt_path)


def regenerate_qc(scene_dir, scene, logger):
    """已有主干成果时，仅补齐质检图/预览。"""
    reflectance = os.path.join(scene_dir, scene + OUT_REFLECTANCE_SUFFIX)
    wavelengths_path = os.path.join(scene_dir, OUT_WAVELENGTHS_NAME)
    stats_path = os.path.join(scene_dir, OUT_STATS_NAME)
    qc_path = os.path.join(scene_dir, OUT_QC_NAME)
    preview_path = os.path.join(scene_dir, OUT_PREVIEW_GEO_NAME)
    try:
        wavelengths = np.load(wavelengths_path)
        band_stats = read_band_stats_csv(stats_path) if os.path.isfile(stats_path) else None
        params = {}
        params_path = os.path.join(scene_dir, OUT_PARAMS_JSON)
        if os.path.isfile(params_path):
            params = load_georef_json(params_path)
        preview_tif = None
        if QC_GEOREF_PREVIEW:
            try:
                prepare_qc_preview(reflectance, wavelengths, preview_path, logger)
                preview_tif = preview_path
            except Exception:
                logger.warning("补齐带坐标预览失败，质检图回退无坐标降采样",
                               exc_info=True)
        write_quality_figure(reflectance, wavelengths, qc_path, scene,
                             params_summary_text(params), logger,
                             band_stats=band_stats, preview_tif=preview_tif)
    except Exception:
        logger.warning("补齐质检图失败: %s", scene_dir, exc_info=True)


# =================================================================
# 单景处理
# =================================================================

def resolve_py6s_version():
    try:
        import Py6S
        return getattr(Py6S, "__version__", None) or \
            getattr(getattr(Py6S, "SixS", None), "__version__", None)
    except Exception:
        return None


def process_scene(scene_dir, out_root, args, logger):
    """处理一景：读步骤1 成果 → 6S 系数 → 逐波段校正 → 参数/统计/质检。"""
    scene = os.path.basename(os.path.normpath(scene_dir))
    in_radiance = resolve_radiance_path(scene_dir)
    in_wavelengths = os.path.join(scene_dir, IN_WAVELENGTHS_NAME)
    in_georef = os.path.join(scene_dir, IN_GEOREF_NAME)

    out_dir = os.path.join(out_root, scene)
    final_tif = os.path.join(out_dir, scene + OUT_REFLECTANCE_SUFFIX)
    legacy_tifs = [os.path.join(out_dir, name)
                   for name in LEGACY_REFLECTANCE_NAMES]
    final_wl = os.path.join(out_dir, OUT_WAVELENGTHS_NAME)
    final_georef = os.path.join(out_dir, OUT_GEOREF_NAME)
    final_coef = os.path.join(out_dir, OUT_COEF_NAME)
    final_stats = os.path.join(out_dir, OUT_STATS_NAME)
    final_params_json = os.path.join(out_dir, OUT_PARAMS_JSON)
    final_params_txt = os.path.join(out_dir, OUT_PARAMS_TXT)
    final_qc = os.path.join(out_dir, OUT_QC_NAME)
    final_preview = os.path.join(out_dir, OUT_PREVIEW_GEO_NAME)

    legacy_tif = next((path for path in legacy_tifs if os.path.isfile(path)),
                      None)
    if (not args.overwrite and not os.path.isfile(final_tif)
            and legacy_tif and os.path.isfile(final_coef)
            and os.path.isfile(final_params_json)):
        os.replace(legacy_tif, final_tif)
        logger.info("复用已有完整成果并更新反射率文件名：%s", final_tif)

    have_main = all(os.path.isfile(path)
                    for path in (final_tif, final_coef, final_params_json))
    if not args.overwrite and have_main:
        # 补齐从步骤1 复制的光谱/几何小文件，避免重复计算
        try:
            if not os.path.isfile(final_wl) and os.path.isfile(in_wavelengths):
                copy_file_atomic(in_wavelengths, final_wl)
            if not os.path.isfile(final_georef) and os.path.isfile(in_georef):
                copy_file_atomic(in_georef, final_georef)
        except OSError:
            logger.warning("补齐 wavelengths_clean.npy/georef.json 失败", exc_info=True)
        if MAKE_QC and not args.no_qc and not os.path.isfile(final_qc):
            regenerate_qc(out_dir, scene, logger)
        logger.info("跳过（成果已存在）: %s", out_dir)
        return {"scene": scene, "status": "跳过", "bands": None,
                "elapsed": 0.0, "message": "成果已存在"}

    started = time.time()
    current_stage = 1
    stage_description = "校验步骤1输入与正射网格"
    created_dir = not os.path.isdir(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    tmp_tif = os.path.join(out_dir, scene + "_cal_badband_ortho_atmcorr.tmp.tif")
    tmp_coef = os.path.join(out_dir, "py6s_coefficients.tmp.npy")
    tmp_wl = os.path.join(out_dir, "wavelengths_clean.tmp.npy")
    written_targets = []
    main_complete = False
    source_ds = None
    try:
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        if in_radiance is None:
            raise FileNotFoundError(
                "找不到步骤1正射辐亮度（*_cal_badband_ortho.tif 或旧格式 %s）: %s"
                % (LEGACY_RADIANCE_NAME, scene_dir))
        if not os.path.isfile(in_wavelengths):
            raise FileNotFoundError("找不到步骤1 波长文件: %s" % in_wavelengths)

        wavelengths = np.load(in_wavelengths, allow_pickle=False)
        wavelengths = np.asarray(wavelengths, dtype=np.float64).ravel()
        if wavelengths.size == 0 or not np.all(np.isfinite(wavelengths)) \
                or np.any(wavelengths <= 0):
            raise ValueError("波长数据异常: %s" % in_wavelengths)

        # georef.json 仅作元数据来源（足迹/源 tar/目标 CRS/分辨率），本步不做几何校正
        georef = load_georef_json(in_georef) if os.path.isfile(in_georef) else None
        if georef is None:
            logger.warning("缺少 georef.json，足迹将回退产品 XML；"
                           "不影响正射反射率网格（以栅格自带地理参考为准）")

        source_ds = gdal.Open(in_radiance)
        if source_ds is None:
            raise RuntimeError("无法打开正射辐亮度: %s" % in_radiance)
        if source_ds.RasterCount != wavelengths.size:
            raise ValueError("正射辐亮度波段数 %d != 波长数 %d"
                             % (source_ds.RasterCount, wavelengths.size))
        if georef and georef.get("raster_size"):
            expected = [int(georef["raster_size"][0]), int(georef["raster_size"][1])]
            actual = [source_ds.RasterXSize, source_ds.RasterYSize]
            if expected != actual:
                logger.warning("正射辐亮度尺寸 %s 与 georef.json 记录 %s 不一致",
                               actual, expected)

        # 输出网格完全继承输入正射栅格：读取 CRS/分辨率并做一致性校验（不重采样）
        source_gt = source_ds.GetGeoTransform()
        source_srs = osr.SpatialReference()
        source_srs.ImportFromWkt(source_ds.GetProjection())
        source_epsg = source_srs.GetAuthorityCode(None)
        source_res = abs(source_gt[1]) if source_gt and source_gt[1] else None
        if source_epsg is None or source_res is None:
            logger.warning("正射辐亮度缺少有效 CRS/仿射变换，输出将沿用该空地理参考")
        else:
            logger.info("正射网格: EPSG:%s @ %g m, %d×%d（不做重采样）",
                        source_epsg, source_res,
                        source_ds.RasterXSize, source_ds.RasterYSize)
            if georef and georef.get("target_crs_epsg") is not None:
                try:
                    if int(georef["target_crs_epsg"]) != int(source_epsg):
                        logger.warning("正射辐亮度 CRS (EPSG:%s) 与 georef.json "
                                       "target_crs_epsg=%s 不一致，以栅格为准",
                                       source_epsg, georef["target_crs_epsg"])
                except (TypeError, ValueError):
                    pass
            if georef and georef.get("target_res") is not None:
                try:
                    if abs(float(georef["target_res"]) - source_res) > 1e-6:
                        logger.warning("正射辐亮度分辨率 %g m 与 georef.json "
                                       "target_res=%s 不一致，以栅格为准",
                                       source_res, georef["target_res"])
                except (TypeError, ValueError):
                    pass

        logger.info("正射辐亮度: %d 波段, %d×%d, 波长 %.4f~%.4f μm",
                    source_ds.RasterCount, source_ds.RasterXSize,
                    source_ds.RasterYSize, wavelengths.min(), wavelengths.max())
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        # 1) 读取观测几何、日期、中心位置和足迹
        current_stage = 2
        stage_description = "解析成像时间、观测几何和场景足迹"
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        meta = resolve_metadata(scene, scene_dir, georef, TAR_DIR, logger)
        ring = footprint_ring(meta["footprint"])
        logger.info("观测几何: 太阳 %.2f°/%.2f°, 卫星 %.2f°/%.2f°, 日期 %s",
                    meta["solar_z"], meta["solar_a"], meta["view_z"],
                    meta["view_a"], meta["date"].strftime("%Y-%m-%d"))
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        # 2) 确定 AOT550、海拔、大气模式和气溶胶类型
        current_stage = 3
        stage_description = "确定 AOT550、DEM 海拔及大气/气溶胶模式"
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        aot550, aot_detail = resolve_aot550(meta, ring, args, logger)
        elevation_km, elevation_detail = resolve_elevation(meta, ring, args, logger)
        atmos_profile, atmos_source = resolve_atmos_profile(meta, args, logger)
        aero_profile, aero_source = resolve_aero_profile(meta, args, logger)
        logger.info("6S 参数: 大气=%s, 气溶胶=%s, AOT550=%.4f, 海拔=%.4f km",
                    atmos_profile, aero_profile, aot550, elevation_km)
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        # 3) 逐波段计算 6S 系数
        current_stage = 4
        stage_description = "运行 6S 并计算各波段校正系数"
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        sixs_exe = resolve_sixs_exe(logger)
        sixs = build_sixs(sixs_exe, meta, atmos_profile, aero_profile,
                          aot550, elevation_km, logger)
        coefficients, sixs_version, coef_elapsed, perturbed = \
            compute_sixs_coefficients(sixs, wavelengths, logger)
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        # 4) 逐波段写出反射率，保持输入正射网格且不写入 RPC
        current_stage = 5
        stage_description = "逐波段计算并写出地表反射率"
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        stats = write_reflectance(source_ds, coefficients, wavelengths,
                                  tmp_tif, logger)
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        # 5) 写入成果元数据
        current_stage = 6
        stage_description = "写入影像元数据、参数记录和波段统计"
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        update_ds = gdal.Open(tmp_tif, gdal.GA_Update)
        update_ds.SetMetadataItem("SOURCE_RADIANCE", os.path.basename(in_radiance))
        update_ds.SetMetadataItem("SOURCE_TAR",
                                  os.path.basename(meta["tar"])
                                  if meta["tar"] else "MISSING")
        update_ds.SetMetadataItem("BAND_COUNT", str(source_ds.RasterCount))
        update_ds.SetMetadataItem(
            "IMAGING_DATE", meta["date"].strftime("%Y-%m-%d"))
        update_ds.SetMetadataItem("ATMOSPHERIC_CORRECTION", "6S (Py6S)")
        update_ds.SetMetadataItem("SIXS_EXE", os.path.basename(sixs_exe))
        update_ds.SetMetadataItem("SIXS_VERSION", str(sixs_version or "unknown"))
        update_ds.SetMetadataItem("ATMOS_PROFILE", str(atmos_profile))
        update_ds.SetMetadataItem("AERO_PROFILE", str(aero_profile))
        update_ds.SetMetadataItem("AOT550", "%.6f" % aot550)
        update_ds.SetMetadataItem("ELEVATION_KM", "%.6f" % elevation_km)
        update_ds.SetMetadataItem("AOT_SOURCE", str(aot_detail.get("source", "")))
        update_ds.SetMetadataItem("ELEVATION_SOURCE",
                                  str(elevation_detail.get("source", "")))
        update_ds.SetMetadataItem(
            "REFLECTANCE_FORMULA",
            "rho = (xa*L - xb) / (1 + xc*(xa*L - xb))")
        update_ds.SetMetadataItem(
            "GEOMETRY",
            "orthorectified grid from step1 (same grid, no resampling, no RPC)")
        res_text = ("%g" % source_res) if source_res else "?"
        update_ds.SetMetadataItem(
            "OUTPUT_GRID", "EPSG:%s @ %s m (inherited from %s)"
            % (source_epsg or "?", res_text, IN_RADIANCE_SUFFIX))
        update_ds.SetMetadataItem("PY6S_COEFFICIENTS", OUT_COEF_NAME)
        if georef is not None:
            # 以下仅为步骤1 元数据留档，本步几何完全由栅格自身地理参考决定
            update_ds.SetMetadataItem("GEOREF_PARAMS", IN_GEOREF_NAME)
            update_ds.SetMetadataItem("GEOREF_TARGET_CRS",
                                      "EPSG:%s" % georef.get("target_crs_epsg", "?"))
            update_ds.SetMetadataItem("GEOREF_TARGET_RES",
                                      "%g" % float(georef.get("target_res", 0.0)))
        update_ds = None

        # 6) 原子替换主成果并写出参数和统计
        np.save(tmp_coef, coefficients)
        np.save(tmp_wl, wavelengths)
        os.replace(tmp_tif, final_tif)
        written_targets.append(final_tif)
        os.replace(tmp_coef, final_coef)
        written_targets.append(final_coef)
        os.replace(tmp_wl, final_wl)
        written_targets.append(final_wl)
        georef_local = None
        if georef is not None and os.path.isfile(in_georef):
            copy_file_atomic(in_georef, final_georef)
            written_targets.append(final_georef)
            georef_local = final_georef
        write_band_stats_csv(stats, final_stats)

        payload = {
            "scene": scene,
            "source_radiance": in_radiance,
            "source_tar": meta["tar"],
            "georef_json": in_georef if georef is not None else None,
            "georef_json_local": georef_local,
            "band_count": int(source_ds.RasterCount),
            "raster_size": [int(source_ds.RasterXSize),
                            int(source_ds.RasterYSize)],
            "date": meta["date"].strftime("%Y-%m-%d %H:%M:%S")
            if meta.get("date") else None,
            "satellite": meta.get("satellite"),
            "sensor": meta.get("sensor"),
            "solar_zenith": meta["solar_z"],
            "solar_azimuth": meta["solar_a"],
            "view_zenith": meta["view_z"],
            "view_azimuth": meta["view_a"],
            "center_lat": meta.get("center_lat"),
            "center_lon": meta.get("center_lon"),
            "footprint_wgs84": {key: [float(value[0]), float(value[1])]
                                for key, value in meta["footprint"].items()},
            "atmos_profile": str(atmos_profile),
            "atmos_source": atmos_source,
            "aero_profile": str(aero_profile),
            "aero_source": aero_source,
            "aot550": aot550,
            "aot_source": aot_detail.get("source", ""),
            "aot_detail": aot_detail,
            "elevation_km": elevation_km,
            "elevation_source": elevation_detail.get("source", ""),
            "elevation_detail": elevation_detail,
            "water_vapor": ("water=%.3f g/cm2, ozone=%.3f cm-atm"
                            % (float(WATER_VAPOR_CM),
                               float(OZONE_CM if OZONE_CM is not None
                                     else OZONE_CM_DEFAULT)))
            if WATER_VAPOR_CM is not None else None,
            "sixs_path": sixs_exe,
            "sixs_version": sixs_version,
            "py6s_version": resolve_py6s_version(),
            "sixs_bands": int(len(wavelengths)),
            "sixs_elapsed_s": round(float(coef_elapsed), 2),
            "sixs_perturbed_bands": [
                {"band": int(band),
                 "wavelength_used_um": round(float(value[0]), 6),
                 "offset_nm": round(float(value[1]) * 1000.0, 3)}
                for band, value in sorted(perturbed.items())
            ],
            "formula": "rho = (xa*L - xb) / (1 + xc*(xa*L - xb))",
            "reflectance_clip": list(REFLECTANCE_CLIP) if REFLECTANCE_CLIP else None,
            "nodata": NODATA_OUT,
            # 输出网格直接来自正射辐亮度栅格；georef.json 仅作元数据参考
            "output_grid_crs_epsg": source_epsg,
            "output_grid_res": source_res,
            "output_grid_size": [int(source_ds.RasterXSize),
                                 int(source_ds.RasterYSize)],
            "geometry": "orthorectified grid from step1; no resampling, no RPC",
            "target_crs_epsg": (georef.get("target_crs_epsg")
                                if georef else source_epsg),
            "target_res": (georef.get("target_res") if georef else source_res),
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
        }
        write_parameters(payload, final_params_json, final_params_txt, logger)
        main_complete = True
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        # 7) 生成带坐标预览和质检图
        current_stage = 7
        stage_description = "生成降采样预览和质检图"
        stage_started = log_stage_start(logger, current_stage, 7,
                                        stage_description)
        preview_tif = None
        if MAKE_QC and not args.no_qc:
            if QC_GEOREF_PREVIEW:
                try:
                    prepare_qc_preview(final_tif, wavelengths, final_preview,
                                       logger)
                    preview_tif = final_preview
                except Exception:
                    logger.warning("带坐标预览生成失败，质检图回退无坐标降采样",
                                   exc_info=True)
            try:
                write_quality_figure(final_tif, wavelengths, final_qc, scene,
                                     params_summary_text(payload), logger,
                                     band_stats=stats, preview_tif=preview_tif)
            except Exception:
                logger.warning("质检图生成失败: %s", final_qc, exc_info=True)
        else:
            logger.info("按配置跳过质检图生成")
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)

        elapsed = time.time() - started
        valid_ratio = 100.0 * float(np.mean([item["valid"] for item in stats])) \
            / max(int(source_ds.RasterXSize) * int(source_ds.RasterYSize), 1)
        logger.info("完成: %d 波段，%d×%d，平均有效像元 %.1f%%，"
                    "6S %s，总计 %.1f s，大小 %s",
                    len(wavelengths), source_ds.RasterXSize,
                    source_ds.RasterYSize, valid_ratio,
                    fmt_elapsed(coef_elapsed), elapsed,
                    human_size(os.path.getsize(final_tif)))
        logger.info("  成果目录: %s", out_dir)
        return {"scene": scene, "status": "成功", "bands": int(len(wavelengths)),
                "elapsed": elapsed,
                "message": "%d 波段，AOT=%.3f，海拔=%.3f km，%.1f s"
                           % (len(wavelengths), aot550, elevation_km, elapsed)}
    except Exception:
        for path in (tmp_tif, tmp_coef, tmp_wl):
            try:
                if os.path.exists(path):
                    os.remove(path)
            except OSError:
                pass
        # 本次未完整写出参数记录时，清掉已替换的半成品，避免下次误判为已完成
        if not main_complete:
            for path in written_targets:
                try:
                    if os.path.exists(path):
                        os.remove(path)
                except OSError:
                    pass
        if created_dir:
            try:
                if not os.listdir(out_dir):
                    os.rmdir(out_dir)
            except OSError:
                pass
        logger.error("场景处理在步骤 %d/7 失败：%s", current_stage,
                     stage_description)
        raise
    finally:
        source_ds = None


# =================================================================
# 批处理主流程
# =================================================================

def main(argv=None):
    args = parse_args(argv)
    out_root = os.path.abspath(args.out_root)
    os.makedirs(out_root, exist_ok=True)
    log_path = args.log or os.path.join(out_root, "step2_batch.log")
    logger = setup_logging(log_path, args.verbose)
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)

    in_root = os.path.abspath(args.in_root)
    logger.info("开始扫描步骤1成果目录：%s", in_root)
    try:
        scene_dirs = discover_scene_dirs(in_root)
    except Exception as exc:
        logger.error("扫描步骤1 成果失败: %s", exc)
        return 2

    if args.scene:
        patterns = list(args.scene)
        scene_dirs = [path for path in scene_dirs
                      if any(fnmatch.fnmatch(os.path.basename(path), pattern)
                             for pattern in patterns)]
    if args.limit > 0:
        scene_dirs = scene_dirs[: args.limit]
    if not scene_dirs:
        if args.scene:
            logger.error("在 %s 下没有匹配 --scene %s 的步骤1 成果（%s）",
                         in_root, ", ".join(args.scene),
                         "*" + IN_RADIANCE_SUFFIX)
        else:
            logger.error("在 %s 下未找到可处理的步骤1 成果（%s）",
                         in_root, "*" + IN_RADIANCE_SUFFIX)
        return 2

    logger.info("=" * 72)
    logger.info(
        "ZY1E AHSI 步骤2：6S 大气校正（正射格网，输入 *%s）",
        IN_RADIANCE_SUFFIX)
    logger.info("步骤1 成果: %s", in_root)
    logger.info("输出根目录: %s", out_root)
    logger.info("待处理    : %d 景", len(scene_dirs))
    logger.info("辅助数据  : MODIS=%s, DEM=%s", MODIS_DIR, DEM_DIR)
    logger.info("处理选项  : AOT来源=%s, 覆盖已有成果=%s, 质检图=%s",
                args.aot_source or AOT_SOURCE, args.overwrite,
                MAKE_QC and not args.no_qc)
    logger.info("日志文件  : %s", log_path)
    logger.info("=" * 72)

    if args.dry_run:
        for index, scene_dir in enumerate(scene_dirs, start=1):
            scene = os.path.basename(os.path.normpath(scene_dir))
            logger.info("[%d/%d] %s", index, len(scene_dirs), scene)
            logger.info("        <- %s", resolve_radiance_path(scene_dir))
            logger.info("        -> %s", os.path.join(
                out_root, scene, scene + OUT_REFLECTANCE_SUFFIX))
        logger.info("dry-run 结束，未执行任何处理")
        return 0

    results = []
    started = time.time()
    for index, scene_dir in enumerate(scene_dirs, start=1):
        scene = os.path.basename(os.path.normpath(scene_dir))
        logger.info("-" * 72)
        logger.info("开始处理场景 (%d/%d)：%s", index, len(scene_dirs), scene)
        scene_started = time.time()
        try:
            result = process_scene(scene_dir, out_root, args, logger)
            results.append(result)
            logger.info("场景结束 (%d/%d)：%s，状态=%s，用时 %s",
                        index, len(scene_dirs), scene, result["status"],
                        fmt_elapsed(time.time() - scene_started))
        except Exception as exc:
            logger.exception("处理失败: %s", scene)
            logger.error("场景结束 (%d/%d)：%s，状态=失败，用时 %s",
                         index, len(scene_dirs), scene,
                         fmt_elapsed(time.time() - scene_started))
            results.append({"scene": scene, "status": "失败", "bands": None,
                            "elapsed": 0.0, "message": str(exc)})

    elapsed_total = time.time() - started
    success = [item for item in results if item["status"] == "成功"]
    skipped = [item for item in results if item["status"].startswith("跳过")]
    failed = [item for item in results if item["status"] == "失败"]

    logger.info("=" * 72)
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %s",
                len(success), len(skipped), len(failed),
                fmt_elapsed(elapsed_total))
    for result in results:
        note = result.get("message") or ""
        logger.info("  [%s] %s%s", result["status"], result["scene"],
                    (" - " + note) if note else "")

    if failed:
        logger.error("有 %d 景处理失败，详见日志: %s", len(failed), log_path)
        return 1
    logger.info("全部完成，成果位于: %s", out_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
