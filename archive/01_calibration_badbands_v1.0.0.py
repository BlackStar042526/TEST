# -*- coding: utf-8 -*-
"""
01_calibration_badbands_batch1.0.0.py
ZY1E AHSI 步骤1：辐射定标 + 坏波段剔除 + RPC 正射 + 质检（）

本脚本把预处理链的前三步合并成一个批处理程序，输入 L1A 压缩包，直接输出
**正射到统一分析网格的辐亮度立方体**：

    定标 L = Gain·DN + Bias（166 波段）
        → 水汽/坏波段剔除（24 个 → 142 波段）
        → RPC 正射（EPSG:32650 / 30 m，双线性重采样，落在辐亮度上）
        → radiance_ortho.tif

为什么正射放在大气校正之前
--------------------------
后续步骤要在反射率上做 2.2 μm 附近的窄吸收诊断。若先大气校正、后正射，
几何重采样（双线性/三次卷积）会把反射率上的窄吸收峰“磨钝”；而辐亮度是
线性量，先把重采样做掉、再在大气校正中逐像元处理，对吸收形状的破坏最小。
因此顺序固定为：定标 → 剔坏波段 → 正射 → 大气校正（下一步 02 脚本）。

每景输出（OUT_ROOT/<源文件同名>/）
----------------------------------
    radiance_ortho.tif     正射辐亮度（142 波段, Float32, EPSG:32650, 30 m, nodata=-9999）
    wavelengths_clean.npy  保留波段中心波长（μm，142 个）
    georef.json            本景几何/波段信息（供步骤2/3 读取：网格、足迹、RPC 校正量、
                           剔除波段表、源压缩包等）
    band_stats.csv         逐波段辐射统计（在原始格网上全分辨率统计，值不受重采样影响）
    preview_quality.png    质检成果图（真彩色/假彩色/2.2 μm/均值光谱，带 UTM 坐标）
    radiance_native.tif    原始格网辐亮度（可选，KEEP_NATIVE_RADIANCE=True 时保留；
                           默认处理完即删，只留正射成果）

处理要点
--------
1. 批量：扫描 SRC_DIR 下所有压缩包，逐景处理；单景失败不影响整批。
2. 定标：逐波段流式 L = Gain * DN + Bias，float64 计算、float32 写盘；支持
   Gain/Bias 两列、仅 Gain、自动识别三种定标文件格式。
3. 坏波段：按“原始 166 波段 1-based 编号表”剔除 24 个（见 REMOVE_BAND_GROUPS），
   剩余 142 个；编号表与波长交叉校验，表为空时可退回按 WATER_VAPOR_BANDS 自动剔除。
4. 几何：ZY1E 的 .rpb 与交付影像格网存在约 2.4% 的缩放/平移偏差，直接用原始
   RPC 会偏移约 1 km。本步先用 XML 四角足迹反算并保存校正后的 RPC（残差通常
   <0.01 像元），再以该 RPC 做 gdal.Warp 正射；目标网格 EPSG:32650 / 30 m，
   按四角足迹裁掉黑边。
5. 质检：全分辨率逐波段统计 CSV + 带 UTM 坐标的成果图（matplotlib 可选）。
6. 稳健性：原子写出、跳过已完成、BOM/CRLF/逗号兼容、定标列数自动识别、
   VNIR/SWIR 尺寸与波段数校验、RPC/足迹缺失时明确报错、正射结果网格/波段校验。

用法
----
    python 01_calibration_badbands_batch.py
    python 01_calibration_badbands_batch.py --src "D:\\ZY1E_data" --out "D:\\ZY1E_step1"
    python 01_calibration_badbands_batch.py --dry-run
    python 01_calibration_badbands_batch.py --overwrite --limit 2
    python 01_calibration_badbands_batch.py --keep-native      # 同时保留原始格网辐亮度
    python 01_calibration_badbands_batch.py --no-qc            # 跳过质检图

依赖：GDAL >= 3.0（/vsitar/、RPC、Warp）、numpy；
      质检图可选依赖 matplotlib（缺失时只跳过 PNG，不影响主流程）。
      建议用装好 GDAL 的 Anaconda 解释器，例如 C:\\ProgramData\\anaconda3\\python.exe

后续约定
--------
    步骤2（02_atmospheric_correction_batch.py）：读取 radiance_ortho.tif +
    wavelengths_clean.npy + georef.json，在**同一正射网格**上做 6S 大气校正，
    输出 reflectance_atm.tif（142 波段，EPSG:32650，30 m，无重采样）。
    步骤3：读取步骤2 反射率，按研究区裁剪/镶嵌（不再做 RPC 正射）。
"""

import argparse
import csv
import glob
import json
import logging
import os
import re
import sys
import tarfile
import time

import numpy as np
from osgeo import gdal, osr

gdal.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")

# ╔══════════════════════════════════════════════════════════════╗
# ║  【配置区】按需修改                                            ║
# ╚══════════════════════════════════════════════════════════════╝

# ────────────【原始数据】存放 .tar / .tar.gz 的文件夹 ────────────
SRC_DIR = r"C:\Users\yangz\Desktop\test\ZY1E_data"

# ────────────【输出数据】成果根目录（每景一个同名子文件夹） ────────────
OUT_ROOT = r"C:\Users\yangz\Desktop\test\ZY1E_output_step1"

# 源文件通配符（可用 --pattern 覆盖，支持子目录写法如 "**/*.tar"）
TAR_PATTERNS = ("*.tar", "*.tar.gz", "*.tgz")

# ────────────【参数】ZY1E AHSI 波段结构（166 = 76 VNIR + 90 SWIR） ────────────
N_VNIR = 76
N_SWIR = 90

# ────────────【参数】坏波段剔除表（**原始 166 波段，1-based 编号**） ────────────
# 依据《Al-OH 代码重构规划》阶段 0.2（明舒井论文页 77）：共 24 个 → 剩 142 波段。
# 每组 = (说明, 1-based 波段号元组, 期望起始波长 nm, 期望结束波长 nm)
# 期望波长只用于与光谱响应文件交叉校验，越界只告警、不中断（换传感器时请重填）。
REMOVE_BAND_GROUPS = [
    ("VNIR/SWIR 重叠低信噪", (77, 78, 79), 1006, 1039),
    ("1.4 μm 水汽吸收", (98, 99, 100, 101), 1358, 1408),
    ("1.4 μm 邻近低信噪", (102, 103, 104, 105), 1425, 1476),
    ("1.812 μm 水汽吸收", (125,), 1812, 1812),
    ("1.9 μm 水汽吸收", (126, 127, 128, 129, 130, 131), 1829, 1913),
    ("1.9 μm 邻近低信噪", (132, 133, 134), 1930, 1964),
    ("SWIR 末段低信噪", (164, 165, 166), 2468, 2501),
]

# ────────────【参数】水汽吸收带（μm）────────────────────────────
# 仅用于日志/质检图阴影，以及在 REMOVE_BAND_GROUPS 为空时自动剔除。
WATER_VAPOR_BANDS = [
    (1.35, 1.43),   # 1.4 μm 水汽吸收带
    (1.80, 1.96),   # 1.9 μm 水汽吸收带
]

# ────────────【参数】定标系数格式 ────────────
# True: 每行 "Gain Bias"；False: 每行只有 "Gain"；None: 自动识别
CAL_HAS_BIAS = True

# ────────────【参数】无效像元处理 ────────────
# 默认 False：保持旧版行为（DN=0 定标后为 Bias，不写 -9999）。
# 设为 True 时，DN==NODATA_IN 的像元在成果中写成 NODATA_OUT。
MASK_NODATA = False
NODATA_IN = 0
NODATA_OUT = -9999.0

# ────────────【参数】正射目标网格 ────────────
TARGET_EPSG = 32650        # WGS84 / UTM zone 50N（本区 E117.2，与后续分析网格一致）
TARGET_RES = 30.0          # 目标分辨率（米）
RESAMPLING = "bilinear"    # 辐亮度重采样方法："bilinear" / "cubic" / "cubicspline"
RPC_FIT_MAX_RESIDUAL = 2.0 # 足迹反算 RPC 的允许残差（像元），超过则报错
REQUIRE_GEOREF = True      # True: 缺少 RPC/XML 足迹时报错；False: 仅告警并继续
WRITE_GEOREF = True        # 计算并保存 georef.json（步骤2/3 需要）
KEEP_NATIVE_RADIANCE = False   # True: 额外保留原始格网辐亮度（radiance_native.tif）

# ────────────【参数】批处理行为 ────────────
OVERWRITE = False          # True: 已有成果也重做；False: 跳过
GDAL_CACHE_MB = 1024       # GDAL 块缓存（MB）

# 中间/成果 GeoTIFF 创建选项
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]
# 高光谱立方体用 BSQ（按波段存储），逐波段写入效率高；DEFLATE+PREDICTOR=3 压缩率好
GEOTIFF_COMPRESS = "DEFLATE"
GEOTIFF_COMPRESS_LEVEL = 6

# ────────────【参数】质检成果图 ────────────
MAKE_QC = True             # 生成 preview_quality.png（需要 matplotlib）
QC_MAX_SIZE = 1400         # 质检图采样后长边像素数
QC_DPI = 150

OUT_ORTHO_NAME = "radiance_ortho.tif"
OUT_NATIVE_NAME = "radiance_native.tif"
OUT_WAVELENGTHS_NAME = "wavelengths_clean.npy"
OUT_GEOREF_NAME = "georef.json"
OUT_STATS_NAME = "band_stats.csv"
OUT_QC_NAME = "preview_quality.png"

LOG_NAME = "zy1e.step1"


def _flatten_remove_bands(groups):
    bands = []
    for _, indexes, _, _ in groups:
        bands.extend(int(index) for index in indexes)
    return bands


REMOVE_BANDS_1BASED = _flatten_remove_bands(REMOVE_BAND_GROUPS)


# ╔══════════════════════════════════════════════════════════════╗
# ║  通用工具                                                     ║
# ╚══════════════════════════════════════════════════════════════╝

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="ZY1E AHSI 步骤1：定标 + 水汽/坏波段剔除 + RPC 正射（一次完成）")
    parser.add_argument("--src", default=SRC_DIR,
                        help="存放压缩包的源目录（默认: %(default)s）")
    parser.add_argument("--out", default=OUT_ROOT,
                        help="成果根目录，每景生成同名子文件夹（默认: %(default)s）")
    parser.add_argument("--pattern", action="append", default=None,
                        help="源文件通配符，可多次指定（默认: %s）" % " ".join(TAR_PATTERNS))
    parser.add_argument("--scene", action="append", default=None,
                        help="只处理匹配的场景名，可多次指定，支持通配符")
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE,
                        help="已有成果时覆盖重做（默认跳过）")
    parser.add_argument("--limit", type=int, default=0,
                        help="只处理前 N 景，便于调试（默认 0 = 全部）")
    parser.add_argument("--keep-native", action="store_true", default=KEEP_NATIVE_RADIANCE,
                        help="额外保留原始格网辐亮度 radiance_native.tif（默认只留正射成果）")
    parser.add_argument("--no-qc", action="store_true",
                        help="跳过质检图生成（保留统计 CSV）")
    parser.add_argument("--dry-run", action="store_true",
                        help="只列出待处理场景与目标目录，不执行")
    parser.add_argument("--log", default=None,
                        help="日志文件路径（默认 <out>/step1_batch.log）")
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


def discover_scenes(src_dir, patterns):
    if not os.path.isdir(src_dir):
        raise NotADirectoryError("源目录不存在: %s" % src_dir)
    found = set()
    for pattern in patterns:
        found.update(glob.glob(os.path.join(src_dir, pattern), recursive=True))
    return sorted(f for f in found if os.path.isfile(f))


def scene_name(tar_path):
    name = os.path.basename(tar_path)
    lower = name.lower()
    for ext in (".tar.gz", ".tgz", ".tar"):
        if lower.endswith(ext):
            return name[: -len(ext)]
    return os.path.splitext(name)[0]


def find_tar_member(names, required, suffixes, label, logger):
    hits = []
    for name in names:
        base = os.path.basename(name).lower()
        if not base.endswith(suffixes):
            continue
        if all(key.lower() in base for key in required):
            hits.append(name)
    if not hits:
        raise FileNotFoundError(
            "压缩包中未找到 %s（关键词 %s，扩展名 %s）"
            % (label, "/".join(required) or "任意", "/".join(suffixes)))
    hits.sort()
    if len(hits) > 1:
        logger.warning("%s 匹配到 %d 个文件，使用第一个: %s",
                       label, len(hits), os.path.basename(hits[0]))
    return hits[0]


def decode_text(raw):
    text = raw.decode("utf-8", errors="ignore")
    return text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")


def read_member_text(tar, member, label):
    handle = tar.extractfile(member)
    if handle is None:
        raise FileNotFoundError("%s 无法从压缩包读取: %s" % (label, member))
    with handle:
        return decode_text(handle.read())


# ╔══════════════════════════════════════════════════════════════╗
# ║  产品元数据解析                                               ║
# ╚══════════════════════════════════════════════════════════════╝

def parse_footprint(xml_text, logger):
    """从产品 XML 解析四角经纬度（WGS84），返回 {'TL': (lon, lat), ...}。"""
    footprint = {}
    corner_map = {"TL": "TopLeft", "TR": "TopRight",
                  "BR": "BottomRight", "BL": "BottomLeft"}
    for abbr, corner in corner_map.items():
        lat_match = re.search(r"<%sLatitude>\s*([-+0-9.eE]+)\s*</%sLatitude>"
                              % (corner, corner), xml_text, re.IGNORECASE)
        lon_match = re.search(r"<%sLongitude>\s*([-+0-9.eE]+)\s*</%sLongitude>"
                              % (corner, corner), xml_text, re.IGNORECASE)
        if not lat_match or not lon_match:
            raise ValueError("产品 XML 缺少 %s 经纬度" % corner)
        footprint[abbr] = (float(lon_match.group(1)), float(lat_match.group(1)))
    return footprint


def parse_wavelengths(text, n_bands, label, logger):
    """从 Spectralresponse.raw 提取中心波长（每行首列为 nm），返回 μm。"""
    values = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        try:
            values.append(float(line.replace(",", " ").split()[0]))
        except (ValueError, IndexError):
            continue
    if len(values) < n_bands:
        raise ValueError(
            "%s Spectralresponse 只有 %d 行有效数据，少于 %d 个波段"
            % (label, len(values), n_bands))
    if len(values) > n_bands:
        logger.debug("%s Spectralresponse 有 %d 行，截取前 %d 行",
                     label, len(values), n_bands)
    wavelengths = np.asarray(values[:n_bands], dtype=np.float64) / 1000.0
    if not np.all(np.isfinite(wavelengths)) or np.any(wavelengths <= 0):
        raise ValueError("%s 波长数据异常（存在非正数或 NaN）" % label)
    logger.info("%s 波长: %d 波段, %.4f ~ %.4f μm",
                label, wavelengths.size, wavelengths.min(), wavelengths.max())
    return wavelengths


def parse_calibration(text, n_bands, has_bias, label, logger):
    """解析定标系数，返回 (gain, bias) 两个 float64 数组。"""
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = re.split(r"[,\s]+", line)
        numbers = []
        for token in tokens:
            try:
                numbers.append(float(token))
            except ValueError:
                break
        if numbers:
            rows.append(numbers[:2])
    if not rows:
        raise ValueError("%s 定标系数文件为空" % label)
    if has_bias is None:
        has_bias = all(len(row) >= 2 for row in rows)
        logger.info("%s 自动识别定标格式: %s", label,
                    "Gain + Bias" if has_bias else "仅 Gain")
    if has_bias:
        bad_rows = [i for i, row in enumerate(rows) if len(row) < 2]
        if bad_rows:
            raise ValueError(
                "%s 定标文件第 %d 行缺少 Bias 列，请检查 CAL_HAS_BIAS"
                % (label, bad_rows[0] + 1))
        data = np.asarray([row[:2] for row in rows], dtype=np.float64)
    else:
        data = np.asarray([[row[0], 0.0] for row in rows], dtype=np.float64)
    if len(data) < n_bands:
        raise ValueError(
            "%s 定标系数只有 %d 行，少于 %d 个波段；请检查 CAL_HAS_BIAS 或文件格式"
            % (label, len(data), n_bands))
    if len(data) > n_bands:
        logger.debug("%s 定标系数 %d 行，截取前 %d 行", label, len(data), n_bands)
    gain = data[:n_bands, 0]
    bias = data[:n_bands, 1]
    if not np.all(np.isfinite(gain)) or not np.all(np.isfinite(bias)):
        raise ValueError("%s 定标系数包含 NaN/Inf" % label)
    if np.any(gain == 0):
        logger.warning("%s 有 %d 个波段 Gain = 0", label, int(np.sum(gain == 0)))
    logger.info("%s 定标: Gain %.4f ~ %.4f, Bias %.4f ~ %.4f",
                label, gain.min(), gain.max(), bias.min(), bias.max())
    return gain, bias


# ╔══════════════════════════════════════════════════════════════╗
# ║  坏波段剔除                                                   ║
# ╚══════════════════════════════════════════════════════════════╝

def detect_water_vapor_bands(wavelengths, band_ranges, logger):
    """按波长定位水汽带，返回 (连续区间列表, 剔除的 0-based 索引列表)。"""
    removed = [idx for idx, wl in enumerate(wavelengths)
               if any(lo <= wl <= hi for lo, hi in band_ranges)]
    ranges = []
    if removed:
        start = prev = removed[0]
        for idx in removed[1:]:
            if idx == prev + 1:
                prev = idx
            else:
                ranges.append((start, prev))
                start = prev = idx
        ranges.append((start, prev))
    logger.info("水汽带自动识别: 剔除 %d 个波段，合并为 %d 段",
                len(removed), len(ranges))
    for start, end in ranges:
        logger.info("    波段 %3d-%3d  (%.4f ~ %.4f μm)",
                    start, end, wavelengths[start], wavelengths[end])
    return ranges, removed


def resolve_band_selection(wavelengths, logger):
    """
    决定剔除哪些波段，返回 (keep 索引数组, removed 0-based 列表, 剔除组日志信息)。

    优先使用 REMOVE_BAND_GROUPS（显式 1-based 编号表）；表为空时按
    WATER_VAPOR_BANDS 自动剔除。
    """
    n_total = int(wavelengths.size)
    if REMOVE_BANDS_1BASED:
        invalid = [b for b in REMOVE_BANDS_1BASED if not (1 <= b <= n_total)]
        if invalid:
            raise ValueError(
                "REMOVE_BAND_GROUPS 中的波段号超出 1~%d 范围: %s"
                % (n_total, invalid[:10]))
        duplicates = sorted({b for b in REMOVE_BANDS_1BASED
                             if REMOVE_BANDS_1BASED.count(b) > 1})
        if duplicates:
            logger.warning("坏波段表存在重复编号，已去重: %s", duplicates)
        removed = sorted({int(b) - 1 for b in REMOVE_BANDS_1BASED})
        groups = []
        for label, indexes, wl_lo, wl_hi in REMOVE_BAND_GROUPS:
            idx0 = sorted({int(b) - 1 for b in indexes})
            if not idx0:
                continue
            actual_lo = float(wavelengths[idx0[0]]) * 1000.0
            actual_hi = float(wavelengths[idx0[-1]]) * 1000.0
            groups.append({
                "label": label,
                "bands_1based": [int(b) for b in indexes],
                "wavelength_nm": [actual_lo, actual_hi],
            })
            if actual_hi < wl_lo - 6.0 or actual_lo > wl_hi + 6.0:
                logger.warning(
                    "坏波段组「%s」编号与期望波长不符: 期望 %d~%d nm，"
                    "实际 %.1f~%.1f nm；请核对波段表",
                    label, wl_lo, wl_hi, actual_lo, actual_hi)
            logger.info("  剔除「%s」: 波段 %s (%.1f~%.1f nm)",
                        label, ",".join(str(b) for b in indexes),
                        actual_lo, actual_hi)
    else:
        removed = detect_water_vapor_bands(wavelengths, WATER_VAPOR_BANDS, logger)[1]
        groups = [{"label": "按 WATER_VAPOR_BANDS 自动剔除",
                   "bands_1based": [i + 1 for i in removed],
                   "wavelength_nm": [float(wavelengths[removed[0]]) * 1000.0,
                                     float(wavelengths[removed[-1]]) * 1000.0]}
                  if removed else {"label": "未剔除", "bands_1based": [],
                                   "wavelength_nm": []}]

    removed_set = set(removed)
    keep = np.asarray([i for i in range(n_total) if i not in removed_set],
                      dtype=np.int64)
    if keep.size <= 0:
        raise ValueError("剔除后没有剩余波段，请检查 REMOVE_BAND_GROUPS")
    logger.info("坏波段剔除: %d → %d 波段（剔除 %d 个）",
                n_total, keep.size, len(removed))
    return keep, removed, groups


# ╔══════════════════════════════════════════════════════════════╗
# ║  几何：RPC 格网校正 + 正射裁剪线                                ║
# ╚══════════════════════════════════════════════════════════════╝

def open_vsitar_dataset(tar_path, internal_path, label):
    """通过 /vsitar/ 直接打开压缩包内的影像（不解压、不落盘）。"""
    vsi_path = "/vsitar/" + os.path.abspath(tar_path).replace("\\", "/") + "/" + internal_path
    dataset = gdal.Open(vsi_path)
    if dataset is None:
        raise RuntimeError(
            "无法通过 /vsitar/ 打开 %s: %s\n"
            "请确认 GDAL 版本支持 /vsitar/（GDAL >= 3.0）" % (label, vsi_path))
    return dataset


def compute_rpc_grid_correction(ds_vn, footprint, logger):
    """
    用官方四角足迹反算 RPC 格网校正量。

    ZY1E 交付的 L1A tif 与 .rpb 的格网并不一致（交付影像经过重采样），
    实测差异是一个纯缩放+平移：raw_pixel = scale * delivered_pixel + shift。
    这里通过 RPC 反变换把四角经纬度换算成 RPC 格网像素，再最小二乘拟合
    scale/shift，得到可写回 RPC 的修正参数（只改 LINE_OFF/SAMP_OFF/SCALE）。

    返回 (修正后的 RPC 元数据字典, 校正信息字典)。
    """
    rpc = ds_vn.GetMetadata("RPC")
    if not rpc:
        raise ValueError("影像缺少 RPC 元数据（.rpb），无法生成地理参考参数")
    transformer = gdal.Transformer(ds_vn, None, ["METHOD=RPC"])
    width, height = ds_vn.RasterXSize, ds_vn.RasterYSize
    pixel_corners = np.array([[0.0, 0.0], [width, 0.0],
                              [width, height], [0.0, height]])
    raw = []
    for corner in ("TL", "TR", "BR", "BL"):
        lon, lat = footprint[corner]
        ok, point = transformer.TransformPoint(1, lon, lat)
        if not ok:
            raise ValueError("RPC 反变换失败: 足迹角点 %s (%.6f, %.6f)"
                             % (corner, lon, lat))
        raw.append([point[0], point[1]])
    raw = np.asarray(raw, dtype=np.float64)

    (sx, tx) = np.linalg.lstsq(
        np.column_stack([pixel_corners[:, 0], np.ones(4)]), raw[:, 0],
        rcond=None)[0]
    (sy, ty) = np.linalg.lstsq(
        np.column_stack([pixel_corners[:, 1], np.ones(4)]), raw[:, 1],
        rcond=None)[0]
    predicted = np.column_stack([pixel_corners[:, 0] * sx + tx,
                                 pixel_corners[:, 1] * sy + ty])
    residual = float(np.max(np.hypot(predicted[:, 0] - raw[:, 0],
                                     predicted[:, 1] - raw[:, 1])))
    logger.info("RPC 格网校正: sx=%.6f tx=%.1f, sy=%.6f ty=%.1f, 拟合残差 %.4f 像元",
                sx, tx, sy, ty, residual)
    if residual > RPC_FIT_MAX_RESIDUAL:
        raise ValueError(
            "足迹与 RPC 格网不一致（残差 %.2f 像元 > %.1f），"
            "请检查 .rpb / XML 是否匹配" % (residual, RPC_FIT_MAX_RESIDUAL))

    corrected = dict(rpc)
    corrected["SAMP_OFF"] = str((float(rpc["SAMP_OFF"]) - tx) / sx)
    corrected["SAMP_SCALE"] = str(float(rpc["SAMP_SCALE"]) / sx)
    corrected["LINE_OFF"] = str((float(rpc["LINE_OFF"]) - ty) / sy)
    corrected["LINE_SCALE"] = str(float(rpc["LINE_SCALE"]) / sy)
    info = {"sx": sx, "tx": tx, "sy": sy, "ty": ty, "residual": residual}
    return corrected, info


def write_cutline_vsimem(footprint, logger):
    """用四角足迹构造 WGS84 多边形裁剪线，写入 /vsimem。"""
    ring = [footprint[key] for key in ("TL", "TR", "BR", "BL")]
    ring.append(ring[0])
    geojson = {
        "type": "FeatureCollection",
        "features": [{
            "type": "Feature",
            "properties": {},
            "geometry": {
                "type": "Polygon",
                "coordinates": [[[float(lon), float(lat)] for lon, lat in ring]],
            },
        }],
    }
    path = "/vsimem/zy1e_step1_footprint.geojson"
    gdal.FileFromMemBuffer(path, json.dumps(geojson).encode("utf-8"))
    logger.debug("裁剪线已创建: %s", path)
    return path


def build_gtiff_creation_options(band_interleave=True):
    creation = list(GTIFF_OPTIONS)
    if band_interleave:
        creation.append("INTERLEAVE=BAND")
    if GEOTIFF_COMPRESS:
        creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS,
                         "PREDICTOR=3",
                         "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL])
    return creation


def orthorectify_radiance(native_path, footprint, ortho_tmp, logger):
    """
    把原始格网辐亮度（带校正后 RPC）RPC 正射到 TARGET_EPSG/TARGET_RES，
    按足迹裁掉黑边。返回正射成果的网格信息 dict。
    """
    cutline = write_cutline_vsimem(footprint, logger)
    try:
        options = gdal.WarpOptions(
            format="GTiff",
            dstSRS="EPSG:%d" % TARGET_EPSG,
            xRes=TARGET_RES,
            yRes=TARGET_RES,
            resampleAlg=RESAMPLING,
            dstNodata=NODATA_OUT,
            targetAlignedPixels=True,
            multithread=True,
            warpOptions=["NUM_THREADS=ALL_CPUS"],
            creationOptions=build_gtiff_creation_options(),
            rpc=True,
            cutlineDSName=cutline,
            cropToCutline=True,
            outputType=gdal.GDT_Float32,
        )
        if os.path.exists(ortho_tmp):
            os.remove(ortho_tmp)
        logger.info("RPC 正射: EPSG:%d @ %g m（%s 重采样，按足迹裁剪）...",
                    TARGET_EPSG, TARGET_RES, RESAMPLING)
        result = gdal.Warp(ortho_tmp, native_path, options=options)
        if result is None:
            raise RuntimeError("RPC 正射失败: %s" % ortho_tmp)
        width, height = result.RasterXSize, result.RasterYSize
        geo_transform = result.GetGeoTransform()
        projection = result.GetProjection()
        band_count = result.RasterCount
        result = None
    finally:
        try:
            gdal.Unlink(cutline)
        except Exception:
            pass

    if width <= 0 or height <= 0:
        raise RuntimeError("正射结果尺寸异常: %dx%d" % (width, height))
    x_res = abs(geo_transform[1])
    y_res = abs(geo_transform[5])
    if abs(x_res - TARGET_RES) > 1e-3 or abs(y_res - TARGET_RES) > 1e-3:
        raise RuntimeError("正射结果分辨率异常: %.6f x %.6f" % (x_res, y_res))
    srs = osr.SpatialReference()
    srs.ImportFromWkt(projection)
    epsg = srs.GetAuthorityCode(None)
    if epsg is None or int(epsg) != TARGET_EPSG:
        logger.warning("正射结果 CRS 非 EPSG:%d（实际 %s）", TARGET_EPSG, epsg)
    bounds = [geo_transform[0],
              geo_transform[3] + geo_transform[5] * height,
              geo_transform[0] + geo_transform[1] * width,
              geo_transform[3]]
    logger.info("正射成果: %d×%d×%d 波段，%.4f×%.4f m，范围 %.0f,%.0f ~ %.0f,%.0f",
                width, height, band_count, x_res, y_res, *bounds)
    return {
        "width": int(width),
        "height": int(height),
        "band_count": int(band_count),
        "geotransform": [float(value) for value in geo_transform],
        "projection_wkt": projection,
        "crs_epsg": int(epsg) if epsg is not None and epsg.isdigit() else TARGET_EPSG,
        "bounds": [float(value) for value in bounds],
    }


# ╔══════════════════════════════════════════════════════════════╗
# ║  元数据 / 统计输出                                             ║
# ╚══════════════════════════════════════════════════════════════╝

def write_georef_json(output_path, scene, tar_path, vn_tif, sw_tif,
                      native_size, ortho_info, band_count, bands_source,
                      footprint, rpc_correction, wavelengths, keep, removed,
                      remove_groups, logger):
    """把步骤2/3 所需的几何、网格与波段信息写成 JSON。"""
    created = time.strftime("%Y-%m-%d %H:%M:%S")
    payload = {
        "scene": scene,
        "source_tar": os.path.basename(tar_path),
        "source_tiffs": [os.path.basename(vn_tif), os.path.basename(sw_tif)],
        "created": created,
        "geometry": "RPC orthorectified (corrected RPC) → EPSG:%d @ %g m"
                    % (TARGET_EPSG, TARGET_RES),
        "target_crs_epsg": TARGET_EPSG,
        "target_res": TARGET_RES,
        "resampling": RESAMPLING,
        "native_raster_size": [int(native_size[0]), int(native_size[1])],
        "raster_size": [int(ortho_info["width"]), int(ortho_info["height"])],
        "band_count": int(band_count),
        "bands_source": int(bands_source),
        "geotransform": [float(value) for value in ortho_info["geotransform"]],
        "bounds": [float(value) for value in ortho_info["bounds"]],
        "crs_wkt": ortho_info["projection_wkt"],
        "footprint_wgs84": {key: [float(value[0]), float(value[1])]
                            for key, value in footprint.items()},
        "rpc_grid_correction": {
            "sx": rpc_correction["sx"],
            "tx": rpc_correction["tx"],
            "sy": rpc_correction["sy"],
            "ty": rpc_correction["ty"],
            "residual_px": rpc_correction["residual"],
        },
        "remove_band_groups": remove_groups,
        "removed_bands_1based": [int(index) + 1 for index in removed],
        "removed_bands_0based": [int(index) for index in removed],
        "removed_wavelengths_um": [float(wavelengths[index]) for index in removed],
        "water_vapor_bands_um": [[float(lo), float(hi)]
                                 for lo, hi in WATER_VAPOR_BANDS],
        "radiance_formula": "L = Gain * DN + Bias",
        "nodata": NODATA_OUT,
        "mask_nodata": bool(MASK_NODATA),
        "dem_used": False,
        "note": "本景已正射；步骤2 直接在同一网格上做大气校正，无需再次 warp",
    }
    tmp_path = output_path + ".tmp"
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2)
    os.replace(tmp_path, output_path)
    logger.info("  几何/波段参数: %s", output_path)


def write_band_stats_csv(stats, output_path):
    with open(output_path, "w", newline="", encoding="utf-8") as handle:
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


# ╔══════════════════════════════════════════════════════════════╗
# ║  质检：成果图                                                 ║
# ╚══════════════════════════════════════════════════════════════╝

def _nearest_band(wavelengths, target):
    return int(np.argmin(np.abs(wavelengths - target)))


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


def _stretch(array, valid):
    values = array[valid]
    if values.size == 0:
        return np.zeros_like(array)
    low, high = np.percentile(values, [2, 98])
    if high <= low:
        high = low + 1e-6
    return np.clip((array - low) / (high - low), 0.0, 1.0)


def load_qc_arrays(tif_path, band_map, unique):
    """从正射成果读取 5 个质检波段（按 QC_MAX_SIZE 降采样），返回数组与地理范围。"""
    dataset = gdal.Open(tif_path)
    if dataset is None:
        return None
    width, height = dataset.RasterXSize, dataset.RasterYSize
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))
    buf_w, buf_h = max(1, width // scale), max(1, height // scale)
    arrays = {}
    for name, band_index in band_map.items():
        position = unique.index(band_index) + 1
        arrays[name] = dataset.GetRasterBand(position).ReadAsArray(
            buf_xsize=buf_w, buf_ysize=buf_h).astype(np.float64)
    geo_transform = dataset.GetGeoTransform()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(dataset.GetProjection())
    epsg = srs.GetAuthorityCode(None) or "?"
    dataset = None
    extent = [geo_transform[0] / 1000.0,
              (geo_transform[0] + geo_transform[1] * width) / 1000.0,
              (geo_transform[3] + geo_transform[5] * height) / 1000.0,
              geo_transform[3] / 1000.0]
    info = {
        "extent": extent,
        "xlabel": "Easting (km)",
        "ylabel": "Northing (km)",
        "crs_text": "EPSG:%s @ %g m ortho" % (epsg, abs(geo_transform[1])),
        "size_text": "%d×%d" % (width, height),
    }
    return arrays, info


def write_quality_figure(tif_path, wavelengths, removed_ranges, output_png,
                         scene, logger, band_stats=None):
    """生成质检成果图：真彩色 / 假彩色 / 短波红外 / 均值光谱（带 UTM 坐标）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        logger.warning("matplotlib 不可用，跳过质检图: %s", exc)
        return False

    wavelengths = np.asarray(wavelengths, dtype=np.float64)
    band_map, unique = select_qc_bands(wavelengths)
    prepared = load_qc_arrays(tif_path, band_map, unique)
    if prepared is None:
        logger.warning("无法读取正射成果，跳过质检图: %s", tif_path)
        return False
    arrays, info = prepared
    extent = info["extent"]

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
    swir = arrays["SWIR"].copy()
    swir_valid = np.isfinite(swir) & (swir != NODATA_OUT)
    swir_stretched = _stretch(swir, swir_valid)
    swir_stretched[~swir_valid] = np.nan

    wl_r = wavelengths[band_map["R"]]
    wl_g = wavelengths[band_map["G"]]
    wl_b = wavelengths[band_map["B"]]
    wl_nir = wavelengths[band_map["NIR"]]
    wl_swir = wavelengths[band_map["SWIR"]]

    if band_stats:
        figure, axes = plt.subplots(2, 2, figsize=(13, 10))
        ax_true, ax_false, ax_swir, ax_spec = axes.ravel()
    else:
        figure, axes = plt.subplots(1, 3, figsize=(15, 5.5))
        ax_true, ax_false, ax_swir = axes.ravel()
        ax_spec = None

    ax_true.imshow(true_rgb, extent=extent, origin="upper", interpolation="nearest")
    ax_true.set_title("True color  (R %.2f / G %.2f / B %.2f μm)"
                      % (wl_r, wl_g, wl_b), fontsize=10)
    ax_false.imshow(false_rgb, extent=extent, origin="upper", interpolation="nearest")
    ax_false.set_title("False color  (NIR %.2f / R %.2f / G %.2f μm)"
                       % (wl_nir, wl_r, wl_g), fontsize=10)
    ax_swir.imshow(swir_stretched, extent=extent, origin="upper",
                   cmap="gray", interpolation="nearest")
    ax_swir.set_title("SWIR %.2f μm (band %d, 2–98%% stretch)"
                      % (wl_swir, band_map["SWIR"] + 1), fontsize=10)
    for axis in (ax_true, ax_false, ax_swir):
        axis.set_xlabel(info["xlabel"], fontsize=9)
        axis.set_ylabel(info["ylabel"], fontsize=9)
        axis.tick_params(labelsize=8)

    if ax_spec is not None:
        spectrum_wl = np.array([item["wavelength"] for item in band_stats])
        spectrum_mean = np.array([item["mean"] for item in band_stats])
        spectrum_std = np.array([item["std"] for item in band_stats])
        ax_spec.plot(spectrum_wl, spectrum_mean, color="black", linewidth=1.0)
        ax_spec.fill_between(spectrum_wl,
                             spectrum_mean - spectrum_std,
                             spectrum_mean + spectrum_std,
                             color="0.75", alpha=0.4, linewidth=0)
        for low, high in removed_ranges:
            ax_spec.axvspan(low, high, color="0.88", zorder=0)
        ax_spec.set_xlabel("Wavelength (μm)", fontsize=9)
        ax_spec.set_ylabel("Mean radiance", fontsize=9)
        ax_spec.grid(alpha=0.3)
        ax_spec.set_title("Mean spectrum (native grid, shaded = removed bands)",
                          fontsize=10)

    valid_ratio = float(valid.mean()) if valid.size else 0.0
    figure.suptitle(
        "%s | %d bands | %s | %s | valid %.1f%%"
        % (scene, len(wavelengths), info["size_text"], info["crs_text"],
           100.0 * valid_ratio),
        fontsize=11)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    figure.savefig(output_png, dpi=QC_DPI, bbox_inches="tight")
    plt.close(figure)
    logger.info("质检成果图: %s", output_png)
    return True


# ╔══════════════════════════════════════════════════════════════╗
# ║  单景处理                                                     ║
# ╚══════════════════════════════════════════════════════════════╝

def resolve_scene_paths(scene_dir):
    return {
        "ortho": os.path.join(scene_dir, OUT_ORTHO_NAME),
        "native": os.path.join(scene_dir, OUT_NATIVE_NAME),
        "wavelengths": os.path.join(scene_dir, OUT_WAVELENGTHS_NAME),
        "georef": os.path.join(scene_dir, OUT_GEOREF_NAME),
        "stats": os.path.join(scene_dir, OUT_STATS_NAME),
        "qc": os.path.join(scene_dir, OUT_QC_NAME),
    }


def regenerate_qc(scene_dir, scene, logger):
    """已有主干成果时，仅补齐质检图。"""
    paths = resolve_scene_paths(scene_dir)
    try:
        wavelengths = np.load(paths["wavelengths"])
        stats = []
        if os.path.isfile(paths["stats"]):
            with open(paths["stats"], "r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    stats.append({
                        "band": int(row["band"]),
                        "wavelength": float(row["wavelength_um"]),
                        "mean": float(row["mean"]),
                        "std": float(row["std"]),
                    })
        write_quality_figure(paths["ortho"], wavelengths, WATER_VAPOR_BANDS,
                             paths["qc"], scene, logger,
                             band_stats=stats or None)
    except Exception:
        logger.warning("补齐质检图失败: %s", scene_dir, exc_info=True)


def process_scene(tar_path, out_root, args, logger):
    """处理一景：定标 + 坏波段剔除 + RPC 正射 + 质检。"""
    scene = scene_name(tar_path)
    scene_dir = os.path.join(out_root, scene)
    paths = resolve_scene_paths(scene_dir)
    make_qc = MAKE_QC and not args.no_qc

    have_main = all(os.path.isfile(path) for path in
                    (paths["ortho"], paths["wavelengths"], paths["georef"]))
    if not args.overwrite and have_main:
        if make_qc and not os.path.isfile(paths["qc"]):
            regenerate_qc(scene_dir, scene, logger)
        logger.info("跳过（成果已存在）: %s", scene_dir)
        return {"scene": scene, "status": "跳过", "bands": None,
                "elapsed": 0.0, "message": "成果已存在"}

    started = time.time()
    created_dir = not os.path.isdir(scene_dir)
    os.makedirs(scene_dir, exist_ok=True)
    tmp_native = os.path.join(scene_dir, "radiance_native.tmp.tif")
    tmp_ortho = os.path.join(scene_dir, "radiance_ortho.tmp.tif")
    tmp_wl = os.path.join(scene_dir, "wavelengths_clean.tmp.npy")
    ds_vn = ds_swir = None
    try:
        with tarfile.open(tar_path, "r:*") as tar:
            names = tar.getnames()
            logger.info("压缩包内文件数: %d", len(names))

            vn_tif = find_tar_member(names, ("_vn",), (".tif", ".tiff"), "VNIR 影像", logger)
            sw_tif = find_tar_member(names, ("_sw",), (".tif", ".tiff"), "SWIR 影像", logger)
            vn_cal = find_tar_member(names, ("vnir", "radcal"), (".raw",), "VNIR 定标系数", logger)
            sw_cal = find_tar_member(names, ("swir", "radcal"), (".raw",), "SWIR 定标系数", logger)
            vn_sr = find_tar_member(names, ("vnir", "spectralresponse"), (".raw",), "VNIR 光谱响应", logger)
            sw_sr = find_tar_member(names, ("swir", "spectralresponse"), (".raw",), "SWIR 光谱响应", logger)
            xml_member = find_tar_member(names, (), (".xml",), "产品元数据 XML", logger)

            logger.info("  VNIR 影像: %s", os.path.basename(vn_tif))
            logger.info("  SWIR 影像: %s", os.path.basename(sw_tif))
            logger.info("  VNIR 定标: %s", os.path.basename(vn_cal))
            logger.info("  SWIR 定标: %s", os.path.basename(sw_cal))
            logger.info("  VNIR 光谱响应: %s", os.path.basename(vn_sr))
            logger.info("  SWIR 光谱响应: %s", os.path.basename(sw_sr))
            logger.info("  元数据 XML: %s", os.path.basename(xml_member))

            wavelengths = np.concatenate([
                parse_wavelengths(read_member_text(tar, vn_sr, "VNIR 光谱响应"),
                                  N_VNIR, "VNIR", logger),
                parse_wavelengths(read_member_text(tar, sw_sr, "SWIR 光谱响应"),
                                  N_SWIR, "SWIR", logger),
            ])
            if wavelengths.size != N_VNIR + N_SWIR:
                raise ValueError("波长数 %d != %d"
                                 % (wavelengths.size, N_VNIR + N_SWIR))
            if np.any(np.diff(wavelengths) <= 0):
                logger.warning(
                    "VNIR 与 SWIR 波长存在重叠（VNIR 末波段 %.4f μm ≥ SWIR 首波段 %.4f μm），"
                    "波段顺序按 VNIR→SWIR 拼接保留",
                    wavelengths[N_VNIR - 1], wavelengths[N_VNIR])

            vn_gain, vn_bias = parse_calibration(
                read_member_text(tar, vn_cal, "VNIR 定标系数"),
                N_VNIR, CAL_HAS_BIAS, "VNIR", logger)
            sw_gain, sw_bias = parse_calibration(
                read_member_text(tar, sw_cal, "SWIR 定标系数"),
                N_SWIR, CAL_HAS_BIAS, "SWIR", logger)

            keep, removed, remove_groups = resolve_band_selection(wavelengths, logger)

            footprint = None
            if WRITE_GEOREF:
                try:
                    footprint = parse_footprint(
                        read_member_text(tar, xml_member, "产品元数据 XML"), logger)
                    logger.info("  足迹 TL (%.6f, %.6f)  BR (%.6f, %.6f)",
                                footprint["TL"][0], footprint["TL"][1],
                                footprint["BR"][0], footprint["BR"][1])
                except Exception:
                    if REQUIRE_GEOREF:
                        raise
                    logger.warning("未解析到足迹，地理参考参数将缺失", exc_info=True)

        ds_vn = open_vsitar_dataset(tar_path, vn_tif, "VNIR 影像")
        ds_swir = open_vsitar_dataset(tar_path, sw_tif, "SWIR 影像")

        if ds_vn.RasterCount != N_VNIR:
            raise ValueError("VNIR 影像波段数 %d != %d" % (ds_vn.RasterCount, N_VNIR))
        if ds_swir.RasterCount != N_SWIR:
            raise ValueError("SWIR 影像波段数 %d != %d" % (ds_swir.RasterCount, N_SWIR))
        if (ds_vn.RasterXSize, ds_vn.RasterYSize) != (ds_swir.RasterXSize, ds_swir.RasterYSize):
            raise ValueError(
                "VNIR/SWIR 影像尺寸不一致: VNIR %dx%d, SWIR %dx%d"
                % (ds_vn.RasterXSize, ds_vn.RasterYSize,
                   ds_swir.RasterXSize, ds_swir.RasterYSize))

        geo_transform = ds_vn.GetGeoTransform()
        projection = ds_vn.GetProjection()
        if not np.allclose(geo_transform, ds_swir.GetGeoTransform()):
            logger.warning("VNIR/SWIR GeoTransform 不一致，以 VNIR 为准")
        if projection and ds_swir.GetProjection() and projection != ds_swir.GetProjection():
            logger.warning("VNIR/SWIR 投影不一致，以 VNIR 为准")

        rpc_corrected = None
        rpc_correction = None
        if WRITE_GEOREF and footprint is not None:
            try:
                rpc_corrected, rpc_correction = compute_rpc_grid_correction(
                    ds_vn, footprint, logger)
            except Exception:
                if REQUIRE_GEOREF:
                    raise
                logger.warning("RPC 校正失败，正射将无法进行", exc_info=True)
        if rpc_corrected is None:
            raise ValueError("缺少可用的校正后 RPC，无法正射；请检查 .rpb / XML 足迹")

        n_rows, n_cols = ds_vn.RasterYSize, ds_vn.RasterXSize
        n_original = N_VNIR + N_SWIR
        n_output = int(keep.size)

        band_plan = []
        for global_idx in keep:
            if global_idx < N_VNIR:
                band_plan.append((ds_vn, int(global_idx) + 1,
                                  vn_gain[global_idx], vn_bias[global_idx],
                                  wavelengths[global_idx]))
            else:
                local = int(global_idx) - N_VNIR
                band_plan.append((ds_swir, local + 1,
                                  sw_gain[local], sw_bias[local],
                                  wavelengths[global_idx]))

        # ── 1) 定标 + 坏波段剔除：原始格网逐波段写出（带校正后 RPC） ──
        driver = gdal.GetDriverByName("GTiff")
        creation = build_gtiff_creation_options()
        if os.path.exists(tmp_native):
            os.remove(tmp_native)
        out_ds = driver.Create(tmp_native, n_cols, n_rows, n_output,
                               gdal.GDT_Float32, options=creation)
        if out_ds is None:
            raise RuntimeError("无法创建 GeoTIFF: %s" % tmp_native)
        out_ds.SetGeoTransform(geo_transform)
        out_ds.SetProjection(projection)
        out_ds.SetMetadata(rpc_corrected, "RPC")

        stats = []
        masked_pixels = 0
        for out_idx, (src_ds, src_band, gain, bias, wavelength) in enumerate(
                band_plan, start=1):
            dn = src_ds.GetRasterBand(src_band).ReadAsArray()
            radiance = dn.astype(np.float64) * float(gain) + float(bias)
            if MASK_NODATA:
                invalid = dn == NODATA_IN
                count = int(invalid.sum())
                if count:
                    radiance[invalid] = NODATA_OUT
                    masked_pixels += count
            out_band = out_ds.GetRasterBand(out_idx)
            out_band.WriteArray(radiance.astype(np.float32))
            out_band.SetNoDataValue(NODATA_OUT)
            out_band.SetDescription("%.4f um" % wavelength)

            valid_values = radiance[radiance != NODATA_OUT] if MASK_NODATA \
                else radiance[np.isfinite(radiance)]
            if valid_values.size:
                sample = valid_values[::4]
                p1, p99 = np.percentile(sample, [1, 99])
                stats.append({
                    "band": out_idx,
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
                stats.append({"band": out_idx, "wavelength": float(wavelength),
                              "min": 0.0, "max": 0.0, "mean": 0.0, "std": 0.0,
                              "p1": 0.0, "p99": 0.0, "valid": 0})

            if out_idx == 1 or out_idx % 25 == 0 or out_idx == n_output:
                logger.info("  定标写入 %3d/%d 波段 (%.1f%%)",
                            out_idx, n_output, 100.0 * out_idx / n_output)

        out_ds.FlushCache()
        out_ds = None
        band_plan = None
        if MASK_NODATA and masked_pixels:
            logger.info("  掩膜无效像元累计: %d 个", masked_pixels)

        # ── 2) RPC 正射：辐亮度上做几何重采样 ──
        ortho_info = orthorectify_radiance(tmp_native, footprint, tmp_ortho, logger)
        if ortho_info["band_count"] != n_output:
            raise RuntimeError("正射成果波段数 %d != %d"
                               % (ortho_info["band_count"], n_output))

        # ── 3) 写元数据（正射成果） ──
        update_ds = gdal.Open(tmp_ortho, gdal.GA_Update)
        if update_ds is None:
            raise RuntimeError("无法打开正射成果: %s" % tmp_ortho)
        try:
            update_ds.SetMetadata(None, "RPC")
        except Exception:
            pass
        update_ds.SetMetadataItem("SOURCE_TAR", os.path.basename(tar_path))
        update_ds.SetMetadataItem("SOURCE_TIFFS",
                                  "%s;%s" % (os.path.basename(vn_tif),
                                             os.path.basename(sw_tif)))
        update_ds.SetMetadataItem("BANDS_SOURCE", str(n_original))
        update_ds.SetMetadataItem("BANDS_OUTPUT", str(n_output))
        update_ds.SetMetadataItem(
            "REMOVED_BANDS_1BASED",
            ",".join(str(int(index) + 1) for index in removed))
        update_ds.SetMetadataItem(
            "WATER_VAPOR_BANDS_UM",
            "; ".join("%g-%g" % (lo, hi) for lo, hi in WATER_VAPOR_BANDS))
        update_ds.SetMetadataItem("RADIOMETRIC_FORMULA", "L = Gain * DN + Bias")
        update_ds.SetMetadataItem("MASK_NODATA", "TRUE" if MASK_NODATA else "FALSE")
        update_ds.SetMetadataItem(
            "GEOMETRY", "RPC orthorectified (corrected RPC) → EPSG:%d @ %g m"
                        % (TARGET_EPSG, TARGET_RES))
        update_ds.SetMetadataItem("TARGET_CRS", "EPSG:%d" % TARGET_EPSG)
        update_ds.SetMetadataItem("TARGET_RES", "%g" % TARGET_RES)
        update_ds.SetMetadataItem("RESAMPLING", RESAMPLING)
        update_ds.SetMetadataItem(
            "RPC_GRID_CORRECTION",
            "sx=%.6f tx=%.1f sy=%.6f ty=%.1f residual=%.4f px"
            % (rpc_correction["sx"], rpc_correction["tx"],
               rpc_correction["sy"], rpc_correction["ty"],
               rpc_correction["residual"]))
        update_ds.SetMetadataItem(
            "FOOTPRINT_WGS84_TL_TR_BR_BL",
            "; ".join("%.6f,%.6f" % footprint[key]
                      for key in ("TL", "TR", "BR", "BL")))
        update_ds.SetMetadataItem("STATS_GRID", "native (original) grid, full-res")
        for out_idx, wavelength in enumerate(wavelengths[keep], start=1):
            band = update_ds.GetRasterBand(out_idx)
            band.SetNoDataValue(NODATA_OUT)
            band.SetDescription("%.4f um" % wavelength)
        update_ds.FlushCache()
        update_ds = None

        # ── 4) 原子替换成果，写波长 / 几何参数 / 统计 ──
        np.save(tmp_wl, wavelengths[keep])
        if args.keep_native:
            os.replace(tmp_native, paths["native"])
        else:
            if os.path.exists(tmp_native):
                os.remove(tmp_native)
        os.replace(tmp_ortho, paths["ortho"])
        os.replace(tmp_wl, paths["wavelengths"])
        write_georef_json(paths["georef"], scene, tar_path, vn_tif, sw_tif,
                          (n_cols, n_rows), ortho_info, n_output, n_original,
                          footprint, rpc_correction, wavelengths, keep, removed,
                          remove_groups, logger)
        write_band_stats_csv(stats, paths["stats"])

        # ── 5) 质检 ──
        if make_qc:
            try:
                write_quality_figure(paths["ortho"], wavelengths[keep],
                                     WATER_VAPOR_BANDS, paths["qc"], scene,
                                     logger, band_stats=stats)
            except Exception:
                logger.warning("质检图生成失败: %s", paths["qc"], exc_info=True)

        elapsed = time.time() - started
        logger.info("完成: %d → %d 波段（剔除 %d），正射 %d×%d @ %g m，"
                    "用时 %.1f s，大小 %s",
                    n_original, n_output, len(removed),
                    ortho_info["width"], ortho_info["height"], TARGET_RES,
                    elapsed, human_size(os.path.getsize(paths["ortho"])))
        logger.info("  成果目录: %s", scene_dir)
        return {"scene": scene, "status": "成功", "bands": n_output,
                "elapsed": elapsed,
                "message": "%d 波段（剔除 %d），正射 %d×%d，%.1f s"
                           % (n_output, len(removed), ortho_info["width"],
                              ortho_info["height"], elapsed)}
    except Exception:
        for path in (tmp_native, tmp_ortho, tmp_wl):
            try:
                if os.path.exists(path):
                    os.remove(path)
            except OSError:
                pass
        if created_dir:
            try:
                if not os.listdir(scene_dir):
                    os.rmdir(scene_dir)
            except OSError:
                pass
        raise
    finally:
        ds_vn = ds_swir = None


# ╔══════════════════════════════════════════════════════════════╗
# ║  批处理主流程                                                 ║
# ╚══════════════════════════════════════════════════════════════╝

def main(argv=None):
    args = parse_args(argv)
    out_root = os.path.abspath(args.out)
    os.makedirs(out_root, exist_ok=True)
    log_path = args.log or os.path.join(out_root, "step1_batch.log")
    logger = setup_logging(log_path, args.verbose)
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)

    patterns = tuple(args.pattern) if args.pattern else TAR_PATTERNS
    src_dir = os.path.abspath(args.src)
    try:
        scenes = discover_scenes(src_dir, patterns)
    except Exception as exc:
        logger.error("扫描源目录失败: %s", exc)
        return 2
    if args.scene:
        import fnmatch
        scenes = [path for path in scenes
                  if any(fnmatch.fnmatch(scene_name(path), pattern)
                         for pattern in args.scene)]
    if args.limit > 0:
        scenes = scenes[: args.limit]
    if not scenes:
        logger.error("在 %s 下未找到匹配 %s 的压缩包",
                     src_dir, ", ".join(patterns))
        return 2

    logger.info("=" * 72)
    logger.info("ZY1E AHSI 步骤1：定标 + 水汽/坏波段剔除 + RPC 正射")
    logger.info("源目录    : %s", src_dir)
    logger.info("输出根目录: %s", out_root)
    logger.info("待处理    : %d 景", len(scenes))
    logger.info("坏波段    : %d 个（%d → %d 波段）",
                len(REMOVE_BANDS_1BASED), N_VNIR + N_SWIR,
                N_VNIR + N_SWIR - len(REMOVE_BANDS_1BASED))
    logger.info("正射网格  : EPSG:%d @ %g m（%s）",
                TARGET_EPSG, TARGET_RES, RESAMPLING)
    logger.info("保留原始格网辐亮度: %s", "是" if args.keep_native else "否")
    logger.info("日志文件  : %s", log_path)
    logger.info("=" * 72)

    if args.dry_run:
        for idx, tar_path in enumerate(scenes, start=1):
            logger.info("[%d/%d] %s", idx, len(scenes), os.path.basename(tar_path))
            logger.info("        -> %s", os.path.join(out_root, scene_name(tar_path)))
        logger.info("dry-run 结束，未执行任何处理")
        return 0

    results = []
    seen = {}
    started = time.time()
    for idx, tar_path in enumerate(scenes, start=1):
        scene = scene_name(tar_path)
        logger.info("-" * 72)
        logger.info("(%d/%d) %s", idx, len(scenes), os.path.basename(tar_path))
        if scene in seen:
            logger.error("场景名重复，跳过: %s（已处理 %s）", scene, seen[scene])
            results.append({"scene": scene, "status": "跳过（重名）", "bands": None,
                            "elapsed": 0.0, "message": "场景名重复"})
            continue
        seen[scene] = tar_path
        try:
            results.append(process_scene(tar_path, out_root, args, logger))
        except Exception as exc:
            logger.exception("处理失败: %s", os.path.basename(tar_path))
            results.append({"scene": scene, "status": "失败", "bands": None,
                            "elapsed": 0.0, "message": str(exc)})

    elapsed_total = time.time() - started
    success = [r for r in results if r["status"] == "成功"]
    skipped = [r for r in results if r["status"].startswith("跳过")]
    failed = [r for r in results if r["status"] == "失败"]

    logger.info("=" * 72)
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %.1f s",
                len(success), len(skipped), len(failed), elapsed_total)
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
