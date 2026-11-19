# -*- coding: utf-8 -*-
"""
01_calibration_badbands.py
ZY1E AHSI 步骤1：辐射定标 + 坏波段剔除 + RPC 正射 + 质检（v1.0.4）

流程：读包内 VNIR/SWIR → 解析光谱/定标/足迹/RPC → 逐波段定标 L=Gain*DN+Bias
→ 按剔除表跳过 24 个坏波段 → RPC 正射到 EPSG:32650/30 m → 统计与质检。
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

# =================================================================
# 配置区（按需修改）
# =================================================================
SRC_DIR = r"E:\ZhouYang\test\ZY1E_AHSI"
OUT_ROOT = r"E:\ZhouYang\test\ZY1E_output_step1"
TAR_PATTERNS = ("*.tar", "*.tar.gz", "*.tgz")

# 波段参数
N_VNIR = 76
N_SWIR = 90

# 坏波段剔除表（1-based）
REMOVE_BAND_GROUPS = [
    ("SWIR 首段低信噪", (77, 78, 79), 1006, 1039),
    ("1.4 μm 水汽吸收", (98, 99, 100, 101), 1358, 1408),
    ("1.4 μm 邻近低信噪", (102, 103, 104, 105), 1425, 1476),
    ("1.812 μm 水汽吸收", (125,), 1812, 1812),
    ("1.9 μm 水汽吸收", (126, 127, 128, 129, 130, 131), 1829, 1913),
    ("1.9 μm 邻近低信噪", (132, 133, 134), 1930, 1964),
    ("SWIR 末段低信噪", (164, 165, 166), 2468, 2501),
]
WATER_VAPOR_BANDS = [(1.35, 1.43), (1.80, 1.96)]

CAL_HAS_BIAS = True
# 掩膜设置：如果确认影像中 DN=0 的像元为无效黑边，请将 MASK_NODATA 改为 True
MASK_NODATA = False
NODATA_IN = 0
NODATA_OUT = -9999.0

# 正射参数
TARGET_EPSG = 32650
TARGET_RES = 30.0
RESAMPLING = "bilinear"
RPC_FIT_MAX_RESIDUAL = 2.0
REQUIRE_GEOREF = True
WRITE_GEOREF = True
KEEP_NATIVE_RADIANCE = False

# 批处理参数
OVERWRITE = False
GDAL_CACHE_MB = 1024
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]
GEOTIFF_COMPRESS = "DEFLATE"
GEOTIFF_COMPRESS_LEVEL = 6

# 质检参数
MAKE_QC = True
QC_MAX_SIZE = 1400
QC_DPI = 150

# 输出文件名：以缩写标注定标、坏波段剔除和正射处理
OUT_ORTHO_SUFFIX = "_cal_badband_ortho.tif"
LEGACY_ORTHO_NAME = "radiance_ortho.tif"
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


# =================================================================
# 通用工具
# =================================================================
def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(description="ZY1E AHSI 步骤1 修复版")
    parser.add_argument("--src", default=SRC_DIR)
    parser.add_argument("--out", default=OUT_ROOT)
    parser.add_argument("--pattern", action="append", default=None)
    parser.add_argument("--scene", action="append", default=None)
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--keep-native", action="store_true", default=KEEP_NATIVE_RADIANCE)
    parser.add_argument("--no-qc", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log", default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)

def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    
    # 增加 flush=True 确保终端实时显示
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(fmt)
    logger.addHandler(console)
    
    if log_path:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
            file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(fmt)
            logger.addHandler(file_handler)
        except OSError:
            pass
    logger.propagate = False
    return logger

def human_size(num_bytes):
    """把字节数格式化为 B/KB/MB/GB/TB 文本。"""
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024.0 or unit == "TB":
            return "%.1f %s" % (size, unit)
        size /= 1024.0

def discover_scenes(src_dir, patterns):
    """扫描源目录下匹配后缀的压缩包，返回排序后的文件路径列表。"""
    if not os.path.isdir(src_dir):
        raise NotADirectoryError("源目录不存在: %s" % src_dir)
    found = set()
    for pattern in patterns:
        found.update(glob.glob(os.path.join(src_dir, pattern), recursive=True))
    return sorted(f for f in found if os.path.isfile(f))

def scene_name(tar_path):
    """由压缩包文件名推导景名（去掉 .tar/.tar.gz/.tgz 后缀）。"""
    name = os.path.basename(tar_path)
    lower = name.lower()
    for ext in (".tar.gz", ".tgz", ".tar"):
        if lower.endswith(ext): return name[: -len(ext)]
    return os.path.splitext(name)[0]

def find_tar_member(names, required, suffixes, label, logger):
    """在包内成员名中按“关键词全部包含 + 后缀命中”定位文件，返回首个匹配。"""
    hits = []
    for name in names:
        base = os.path.basename(name).lower()
        if not base.endswith(suffixes): continue
        if all(key.lower() in base for key in required): hits.append(name)
    if not hits:
        raise FileNotFoundError(f"压缩包中未找到 {label}（关键词 {required}，后缀 {suffixes}）")
    hits.sort()
    if len(hits) > 1:
        logger.warning("%s 匹配到 %d 个文件，使用第一个: %s", label, len(hits), os.path.basename(hits[0]))
    return hits[0]

def decode_text(raw):
    """按 UTF-8 宽松解码并统一换行、去掉 BOM。"""
    text = raw.decode("utf-8", errors="ignore")
    return text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")

def read_member_text(tar, member, label):
    """读取 tar 内文本成员并解码；成员不可读时抛 FileNotFoundError。"""
    handle = tar.extractfile(member)
    if handle is None: raise FileNotFoundError(f"{label} 无法读取: {member}")
    with handle: return decode_text(handle.read())


# =================================================================
# 产品元数据解析
# =================================================================
def parse_footprint(xml_text, logger):
    """从产品 XML 读取四个角点经纬度，返回 {"TL","TR","BR","BL": (lon, lat)}。"""
    """从产品 XML 解析四角经纬度，返回 {TL/TR/BR/BL: (lon, lat)}。"""
    footprint = {}
    corner_map = {"TL": "TopLeft", "TR": "TopRight", "BR": "BottomRight", "BL": "BottomLeft"}
    for abbr, corner in corner_map.items():
        lat_match = re.search(r"<%sLatitude>\s*([-+0-9.eE]+)\s*</%sLatitude>" % (corner, corner), xml_text, re.IGNORECASE)
        lon_match = re.search(r"<%sLongitude>\s*([-+0-9.eE]+)\s*</%sLongitude>" % (corner, corner), xml_text, re.IGNORECASE)
        if not lat_match or not lon_match: raise ValueError(f"XML 缺少 {corner} 经纬度")
        footprint[abbr] = (float(lon_match.group(1)), float(lat_match.group(1)))
    return footprint

def parse_rpb(rpb_text):
    """从 .rpb 文件中提取 RPC 参数，兼容 XML 和 key=value 格式。"""
    rpc_dict = {}
    for match in re.finditer(r"([a-zA-Z_]+)\s*=\s*(.*?);", rpb_text, re.DOTALL):
        key = match.group(1).lower()
        val = match.group(2).strip()
        
        key_map = {
            "errbias": "ERR_BIAS", "errrand": "ERR_RAND",
            "lineoffset": "LINE_OFF", "sampleoffset": "SAMP_OFF", "sampoffset": "SAMP_OFF",
            "latoffset": "LAT_OFF", "longoffset": "LONG_OFF", "heightoffset": "HEIGHT_OFF",
            "linescale": "LINE_SCALE", "samplescale": "SAMP_SCALE", "sampscale": "SAMP_SCALE",
            "latscale": "LAT_SCALE", "longscale": "LONG_SCALE", "heightscale": "HEIGHT_SCALE",
            "linenumcoef": "LINE_NUM_COEFF", "linedencoef": "LINE_DEN_COEFF",
            "sampnumcoef": "SAMP_NUM_COEFF", "sampdencoef": "SAMP_DEN_COEFF"
        }
        if key in key_map:
            if "coef" in key:
                val = val.replace("(", "").replace(")", "").replace("\n", " ").replace("\r", " ").replace(",", " ")
                val = re.sub(r"\s+", " ", val).strip()
            rpc_dict[key_map[key]] = val
            
    if not rpc_dict:
        for key in ["ERR_BIAS", "ERR_RAND", "LINE_OFF", "SAMP_OFF", "LAT_OFF", "LONG_OFF", 
                    "HEIGHT_OFF", "LINE_SCALE", "SAMP_SCALE", "LAT_SCALE", "LONG_SCALE", "HEIGHT_SCALE",
                    "LINE_NUM_COEFF", "LINE_DEN_COEFF", "SAMP_NUM_COEFF", "SAMP_DEN_COEFF"]:
            match = re.search(r"<%s>\s*([^<]+)\s*</%s>" % (key, key), rpb_text, re.IGNORECASE)
            if match: rpc_dict[key] = match.group(1).strip()
    return rpc_dict

def parse_wavelengths(text, n_bands, label, logger):
    """解析光谱响应文件的中心波长（nm→μm），返回前 n_bands 个值。"""
    """解析光谱响应文件，取每行第一列（或第二列）中心波长并换算为 μm。"""
    values = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"): continue
        parts = re.split(r"[,\s]+", line)
        if not parts: continue
        try:
            if len(parts) >= 2 and float(parts[0]).is_integer() and float(parts[0]) < 200:
                values.append(float(parts[1]))
            else:
                values.append(float(parts[0]))
        except ValueError: continue
    if len(values) < n_bands: raise ValueError(f"{label} 数据不足 {n_bands} 行")
    wavelengths = np.asarray(values[:n_bands], dtype=np.float64) / 1000.0
    if not np.all(np.isfinite(wavelengths)) or np.any(wavelengths <= 0):
        raise ValueError(f"{label} 波长异常（存在非正数或 NaN）")
    return wavelengths

def parse_calibration(text, n_bands, has_bias, label, logger):
    """解析辐射定标文件，返回 (gain, bias) 两个数组，长度均为 n_bands。"""
    """解析辐射定标系数，返回 (gain 数组, bias 数组)，各取前 n_bands 行。

    has_bias=True 时每行需为 [波段号?, Gain, Bias]；否则 Bias 置 0。
    """
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"): continue
        tokens = re.split(r"[,\s]+", line)
        numbers = []
        for token in tokens:
            try: numbers.append(float(token))
            except ValueError: break
        if numbers:
            if len(numbers) >= 3 and numbers[0].is_integer() and numbers[0] < 200:
                rows.append(numbers[1:3])
            elif len(numbers) == 2 and numbers[0].is_integer() and numbers[0] < 200:
                rows.append([numbers[1]]) 
            else:
                rows.append(numbers[:2])
                
    if not rows: raise ValueError(f"{label} 定标系数为空")
    if has_bias is None:
        has_bias = all(len(row) >= 2 for row in rows)
    if has_bias:
        bad_rows = [i for i, row in enumerate(rows) if len(row) < 2]
        if bad_rows: raise ValueError(f"{label} 定标文件第 {bad_rows[0]+1} 行缺 Bias")
        data = np.asarray([row[:2] for row in rows], dtype=np.float64)
    else:
        data = np.asarray([[row[0], 0.0] for row in rows], dtype=np.float64)
        
    if len(data) < n_bands: raise ValueError(f"{label} 定标系数不足 {n_bands} 行")
    return data[:n_bands, 0], data[:n_bands, 1]


# =================================================================
# 坏波段剔除与几何处理
# =================================================================
def resolve_band_selection(wavelengths, logger):
    """按 REMOVE_BAND_GROUPS（1-based）生成保留波段索引、剔除列表与分组信息。"""
    """按剔除表（0-based 转换）返回 (保留波段索引, 被剔除索引, 分组说明)。"""
    n_total = int(wavelengths.size)
    if REMOVE_BANDS_1BASED:
        invalid = [b for b in REMOVE_BANDS_1BASED if not (1 <= b <= n_total)]
        if invalid: raise ValueError(f"波段号越界: {invalid[:10]}")
        removed = sorted({int(b) - 1 for b in REMOVE_BANDS_1BASED})
        groups = []
        for label, indexes, wl_lo, wl_hi in REMOVE_BAND_GROUPS:
            idx0 = sorted({int(b) - 1 for b in indexes})
            if not idx0: continue
            groups.append({"label": label, "bands_1based": list(indexes),
                           "wavelength_nm": [float(wavelengths[idx0[0]]) * 1000.0, float(wavelengths[idx0[-1]]) * 1000.0]})
    else:
        removed = [idx for idx, wl in enumerate(wavelengths) if any(lo <= wl <= hi for lo, hi in WATER_VAPOR_BANDS)]
        groups = []
    keep = np.asarray([i for i in range(n_total) if i not in set(removed)], dtype=np.int64)
    if keep.size <= 0: raise ValueError("剔除后没有剩余波段")
    return keep, removed, groups

def open_vsitar_dataset(tar_path, internal_path, label):
    """通过 /vsitar/ 直接打开压缩包内影像，避免解压到磁盘。"""
    """用 /vsitar/ 直接打开压缩包内影像，避免整包解压。"""
    vsi_path = "/vsitar/" + os.path.abspath(tar_path).replace("\\", "/") + "/" + internal_path
    dataset = gdal.Open(vsi_path)
    if dataset is None: raise RuntimeError(f"无法通过 /vsitar/ 打开 {label}")
    return dataset

def compute_rpc_grid_correction(ds_vn, footprint, logger):
    """用四角点 RPC 反算行列偏差，修正 RPC 的 OFF/SCALE 参数。

    做法：把四角经纬度经 RPC 反变换得到“真实”像素行列，再对像素角点
    做一次线性拟合，得到行/列方向的缩放 sx/sy 与平移 tx/ty；残差过大
    （> RPC_FIT_MAX_RESIDUAL）说明不是纯平移缩放，直接报错。
    """
    """用产品四角经纬度反算像素坐标，与影像四角做最小二乘拟合。

    返回 (修正后的 RPC 字典, 拟合参数)。修正只作用于 SAMP/LINE 的
    OFF/SCALE（dx = sx*col + tx，同理 dy），残差超过阈值则报错。
    """
    rpc = ds_vn.GetMetadata("RPC")
    if not rpc: raise ValueError("影像缺少 RPC 元数据")
    
    gt = ds_vn.GetGeoTransform()
    if not (abs(abs(gt[1]) - 1) < 1e-6 and abs(abs(gt[5]) - 1) < 1e-6):
        raise ValueError(f"交付影像 GeoTransform 非单位阵，RPC 校正失效: {gt}")

    transformer = gdal.Transformer(ds_vn, None, ["METHOD=RPC"])
    width, height = ds_vn.RasterXSize, ds_vn.RasterYSize
    pixel_corners = np.array([[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]])
    raw = []
    for corner in ("TL", "TR", "BR", "BL"):
        lon, lat = footprint[corner]
        ok, point = transformer.TransformPoint(1, lon, lat)
        if not ok: raise ValueError(f"RPC 反变换失败: {corner}")
        raw.append([point[0], point[1]])
    raw = np.asarray(raw, dtype=np.float64)

    (sx, tx) = np.linalg.lstsq(np.column_stack([pixel_corners[:, 0], np.ones(4)]), raw[:, 0], rcond=None)[0]
    (sy, ty) = np.linalg.lstsq(np.column_stack([pixel_corners[:, 1], np.ones(4)]), raw[:, 1], rcond=None)[0]
    predicted = np.column_stack([pixel_corners[:, 0] * sx + tx, pixel_corners[:, 1] * sy + ty])
    residual = float(np.max(np.hypot(predicted[:, 0] - raw[:, 0], predicted[:, 1] - raw[:, 1])))
    
    if residual > RPC_FIT_MAX_RESIDUAL: raise ValueError(f"RPC 拟合残差过大: {residual}")
    corrected = dict(rpc)
    corrected["SAMP_OFF"] = str((float(rpc["SAMP_OFF"]) - tx) / sx)
    corrected["SAMP_SCALE"] = str(float(rpc["SAMP_SCALE"]) / sx)
    corrected["LINE_OFF"] = str((float(rpc["LINE_OFF"]) - ty) / sy)
    corrected["LINE_SCALE"] = str(float(rpc["LINE_SCALE"]) / sy)
    return corrected, {"sx": sx, "tx": tx, "sy": sy, "ty": ty, "residual": residual}

def write_cutline_vsimem(footprint, logger):
    """把四角点足迹写成 /vsimem/ 下的 GeoJSON，作为 gdal.Warp 的裁剪线。"""
    """把四角足迹写成 /vsimem/ 下的 GeoJSON 裁切线，供 gdal.Warp 使用。"""
    ring = [footprint[key] for key in ("TL", "TR", "BR", "BL")] + [footprint["TL"]]
    geojson = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {},
               "geometry": {"type": "Polygon", "coordinates": [[[float(lon), float(lat)] for lon, lat in ring]]}}]}
    path = "/vsimem/zy1e_step1_footprint.geojson"
    gdal.FileFromMemBuffer(path, json.dumps(geojson).encode("utf-8"))
    return path

def build_gtiff_creation_options(band_interleave=True):
    """组装 GeoTIFF 创建选项（压缩、预测器、分块等）。"""
    """组装 GeoTIFF 创建选项：BIGTIFF、按波段存储、DEFLATE+PREDICTOR=3。"""
    creation = list(GTIFF_OPTIONS)
    if band_interleave: creation.append("INTERLEAVE=BAND")
    if GEOTIFF_COMPRESS: creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS, "PREDICTOR=3", "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL])
    return creation

def orthorectify_radiance(native_path, footprint, ortho_tmp, logger):
    """用 RPC + 足迹裁剪线把辐亮度正射到 TARGET_EPSG/TARGET_RES，返回网格信息。"""
    cutline = write_cutline_vsimem(footprint, logger)
    try:
        options = gdal.WarpOptions(
            format="GTiff", dstSRS="EPSG:%d" % TARGET_EPSG, xRes=TARGET_RES, yRes=TARGET_RES,
            resampleAlg=RESAMPLING, dstNodata=NODATA_OUT, targetAlignedPixels=True,
            multithread=True, creationOptions=build_gtiff_creation_options(),
            rpc=True, cutlineDSName=cutline, cropToCutline=True, outputType=gdal.GDT_Float32)
        result = gdal.Warp(ortho_tmp, native_path, options=options)
        if result is None: raise RuntimeError("RPC 正射失败")
        width, height = result.RasterXSize, result.RasterYSize
        geo_transform = result.GetGeoTransform()
        band_count = result.RasterCount
        result = None
    finally:
        try: gdal.Unlink(cutline)
        except: pass

    bounds = [geo_transform[0], geo_transform[3] + geo_transform[5] * height, geo_transform[0] + geo_transform[1] * width, geo_transform[3]]
    return {"width": int(width), "height": int(height), "band_count": int(band_count),
            "geotransform": [float(v) for v in geo_transform], "bounds": [float(v) for v in bounds]}


# =================================================================
# 统计输出与质检成图
# =================================================================
def write_georef_json(output_path, scene, tar_path, vn_tif, sw_tif, native_size, ortho_info, band_count, bands_source, footprint, rpc_correction, wavelengths, keep, removed, remove_groups, logger):
    """写 georef.json（网格、足迹、RPC 修正、剔除波段），供步骤2/3 追溯。"""
    payload = {
        "scene": scene, "source_tar": os.path.basename(tar_path), "source_tiffs": [os.path.basename(vn_tif), os.path.basename(sw_tif)],
        "geometry": "RPC orthorectified", "target_crs_epsg": TARGET_EPSG, "target_res": TARGET_RES, "resampling": RESAMPLING,
        "native_raster_size": [int(native_size[0]), int(native_size[1])], "raster_size": [int(ortho_info["width"]), int(ortho_info["height"])],
        "band_count": int(band_count), "bands_source": int(bands_source), "geotransform": ortho_info["geotransform"], "bounds": ortho_info["bounds"],
        "footprint_wgs84": footprint, "rpc_grid_correction": rpc_correction, "remove_band_groups": remove_groups,
        "removed_bands_1based": [int(i) + 1 for i in removed], "removed_wavelengths_um": [float(wavelengths[i]) for i in removed],
        "radiance_formula": "L = Gain * DN + Bias", "nodata": NODATA_OUT, "mask_nodata": bool(MASK_NODATA)
    }
    with open(output_path + ".tmp", "w", encoding="utf-8") as f: json.dump(payload, f, ensure_ascii=False, indent=2)
    os.replace(output_path + ".tmp", output_path)

def write_band_stats_csv(stats, output_path):
    """把逐波段统计写成 CSV（band_stats.csv）。"""
    tmp_path = output_path + ".tmp"
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "wavelength_um", "min", "max", "mean", "std", "p1", "p99", "valid_pixels"])
        for item in stats:
            writer.writerow([item["band"], "%.6f" % item["wavelength"], "%.6f" % item["min"], "%.6f" % item["max"],
                             "%.6f" % item["mean"], "%.6f" % item["std"], "%.6f" % item["p1"], "%.6f" % item["p99"], item["valid"]])
    os.replace(tmp_path, output_path)

def _nearest_band(wavelengths, target):
    """返回最接近 target（μm）的波段下标。"""
    return int(np.argmin(np.abs(wavelengths - target)))
def select_qc_bands(wavelengths):
    """为质检图挑选真彩色、假彩色与 SWIR 所用波段下标。"""
    return {"R": _nearest_band(wavelengths, 0.66), "G": _nearest_band(wavelengths, 0.55), "B": _nearest_band(wavelengths, 0.47),
            "NIR": _nearest_band(wavelengths, 0.86), "SWIR": _nearest_band(wavelengths, 2.20)}, None

def _stretch(array, valid):
    """2%–98% 线性拉伸到 0~1，仅用 valid 像元估计分位数。"""
    values = array[valid]
    if values.size == 0: return np.zeros_like(array)
    low, high = np.percentile(values, [2, 98])
    if high <= low: high = low + 1e-6
    return np.clip((array - low) / (high - low), 0.0, 1.0)

def load_qc_arrays(tif_path, band_map, unique):
    """降采样读取质检所需波段，返回 (数组字典, 绘图范围信息)。"""
    dataset = gdal.Open(tif_path)
    if dataset is None: return None
    width, height = dataset.RasterXSize, dataset.RasterYSize
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))
    buf_w, buf_h = max(1, width // scale), max(1, height // scale)
    
    arrays = {}
    for name, band_index in band_map.items():
        position = int(band_index) + 1
        arrays[name] = dataset.GetRasterBand(position).ReadAsArray(buf_xsize=buf_w, buf_ysize=buf_h).astype(np.float64)
        
    geo_transform = dataset.GetGeoTransform()
    extent = [geo_transform[0] / 1000.0, (geo_transform[0] + geo_transform[1] * width) / 1000.0,
              (geo_transform[3] + geo_transform[5] * height) / 1000.0, geo_transform[3] / 1000.0]
    info = {"extent": extent, "xlabel": "Easting (km)", "ylabel": "Northing (km)", "size_text": f"{width}x{height}"}
    return arrays, info

def write_quality_figure(tif_path, wavelengths, removed_ranges, output_png, scene, logger, band_stats=None):
    """生成 preview_quality.png：真彩色 / 假彩色 / SWIR / 均值光谱（灰带为剔除波段）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError: 
        logger.error("Matplotlib 未安装或导入失败，无法生成质检图。")
        return False
    
    band_map, _ = select_qc_bands(wavelengths)
    prepared = load_qc_arrays(tif_path, band_map, None)
    if prepared is None: 
        logger.error("无法读取正射成果用于生成质检图。")
        return False
    arrays, info = prepared
    
    def combine(names):
        channels = [arrays[name] for name in names]
        valid = np.ones_like(channels[0], dtype=bool)
        for ch in channels: valid &= np.isfinite(ch) & (ch != NODATA_OUT)
        rgb = np.ones(channels[0].shape + (3,), dtype=np.float64)
        for i, ch in enumerate(channels): rgb[..., i] = _stretch(ch, valid)
        rgb[~valid] = 1.0
        return rgb, valid

    true_rgb, valid = combine(("R", "G", "B"))
    false_rgb, _ = combine(("NIR", "R", "G"))
    swir = arrays["SWIR"].copy()
    swir_valid = np.isfinite(swir) & (swir != NODATA_OUT)
    swir_stretched = _stretch(swir, swir_valid)
    swir_stretched[~swir_valid] = np.nan

    figure, axes = plt.subplots(2, 2, figsize=(13, 10)) if band_stats else plt.subplots(1, 3, figsize=(15, 5.5))
    axes = axes.ravel()
    axes[0].imshow(true_rgb, extent=info["extent"], origin="upper"); axes[0].set_title("True color")
    axes[1].imshow(false_rgb, extent=info["extent"], origin="upper"); axes[1].set_title("False color")
    axes[2].imshow(swir_stretched, extent=info["extent"], origin="upper", cmap="gray"); axes[2].set_title("SWIR")
    
    if band_stats:
        wl_arr = np.array([i["wavelength"] for i in band_stats])
        mean_arr = np.array([i["mean"] for i in band_stats])
        axes[3].plot(wl_arr, mean_arr, color="black", linewidth=1.0)
        for low, high in removed_ranges: axes[3].axvspan(low, high, color="0.88", zorder=0)
        axes[3].set_title("Mean spectrum (shaded = removed)")
        
    figure.savefig(output_png, dpi=QC_DPI, bbox_inches="tight")
    plt.close(figure)
    return True


# =================================================================
# 单景核心处理逻辑
# =================================================================
def resolve_scene_paths(scene_dir):
    """按场景目录拼出本步各成果文件路径。"""
    scene = os.path.basename(os.path.normpath(scene_dir))
    return {k: os.path.join(scene_dir, v) for k, v in {
        "ortho": scene + OUT_ORTHO_SUFFIX, "native": OUT_NATIVE_NAME, "wavelengths": OUT_WAVELENGTHS_NAME,
        "georef": OUT_GEOREF_NAME, "stats": OUT_STATS_NAME, "qc": OUT_QC_NAME}.items()}

def process_scene(tar_path, out_root, args, logger):
    """处理单景：①元数据 ②打开影像/RPC ③波段规划 ④定标写盘

    ⑤RPC 正射 ⑥原子替换成果 ⑦元数据与统计 ⑧质检图。
    单景异常向上抛出，由 main 记录堆栈并计入失败，不中断批处理。
    """
    scene = scene_name(tar_path)
    scene_dir = os.path.join(out_root, scene)
    paths = resolve_scene_paths(scene_dir)
    make_qc = MAKE_QC and not args.no_qc

    required_files = ["ortho", "wavelengths", "georef", "stats"]
    if make_qc:
        required_files.append("qc")
    legacy_ortho = os.path.join(scene_dir, LEGACY_ORTHO_NAME)
    if (not args.overwrite and not os.path.isfile(paths["ortho"])
            and os.path.isfile(legacy_ortho)
            and all(os.path.isfile(paths[key])
                    for key in required_files if key != "ortho")):
        os.replace(legacy_ortho, paths["ortho"])
        logger.info("复用已有完整成果并更新正射文件名：%s", paths["ortho"])
    have_main = all(os.path.isfile(paths[p]) for p in required_files)

    if not args.overwrite and have_main:
        logger.info("跳过（成果已存在）: %s", scene_dir)
        return {"scene": scene, "status": "跳过", "bands": None, "elapsed": 0.0, "message": "已存在"}

    logger.info("开始处理场景: %s", scene)
    started = time.time()
    created_dir = not os.path.isdir(scene_dir)
    os.makedirs(scene_dir, exist_ok=True)
    tmp_native = os.path.join(scene_dir, "radiance_native.tmp.tif")
    tmp_ortho = os.path.join(scene_dir, scene + "_cal_badband_ortho.tmp.tif")
    tmp_wl = os.path.join(scene_dir, "wavelengths_clean.tmp.npy")
    ds_vn = ds_swir = None
    
    try:
        # ① 读取并解析元数据：光谱响应、定标系数、足迹与坏波段表
        logger.info("  - 正在读取并解析元数据...")
        with tarfile.open(tar_path, "r:*") as tar:
            names = tar.getnames()
            vn_tif = find_tar_member(names, ("_vn",), (".tif", ".tiff"), "VNIR 影像", logger)
            sw_tif = find_tar_member(names, ("_sw",), (".tif", ".tiff"), "SWIR 影像", logger)
            vn_rpb = find_tar_member(names, ("_vn",), (".rpb",), "VNIR RPC 文件", logger)
            vn_cal = find_tar_member(names, ("vnir", "radcal"), (".raw",), "VNIR 定标", logger)
            sw_cal = find_tar_member(names, ("swir", "radcal"), (".raw",), "SWIR 定标", logger)
            vn_sr = find_tar_member(names, ("vnir", "spectralresponse"), (".raw",), "VNIR 光谱", logger)
            sw_sr = find_tar_member(names, ("swir", "spectralresponse"), (".raw",), "SWIR 光谱", logger)
            xml_member = find_tar_member(names, (), (".xml",), "元数据 XML", logger)

            wavelengths = np.concatenate([
                parse_wavelengths(read_member_text(tar, vn_sr, "VNIR 光谱"), N_VNIR, "VNIR", logger),
                parse_wavelengths(read_member_text(tar, sw_sr, "SWIR 光谱"), N_SWIR, "SWIR", logger)])
            vn_gain, vn_bias = parse_calibration(read_member_text(tar, vn_cal, "VNIR 定标"), N_VNIR, CAL_HAS_BIAS, "VNIR", logger)
            sw_gain, sw_bias = parse_calibration(read_member_text(tar, sw_cal, "SWIR 定标"), N_SWIR, CAL_HAS_BIAS, "SWIR", logger)
            keep, removed, remove_groups = resolve_band_selection(wavelengths, logger)
            footprint = parse_footprint(read_member_text(tar, xml_member, "XML"), logger)

            vn_rpc_dict = {}
            try:
                vn_rpb_text = read_member_text(tar, vn_rpb, "VNIR RPC")
                vn_rpc_dict = parse_rpb(vn_rpb_text)
            except Exception as e:
                logger.warning("手动解析外部 .rpb 失败，将依赖 TIFF 内自带的 RPC: %s", e)

        # ② 打开 /vsitar/ 影像；外部 .rpb 解析失败时依赖 TIFF 自带 RPC
        logger.info("  - 正在打开影像并计算 RPC 校正参数...")
        ds_vn = open_vsitar_dataset(tar_path, vn_tif, "VNIR")
        ds_swir = open_vsitar_dataset(tar_path, sw_tif, "SWIR")
        
        if not ds_vn.GetMetadata("RPC") and vn_rpc_dict:
            ds_vn.SetMetadata(vn_rpc_dict, "RPC")
            logger.info("    已注入外部解析的 RPC 参数至数据集")

        rpc_corrected, rpc_correction = compute_rpc_grid_correction(ds_vn, footprint, logger)
        
        n_rows, n_cols = ds_vn.RasterYSize, ds_vn.RasterXSize
        n_original, n_output = N_VNIR + N_SWIR, int(keep.size)

        # ③ 规划输出波段：剔除波段直接跳过，按全局波段号取对应 VNIR/SWIR 定标系数
        band_plan = []
        for global_idx in keep:
            if global_idx < N_VNIR:
                band_plan.append((ds_vn, int(global_idx) + 1, vn_gain[global_idx], vn_bias[global_idx], wavelengths[global_idx]))
            else:
                local = int(global_idx) - N_VNIR
                band_plan.append((ds_swir, local + 1, sw_gain[local], sw_bias[local], wavelengths[global_idx]))

        driver = gdal.GetDriverByName("GTiff")
        if os.path.exists(tmp_native): os.remove(tmp_native)
        out_ds = driver.Create(tmp_native, n_cols, n_rows, n_output, gdal.GDT_Float32, options=build_gtiff_creation_options())
        out_ds.SetGeoTransform(ds_vn.GetGeoTransform())
        out_ds.SetProjection(ds_vn.GetProjection())
        out_ds.SetMetadata(rpc_corrected, "RPC")

        # ④ 逐波段辐射定标：L = Gain * DN + Bias（可选 DN=0 掩膜），并统计
        logger.info("  - 开始逐波段执行辐射定标与坏波段剔除 (共 %d 个波段)...", n_output)
        calibration_started = time.time()
        stats = []
        for out_idx, (src_ds, src_band, gain, bias, wavelength) in enumerate(band_plan, start=1):
            dn = src_ds.GetRasterBand(src_band).ReadAsArray()
            radiance = dn.astype(np.float64) * float(gain) + float(bias)
            if MASK_NODATA:
                radiance[dn == NODATA_IN] = NODATA_OUT
            
            out_band = out_ds.GetRasterBand(out_idx)
            out_band.WriteArray(radiance.astype(np.float32))
            out_band.SetNoDataValue(NODATA_OUT)
            
            valid_values = radiance[radiance != NODATA_OUT] if MASK_NODATA else radiance[np.isfinite(radiance)]
            if valid_values.size:
                p1, p99 = np.percentile(valid_values[::4], [1, 99])
                stats.append({"band": out_idx, "wavelength": float(wavelength), "min": float(valid_values.min()),
                              "max": float(valid_values.max()), "mean": float(valid_values.mean()),
                              "std": float(valid_values.std()), "p1": float(p1), "p99": float(p99), "valid": int(valid_values.size)})
            else:
                stats.append({"band": out_idx, "wavelength": float(wavelength), "min": 0.0, "max": 0.0, "mean": 0.0, "std": 0.0, "p1": 0.0, "p99": 0.0, "valid": 0})
            
            # 每 10 个波段打印一次进度
            if out_idx == 1 or out_idx % 10 == 0 or out_idx == n_output:
                logger.info("    定标与写盘进度: %d/%d (%.1f%%)", out_idx, n_output, 100.0 * out_idx / n_output)
                
        out_ds.FlushCache()
        out_ds = None
        
        # ⑤ RPC 正射：用足迹裁剪线裁掉黑边，输出 EPSG:32650 / 30 m
        logger.info("  - 定标完成（%.1f s），开始执行 RPC 正射校正 (gdal.Warp 处理中，请耐心等待)...",
                    time.time() - calibration_started)
        ortho_info = orthorectify_radiance(tmp_native, footprint, tmp_ortho, logger)
        logger.info("  - 正射完成，成果尺寸: %d x %d", ortho_info["width"], ortho_info["height"])
        
        # ⑥ 原子替换成果文件：波长 → 正射 tif（先写 .tmp 再 os.replace）
        np.save(tmp_wl, wavelengths[keep])
        if args.keep_native: os.replace(tmp_native, paths["native"])
        else: os.remove(tmp_native)
        os.replace(tmp_ortho, paths["ortho"])
        os.replace(tmp_wl, paths["wavelengths"])
        
        # ⑦ 写 georef.json（格网/足迹/剔除信息）与 band_stats.csv
        logger.info("  - 正在写入几何与统计元数据文件...")
        write_georef_json(paths["georef"], scene, tar_path, vn_tif, sw_tif, (n_cols, n_rows), ortho_info, n_output, n_original, footprint, rpc_correction, wavelengths, keep, removed, remove_groups, logger)
        write_band_stats_csv(stats, paths["stats"])

        # ⑧ 质检图：真彩色 / 假彩色 / SWIR / 均值光谱
        if make_qc:
            logger.info("  - 正在生成质检图...")
            removed_ranges = [(g["wavelength_nm"][0]/1000, g["wavelength_nm"][1]/1000) for g in remove_groups if g["bands_1based"]]
            qc_success = write_quality_figure(paths["ortho"], wavelengths[keep], removed_ranges, paths["qc"], scene, logger, band_stats=stats)
            if not qc_success:
                raise RuntimeError("质检图生成失败")

        logger.info("场景处理完成: %s（%d 波段，用时 %.1f s）",
                    scene, n_output, time.time() - started)
        return {"scene": scene, "status": "成功", "bands": n_output, "elapsed": time.time() - started, "message": "OK"}
        
    except BaseException as e:
        for p in (tmp_native, tmp_ortho, tmp_wl):
            try:
                if os.path.exists(p): os.remove(p)
            except OSError: pass
        if created_dir and os.path.exists(scene_dir) and not os.listdir(scene_dir):
            try: os.rmdir(scene_dir)
            except OSError: pass
        raise e
    finally:
        ds_vn = ds_swir = None


# =================================================================
# 批处理主流程
# =================================================================
def main(argv=None):
    """批处理入口：扫描压缩包 → （dry-run 仅打印计划）→ 逐景处理 → 汇总。"""
    args = parse_args(argv)
    out_root = os.path.abspath(args.out)
    os.makedirs(out_root, exist_ok=True)
    log_path = args.log or os.path.join(out_root, "step1_batch.log")
    logger = setup_logging(log_path, args.verbose)
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)

    logger.info("=" * 60)
    logger.info("ZY1E AHSI 步骤1：辐射定标 + 坏波段剔除 + RPC 正射")
    logger.info("源目录  : %s", os.path.abspath(args.src))
    logger.info("输出目录: %s", out_root)
    logger.info("日志文件: %s", log_path)
    logger.info("=" * 60)

    patterns = tuple(args.pattern) if args.pattern else TAR_PATTERNS
    try:
        scenes = discover_scenes(os.path.abspath(args.src), patterns)
    except Exception as e:
        logger.error("扫描失败: %s", e); return 2
        
    if args.scene:
        import fnmatch
        scenes = [p for p in scenes if any(fnmatch.fnmatch(scene_name(p), pat) for pat in args.scene)]
    if args.limit > 0: scenes = scenes[:args.limit]
    if not scenes: logger.error("未找到匹配压缩包"); return 2

    if args.dry_run:
        for idx, p in enumerate(scenes, 1):
            scene = scene_name(p)
            logger.info("[%d/%d] %s -> %s", idx, len(scenes),
                        os.path.basename(p),
                        os.path.join(out_root, scene, scene + OUT_ORTHO_SUFFIX))
        logger.info("dry-run 完成：共 %d 景，未处理文件", len(scenes))
        return 0

    results = []
    seen = {}
    started = time.time()
    for idx, tar_path in enumerate(scenes, start=1):
        scene = scene_name(tar_path)
        logger.info("="*50)
        logger.info("正在处理第 %d/%d 景数据", idx, len(scenes))
        
        if scene in seen:
            logger.error("场景重名跳过: %s", scene); continue
        seen[scene] = tar_path
        try:
            results.append(process_scene(tar_path, out_root, args, logger))
        except Exception:
            logger.exception("处理失败: %s", os.path.basename(tar_path))
            results.append({"scene": scene, "status": "失败", "bands": None, "elapsed": 0.0, "message": "见日志"})

    success = sum(1 for r in results if r["status"] == "成功")
    skipped = sum(1 for r in results if r["status"].startswith("跳过"))
    failed = sum(1 for r in results if r["status"] == "失败")
    logger.info("="*50)
    logger.info("批处理结束: 成功 %d, 跳过 %d, 失败 %d，用时 %.1f s", success, skipped, failed, time.time() - started)
    for result in results:
        logger.info("  [%s] %s（%.1f s）", result["status"], result["scene"], result["elapsed"])
    if failed:
        logger.error("有 %d 景处理失败，详见日志: %s", failed, log_path)
        return 1
    logger.info("全部完成，成果位于: %s", out_root)
    return 0

if __name__ == "__main__":
    sys.exit(main())