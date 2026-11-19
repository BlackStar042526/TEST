# -*- coding: utf-8 -*-
"""
04_aloh_mineral_extraction.py
===================================
ZY1E AHSI 步骤4：Al-OH 蚀变组合提取（USGS splib07a 光谱库 + SAM / SFF / MF）

本步读取步骤3 的“分析就绪”区域镶嵌反射率，把 USGS splib07a 光谱库
（ChapterM_Minerals 中的 Al-OH 蚀变矿物）重采样到 AHSI 的 152 个波段，
逐像元做：
    - 连续统去除 Al-OH 2.2 μm 吸收特征（深度 / 位置 / 宽度）
    - 连续统去除后的 SAM（CR-SAM，逐矿物类；消除反照率与大气残差影响）
    - SFF 连续统去除光谱特征拟合（丰度代理 + 拟合残差）
    - MF 匹配滤波 + MTMF 不可行性（可选）
再按阈值合成 Al-OH 矿物分类图与置信度图。

输入
----
步骤3 输出根目录 IN_ROOT_STEP3 下每区一个子文件夹（本步只认 mosaic/）：
    <IN_ROOT_STEP3>/<区名>/mosaic/mosaic_reflectance.tif   152 波段区域镶嵌反射率
    <IN_ROOT_STEP3>/<区名>/mosaic/mosaic_mask.tif          0=无数据 1=清晰 2=云填充 3=云影填充
    <IN_ROOT_STEP3>/<区名>/mosaic/wavelengths_clean.npy    152 个波段中心波长（μm）
    <IN_ROOT_STEP3>/<区名>/mosaic/mosaic_info.json         网格/参与景信息（可缺）
光谱库（USGS splib07a ASCII 版）：
    <LIB_ROOT>/ChapterM_Minerals/*.txt
      - 首行为记录头，如 " splib07a Record=108: Actinolite HS116.1B  ASDFRb AREF"
      - 其后每行一个反射率值（0~1 绝对反射率 AREF）
    <LIB_ROOT>/splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt
    <LIB_ROOT>/splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt
    <LIB_ROOT>/splib07a_Wavelengths_AVIRIS_1996_0.37-2.5_microns.txt
    AHSI 光谱响应（FWHM）从原始 tar 内的
    ZY1E_AHSI_VNIR/SWIR_Spectralresponse.raw 读取（两列：中心波长 nm、FWHM nm），
    找不到时回退 VNIR 8.76 nm / SWIR 16.26 nm。

每区输出（OUT_ROOT_STEP4/<区名>/）
----------------------------------
    aloh_class.tif                Al-OH 矿物分类（uint8；0=未分类 255=无效）
    aloh_confidence.tif           分类置信度 0~100
    aloh_depth.tif                2.2 μm 吸收深度（连续统去除）
    aloh_position_um.tif          吸收中心位置（μm）
    aloh_width_nm.tif             半深宽（nm）
    sam_angle_deg.tif             逐矿物类 SAM 角（度，多波段）
    sff_abundance.tif             逐矿物类 SFF 丰度（0~1 相对值，多波段）
    sff_rms.tif                   逐矿物类 SFF 拟合残差（多波段）
    mf_abundance.tif              逐矿物类 MF 匹配滤波丰度（多波段，可选）
    mtmf_infeasibility.tif        逐矿物类 MTMF 不可行性（多波段，可选）
    classification_legend.csv     类别编号 / 矿物中英文 / 端元数 / 颜色
    aloh_summary.csv              逐类像元数、面积、占比、平均深度/置信度
    aloh_info.json                参数、库信息、覆盖率、类别统计（可追溯）
    preview_quality.png           质检图（真彩色/分类/深度/吸收位置/置信度/面积）

_library/（全流程共用，只需建一次）
    ahsi_library.npy              [N, 152] 传感器匹配后的库光谱（float32）
    ahsi_library_meta.csv         端元清单（矿物类 / 样本 / 光谱仪 / 源文件）
    ahsi_band_spec.csv            152 波段：中心波长 / FWHM / 是否参与匹配 / 是否 Al-OH 窗口
    library_overview.png          库光谱图（全谱 + 2.1~2.5 μm 放大）
_catalog/
    step4_batch.log               批处理日志
    areas_summary.csv             本区汇总（类别面积/占比/参数）

处理要点
--------
1. 光谱库→传感器重采样：对每个 AHSI 波段中心做高斯卷积（σ=FWHM/2.3548），
   逐波段独立采样。**不能**把 152 个波长直接当单调网格做 np.interp：
   本数据 VNIR 末波段 1.0402 μm ≥ SWIR 首波段 1.0056 μm，波长非单调；
   另有 1.34-1.44、1.78-1.96 μm 水汽带已被步骤1 剔除，波长序列存在缺口。
2. 波长上限：库光谱只到 2.5 μm；本数据 2.35 μm 以后辐亮度只有 1~6 DN、
   6S 强水汽吸收会把噪声放大成“地表反射率 > 表观 TOA 反射率”的伪值，
   因此参与匹配（SAM/SFF/MF）的波段一律限制在 ≤2.30 μm（实测见
   `_catalog/swir_diagnostic.md`）；aloh_* 索引图仍覆盖 2.08~2.30 μm。
3. 掩膜：只用 mask==1 的清晰像元；反射率须在 [REFL_MIN, REFL_MAX] 内；
   连续统 ≤0 的像元判为无效。
4. 分类规则：先从每类端元中取最佳（CR-SAM 角最小 / SFF 残差最小），
   再要求 吸收深度 ≥ MIN_ALOH_DEPTH、CR-SAM 角 ≤ SAM_MAX_DEG、
   SFF 残差 ≤ SFF_MAX_RMS，全部满足才归入该矿物类，否则记为 0（未分类）。
   注意：本数据 2.35 μm 以后反射率快速趋近 0（长波 SWIR 信噪比很低），
   因此 Al-OH 连续统窗口取 2.08~2.30 μm、特征区间 2.12~2.28 μm；
   同时抽稀统计“坏波段”（缺测比例 >2% 的波段整体剔除），单像元允许 1~2 个
   匹配波段缺测并按 0 填充，避免因个别波段 nodata 丢掉一半像元。
5. MF/MTMF：背景均值/协方差从清晰像元抽稀采样估计（151×151），
   加小岭参数求伪逆；MTMF 不可行性用白化空间残差近似，
   只作为辅助判据/质检，不参与默认分类判决。
6. 流式处理：按 BLOCK_ROWS 行分块读取（152 波段 × 块），逐块写出，
   避免整景 1.8 GB 立方体全读入内存。
7. 稳健性：原子写出、跳过已完成、单区失败不影响整批、
   波段数/波长数/库波段一致性校验、缺库/缺响应文件明确报错。

用法
----
    python 04_aloh_mineral_extraction.py --dry-run
    python 04_aloh_mineral_extraction.py --area "2.LSK"            # 试点
    python 04_aloh_mineral_extraction.py --area "2.LSK" --overwrite
    python 04_aloh_mineral_extraction.py --rebuild-library         # 重建库
    python 04_aloh_mineral_extraction.py --max-angle 12 --min-depth 0.02
    python 04_aloh_mineral_extraction.py --no-mf

依赖：GDAL >= 3.0、numpy；matplotlib 可选（质检图）。
      建议用装好 GDAL 的 Anaconda 解释器，例如
      C:\\ProgramData\\anaconda3\\python.exe

后续约定：步骤5 读取本步 aloh_class.tif / sam_angle_deg.tif / sff_abundance.tif /
aloh_depth.tif + 步骤3 mosaic 与 source_scene 索引，做机器学习分类与精度评价；
_library/ahsi_library.npy 可直接复用。
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
import sys
import tarfile
import time

import numpy as np


# ╔══════════════════════════════════════════════════════════════╗
# ║  Windows / conda 的 GDAL DLL 引导                              ║
# ╚══════════════════════════════════════════════════════════════╝
# Python 3.8+ 的扩展模块不再用 PATH 查找依赖 DLL，需要显式
# os.add_dll_directory。conda 环境的 gdal.dll 在 Library\bin；
# 个别环境因 conda 事务中断，sqlite3.dll 只留在 pkgs 缓存里，
# 这里一并补进 DLL 搜索路径（不影响正常环境）。
def _bootstrap_gdal_dll():
    if os.name != "nt":
        return
    prefixes = []
    for key in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "VIRTUAL_ENV"):
        value = os.environ.get(key)
        if value and os.path.isdir(value):
            prefixes.append(value)
    prefixes.append(os.path.dirname(sys.executable))
    prefixes += [r"C:\ProgramData\anaconda3", r"D:\ProgramData\anaconda3"]
    seen = set()
    for prefix in prefixes:
        if not prefix or not os.path.isdir(prefix):
            continue
        candidates = [os.path.join(prefix, "Library", "bin")]
        candidates += glob.glob(os.path.join(prefix, "pkgs", "libsqlite-*", "Library", "bin"))
        for path in candidates:
            if path in seen or not os.path.isdir(path):
                continue
            seen.add(path)
            try:
                os.add_dll_directory(path)
            except OSError:
                pass


_bootstrap_gdal_dll()

try:
    from osgeo import gdal, osr
except ImportError as exc:
    raise SystemExit(
        "当前 Python 解释器缺少 GDAL（osgeo）：%s\n"
        "请用装好 GDAL 的 Anaconda 解释器运行本脚本，例如：\n"
        "    C:\\ProgramData\\anaconda3\\python.exe 04_aloh_mineral_extraction.py"
        % exc)

gdal.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")


# ╔══════════════════════════════════════════════════════════════╗
# ║  【配置区】按需修改                                            ║
# ╚══════════════════════════════════════════════════════════════╝

# ────────────【输入/输出】步骤3 成果目录、本步成果根目录 ────────────
IN_ROOT_STEP3 = r"E:\ZhouYang\test\ZY1E_output_step3"
OUT_ROOT_STEP4 = r"E:\ZhouYang\test\ZY1E_output_step4"

# ────────────【输入】USGS splib07a ASCII 光谱库目录 ────────────
LIB_ROOT = r"E:\ZhouYang\test\ASCIIdata_splib07a"
LIB_CHAPTER_DIR = "ChapterM_Minerals"

# ────────────【输入】原始 tar 目录（用于读 AHSI 光谱响应 FWHM） ────────────
TAR_DIR = r"E:\ZhouYang\test\ZY1E_data"
RESPONSE_SUFFIXES = ("_VNIR_Spectralresponse.raw", "_SWIR_Spectralresponse.raw")

# ────────────【参数】传感器响应（读不到 raw 时的回退值） ────────────
VNIR_FWHM_NM_FALLBACK = 8.76
SWIR_FWHM_NM_FALLBACK = 16.26
VNIR_SWIR_SPLIT_UM = 1.05      # 回退时按该波长分 VNIR/SWIR

# ────────────【参数】Al-OH 矿物清单（库记录首词必须精确等于英文名） ────────────
# 值：(中文名, 绘图颜色)
TARGET_MINERALS = [
    ("Kaolinite",       "高岭石",       "#5b7db1"),
    ("Dickite",         "地开石",       "#e3a857"),
    ("Halloysite",      "埃洛石",       "#8a5a44"),
    ("Illite",          "伊利石",       "#2a9d8f"),
    ("Muscovite",       "白云母/绢云母", "#e76f51"),
    ("Montmorillonite", "蒙脱石",       "#9b5de5"),
    ("Beidellite",      "贝得石",       "#f4a261"),
    ("Pyrophyllite",    "叶蜡石",       "#4d908e"),
    ("Alunite",         "明矾石",       "#f94144"),
    ("Diaspore",        "硬水铝石",     "#577590"),
    ("Topaz",           "黄玉",         "#43aa8b"),
]
MINERAL_EN = [item[0] for item in TARGET_MINERALS]
MINERAL_CN = {item[0]: item[1] for item in TARGET_MINERALS}
MINERAL_COLOR = {item[0]: item[2] for item in TARGET_MINERALS}

# ────────────【参数】光谱库筛选 ────────────
LIB_MAX_SAMPLES_PER_MINERAL = 6     # 每类最多保留几个端元（按光谱仪优先级/样本名排序）
LIB_MIN_ENDMEMBERS = 3              # 建库后至少要有这么多端元，否则报错
# 光谱仪优先顺序（数字越小越优先）：ASD 系列 2151 通道 → BECK 480 通道 → AVIRIS 224
LIB_SPECTROMETER_RANK = {
    "ASDFRa": 0, "ASDFRb": 1, "ASDHRa": 2, "ASDHRb": 3,
    "ASDNGa": 4, "ASDNGb": 5,
    "BECKa": 6, "BECKb": 7, "BECKc": 8, "BECKd": 9,
    "AVIRIS": 10,
}
# 光谱仪族 → (波长文件, 期望通道数)
LIB_GRIDS = {
    "ASD": ("splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt", 2151),
    "BECK": ("splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt", 480),
    "AVIRIS": ("splib07a_Wavelengths_AVIRIS_1996_0.37-2.5_microns.txt", 224),
}
LIB_WAVELENGTH_FILES = {
    "ASD": "splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt",
    "BECK": "splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt",
    "AVIRIS": "splib07a_Wavelengths_AVIRIS_1996_0.37-2.5_microns.txt",
}

# ────────────【参数】库→AHSI 重采样 ────────────
LIB_MAX_WAVELENGTH_UM = 2.50   # 库光谱波长上限；AHSI 末波段 2.5011 μm 不参与匹配
# 参与匹配（SAM/SFF/MF 的背景统计与丰度）的波长上限：本数据 2.35 μm 以后
# 辐亮度只有 1~6 DN、地表反射率大于表观 TOA 反射率（物理上不可能），
# 属于噪声被 6S 强水汽吸收放大的伪值，故匹配一律只用 ≤2.30 μm 的波段。
MATCH_MAX_WAVELENGTH_UM = 2.30

# ────────────【参数】Al-OH 2.2 μm 吸收特征 ────────────
# 连续统锚点（左右肩）：2.08 / 2.30 μm；特征评价范围：2.12~2.28 μm。
# 说明：本数据 2.35 μm 以后辐射/反射率快速趋近 0（长波 SWIR 信噪比很低），
# 若把右肩放到 2.45 μm 会把连续统压到近零、吸收深度失真，故窗口收窄到 2.30。
ALOH_CONTINUUM_RANGE = (2.08, 2.30)
ALOH_FEATURE_RANGE = (2.12, 2.28)
MIN_ALOH_BANDS = 8             # 窗口内波段数下限，不足则报错

# ────────────【参数】像元有效性 ────────────
REFL_MIN = -0.05               # 反射率下限（允许轻微负值，大气校正残留）
REFL_MAX = 1.20                # 反射率上限
NODATA_THRESHOLD = -9000.0     # 小于该值视为 nodata（-9999）
MIN_CONTINUUM_VALUE = 0.002    # 连续统绝对下限（反射率量级），低于此判为无效
USE_CLOUD_FILLED = False       # False: 只用 mask==1 清晰像元
NODATA_OUT = -9999.0
# 本数据 SWIR 个别波段（如 1.44 μm 附近）有约 40% 像元是步骤2 的 nodata，
# 这类“坏波段”整体剔除，而不是把一半像元丢掉。
MAX_BAND_MISSING_FRACTION = 0.02   # 清晰像元中缺测比例超过 2% 的波段不参与匹配
MAX_MISSING_BANDS_MATCH = 2        # 单像元参与匹配的波段中允许缺测的个数（缺测按 0 填充）

# ────────────【参数】分类阈值（先用默认值试点，按直方图/质检图再调） ────────────
MIN_ALOH_DEPTH = 0.005         # 2.2 μm 吸收深度下限（连续统去除后）
SAM_MAX_DEG = 10.0             # 连续统去除光谱角（CR-SAM）上限（度）
SFF_MAX_RMS = 0.30             # SFF 连续统去除拟合残差上限（本数据噪声大，先放宽）
CLASS_BY = "sam"               # 分类依据："sam"（默认）或 "sff"（按残差最小）

# ────────────【参数】MF / MTMF（辅助产物，不参与默认分类判决） ────────────
ENABLE_MF = True
BG_SAMPLE_PIXELS = 150000      # 背景统计抽样像元数
BG_RIDGE = 1.0e-6              # 协方差求逆的岭参数（相对迹）

# ────────────【参数】分块与压缩 ────────────
BLOCK_ROWS = 96
GDAL_CACHE_MB = 1024
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]
GEOTIFF_COMPRESS = "DEFLATE"
GEOTIFF_COMPRESS_LEVEL = 6

# ────────────【参数】质检图 ────────────
MAKE_QC = True
QC_MAX_SIZE = 1500
QC_DPI = 150
TRUE_COLOR_TARGETS = {"R": 0.66, "G": 0.55, "B": 0.47}

# 文件名约定
IN_MOSAIC_NAME = "mosaic_reflectance.tif"
IN_MASK_NAME = "mosaic_mask.tif"
IN_WAVELENGTHS_NAME = "wavelengths_clean.npy"
IN_INFO_NAME = "mosaic_info.json"

LIB_DIR_NAME = "_library"
CATALOG_DIR_NAME = "_catalog"
OUT_CLASS_NAME = "aloh_class.tif"
OUT_CONFIDENCE_NAME = "aloh_confidence.tif"
OUT_DEPTH_NAME = "aloh_depth.tif"
OUT_POSITION_NAME = "aloh_position_um.tif"
OUT_WIDTH_NAME = "aloh_width_nm.tif"
OUT_SAM_NAME = "sam_angle_deg.tif"
OUT_SFF_ABUNDANCE_NAME = "sff_abundance.tif"
OUT_SFF_RMS_NAME = "sff_rms.tif"
OUT_MF_NAME = "mf_abundance.tif"
OUT_MTMF_NAME = "mtmf_infeasibility.tif"
OUT_LEGEND_NAME = "classification_legend.csv"
OUT_SUMMARY_NAME = "aloh_summary.csv"
OUT_INFO_NAME = "aloh_info.json"
OUT_QC_NAME = "preview_quality.png"
CLASS_NODATA = 255

LOG_NAME = "zy1e.step4"


# ╔══════════════════════════════════════════════════════════════╗
# ║  通用工具                                                     ║
# ╚══════════════════════════════════════════════════════════════╝

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="ZY1E AHSI 步骤4：Al-OH 蚀变矿物提取（USGS splib07a + SAM/SFF/MF）")
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP3,
                        help="步骤3 成果根目录（默认: %(default)s）")
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP4,
                        help="本步成果根目录（默认: %(default)s）")
    parser.add_argument("--lib", dest="lib_root", default=LIB_ROOT,
                        help="USGS splib07a ASCII 光谱库目录（默认: %(default)s）")
    parser.add_argument("--area", action="append", default=None,
                        help="只处理匹配的研究区，可多次指定，如 \"2.LSK\" 或 \"*LSK*\"")
    parser.add_argument("--limit", type=int, default=0,
                        help="只处理前 N 个研究区（默认 0 = 全部）")
    parser.add_argument("--overwrite", action="store_true", default=False,
                        help="已有成果也重做（默认跳过）")
    parser.add_argument("--rebuild-library", action="store_true", default=False,
                        help="强制重建 AHSI 匹配库")
    parser.add_argument("--library-only", action="store_true", default=False,
                        help="只构建/检查 AHSI 匹配库，不做区域提取")
    parser.add_argument("--dry-run", action="store_true",
                        help="只打印处理计划，不执行")
    parser.add_argument("--no-mf", action="store_true", default=False,
                        help="关闭 MF/MTMF 计算（更快，产物更少）")
    parser.add_argument("--max-angle", type=float, default=SAM_MAX_DEG,
                        help="SAM 角上限（度，默认: %(default)s）")
    parser.add_argument("--min-depth", type=float, default=MIN_ALOH_DEPTH,
                        help="2.2 μm 吸收深度下限（默认: %(default)s）")
    parser.add_argument("--max-rms", type=float, default=SFF_MAX_RMS,
                        help="SFF 拟合残差上限（默认: %(default)s）")
    parser.add_argument("--class-by", choices=("sam", "sff"), default=CLASS_BY,
                        help="分类依据（默认: %(default)s）")
    parser.add_argument("--log", default=None,
                        help="日志文件路径（默认 <out>/_catalog/step4_batch.log）")
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
            handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
            handler.setLevel(logging.DEBUG)
            handler.setFormatter(fmt)
            logger.addHandler(handler)
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


def read_text_lines(path):
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:
        text = handle.read().lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")
    return text.split("\n")


# ╔══════════════════════════════════════════════════════════════╗
# ║  传感器波段规格（波长 + FWHM）                                  ║
# ╚══════════════════════════════════════════════════════════════╝

def load_sensor_fwhm(tar_dir, band_wavelengths, logger):
    """
    用 tar 内 ZY1E_AHSI_VNIR/SWIR_Spectralresponse.raw 给每个输出波段匹配 FWHM。
    返回 (fwhm_um, 来源说明)；找不到时用回退常量。
    """
    centers_nm = np.array([], np.float64)
    fwhms_nm = np.array([], np.float64)
    source = "回退常量"
    try:
        tar_list = sorted(glob.glob(os.path.join(tar_dir, "*.tar")) +
                          glob.glob(os.path.join(tar_dir, "*.tar.gz")) +
                          glob.glob(os.path.join(tar_dir, "*.tgz")))
        if tar_list:
            with tarfile.open(tar_list[0], "r:*") as tar:
                names = tar.getnames()
                for suffix in RESPONSE_SUFFIXES:
                    hits = [name for name in names if name.lower().endswith(suffix.lower())]
                    if not hits:
                        continue
                    handle = tar.extractfile(hits[0])
                    if handle is None:
                        continue
                    with handle:
                        raw = handle.read()
                    text = raw.decode("utf-8", errors="ignore")
                    for line in text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n").split("\n"):
                        line = line.strip()
                        if not line:
                            continue
                        parts = re.split(r"[,\s]+", line)
                        if len(parts) < 2:
                            continue
                        try:
                            centers_nm = np.append(centers_nm, float(parts[0]))
                            fwhms_nm = np.append(fwhms_nm, float(parts[1]))
                        except ValueError:
                            continue
            if centers_nm.size:
                source = "%s（响应曲线 %d 个中心）" % (
                    os.path.basename(tar_list[0]), centers_nm.size)
                logger.info("传感器响应: %s, FWHM %.2f~%.2f nm",
                            source, float(fwhms_nm.min()), float(fwhms_nm.max()))
    except Exception as exc:
        logger.warning("读取光谱响应失败，使用回退 FWHM: %s", exc, exc_info=True)

    band_wl_nm = np.asarray(band_wavelengths, np.float64) * 1000.0
    if centers_nm.size:
        idx = np.abs(band_wl_nm[:, None] - centers_nm[None, :]).argmin(axis=1)
        dist = np.abs(band_wl_nm - centers_nm[idx])
        fwhm_nm = fwhms_nm[idx].copy()
        bad = dist > 0.5  # 0.5 nm 以上认为没匹配上（理论上应精确命中）
        if np.any(bad):
            fwhm_nm[bad] = SWIR_FWHM_NM_FALLBACK
            logger.warning("%d 个波段未在响应曲线中找到对应中心，按 SWIR FWHM 回退",
                           int(bad.sum()))
    else:
        fwhm_nm = np.where(band_wl_nm < VNIR_SWIR_SPLIT_UM * 1000.0,
                           VNIR_FWHM_NM_FALLBACK, SWIR_FWHM_NM_FALLBACK)
        logger.info("传感器响应: %s（VNIR %.2f nm / SWIR %.2f nm）",
                    source, VNIR_FWHM_NM_FALLBACK, SWIR_FWHM_NM_FALLBACK)
    return fwhm_nm / 1000.0, source


# ╔══════════════════════════════════════════════════════════════╗
# ║  光谱库：扫描 / 解析 / 传感器重采样                             ║
# ╚══════════════════════════════════════════════════════════════╝

def parse_record_header(line):
    """解析 splib07a 记录头，返回 (样本名, 光谱仪, 类型) 或 None。"""
    line = line.strip().lstrip("\ufeff")
    match = re.match(r"^splib07a\s+Record=\d+:\s*(.*)$", line, re.IGNORECASE)
    if not match:
        return None
    rest = match.group(1).strip()
    parts = rest.rsplit(None, 2)
    if len(parts) != 3:
        return None
    sample, spectrometer, rtype = parts
    return sample.strip(), spectrometer.strip(), rtype.strip().upper()


def spectrometer_family(spectrometer):
    upper = spectrometer.upper()
    if upper.startswith("ASDFR") or upper.startswith("ASDHR") or upper.startswith("ASDNG"):
        return "ASD"
    if upper.startswith("BECK"):
        return "BECK"
    if upper.startswith("AVIRIS"):
        return "AVIRIS"
    return None


def load_lib_wavelength_grid(lib_root, family, logger, cache):
    if family in cache:
        return cache[family]
    filename = LIB_WAVELENGTH_FILES[family]
    path = os.path.join(lib_root, filename)
    if not os.path.isfile(path):
        raise FileNotFoundError("缺少光谱库波长文件: %s" % path)
    values = []
    for line in read_text_lines(path)[1:]:
        line = line.strip()
        if not line:
            continue
        try:
            values.append(float(line.split()[0]))
        except ValueError:
            continue
    grid = np.asarray(values, np.float64)
    expected = LIB_GRIDS[family][1]
    if grid.size != expected:
        raise ValueError("光谱库波长文件 %s 有 %d 个通道，预期 %d"
                         % (os.path.basename(path), grid.size, expected))
    logger.info("库波长网格 %s: %d 通道, %.4f ~ %.4f μm",
                family, grid.size, grid.min(), grid.max())
    cache[family] = grid
    return grid


def scan_library(lib_root, logger):
    """扫描矿物章，返回 Al-OH 目标端元候选记录列表。"""
    chapter = os.path.join(lib_root, LIB_CHAPTER_DIR)
    if not os.path.isdir(chapter):
        raise NotADirectoryError("光谱库矿物章目录不存在: %s" % chapter)
    paths = sorted(glob.glob(os.path.join(chapter, "*.txt")))
    logger.info("光谱库矿物章文件数: %d", len(paths))

    skipped_type = skipped_family = parse_fail = 0
    candidates = []
    for path in paths:
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                head = handle.readline()
        except OSError:
            parse_fail += 1
            continue
        parsed = parse_record_header(head)
        if parsed is None:
            parse_fail += 1
            continue
        sample, spectrometer, rtype = parsed
        if rtype != "AREF":
            skipped_type += 1
            continue
        token = sample.split()[0]
        if token not in MINERAL_CN:
            continue
        family = spectrometer_family(spectrometer)
        if family is None:
            skipped_family += 1
            continue
        rank = LIB_SPECTROMETER_RANK.get(spectrometer, 99)
        candidates.append({
            "path": path,
            "file": os.path.basename(path),
            "sample": sample,
            "mineral": token,
            "spectrometer": spectrometer,
            "family": family,
            "rank": rank,
        })
    logger.info("Al-OH 候选端元(AREF): %d；跳过 非AREF %d / 不支持光谱仪 %d / 头部解析失败 %d",
                len(candidates), skipped_type, skipped_family, parse_fail)

    by_mineral = {}
    for record in candidates:
        by_mineral.setdefault(record["mineral"], []).append(record)
    selected = []
    for mineral in MINERAL_EN:
        group = by_mineral.get(mineral, [])
        group.sort(key=lambda item: (item["rank"], item["sample"], item["file"]))
        chosen = group[:LIB_MAX_SAMPLES_PER_MINERAL]
        if not chosen:
            logger.warning("库中缺少矿物 %s 的 AREF 端元", mineral)
        else:
            selected.extend(chosen)
        logger.info("  %-16s 候选 %2d → 选用 %d（%s）", mineral, len(group),
                    len(chosen), ",".join(sorted({item["spectrometer"] for item in chosen})) or "-")
    if len(selected) < LIB_MIN_ENDMEMBERS:
        raise ValueError("Al-OH 端元数 %d < %d，请检查光谱库路径/清单"
                         % (len(selected), LIB_MIN_ENDMEMBERS))
    return selected


def load_spectrum_values(path):
    """
    读 splib07a ASCII 数值列（跳过首行记录头）。

    注意：splib07a ASCII 用 -1.23e34 作为“缺测”哨兵值（多见于 BECK 光谱仪的
    水汽带等位置），这里统一转成 NaN，交给后续卷积与插值处理。
    """
    values = []
    for line in read_text_lines(path)[1:]:
        line = line.strip()
        if not line:
            continue
        try:
            values.append(float(line.split()[0]))
        except (ValueError, IndexError):
            continue
    array = np.asarray(values, np.float64)
    if array.size:
        invalid = ~np.isfinite(array) | (np.abs(array) > 1.0e3)
        array[invalid] = np.nan
    return array


def _fill_spectrum_nan(spectrum):
    """按波长顺序对 NaN 做线性插值填补（只用于库光谱的少量缺测点）。"""
    spectrum = np.array(spectrum, np.float64, copy=True)
    valid = np.isfinite(spectrum)
    if valid.sum() < 5 or valid.all():
        return spectrum
    index = np.arange(spectrum.size)
    spectrum[~valid] = np.interp(index[~valid], index[valid], spectrum[valid])
    return spectrum


def gaussian_resample(wl_lib_um, spectra, band_wl_um, band_fwhm_um):
    """
    把库光谱（N, nlib）按 AHSI 波段中心做高斯卷积，返回 (N, nband) float64。

    逐波段独立采样：本数据 AHSI 波长非单调（1.0402 → 1.0056），不能整体插值。
    每个波段只取 ±3σ 内的库采样点并归一化权重；库光谱中的 NaN（缺测哨兵）
    不参与加权；窗口内没有有效采样点时返回 NaN。
    """
    spectra = np.asarray(spectra, np.float64)
    n_band = band_wl_um.size
    out = np.full((spectra.shape[0], n_band), np.nan, np.float64)
    coverage = np.zeros(n_band, np.float64)
    sigma = band_fwhm_um / 2.354820045
    spacing = float(np.median(np.abs(np.diff(wl_lib_um))))
    for index in range(n_band):
        center = band_wl_um[index]
        width = sigma[index]
        if not np.isfinite(center) or not np.isfinite(width) or width <= 0:
            continue
        mask = np.abs(wl_lib_um - center) <= 3.0 * width
        count = int(mask.sum())
        expected = max(1.0, (6.0 * width) / spacing)
        coverage[index] = min(1.0, count / expected)
        if count < 1:
            continue
        delta = (wl_lib_um[mask] - center) / width
        weights = np.exp(-0.5 * delta * delta)
        subset = spectra[:, mask]
        finite = np.isfinite(subset)
        weight_sum = finite @ weights
        weighted = np.where(finite, subset, 0.0) @ weights
        good = weight_sum > 0
        out[good, index] = weighted[good] / weight_sum[good]
    return out, coverage


def build_library(lib_root, band_wl, band_fwhm, out_library_dir, logger, force=False):
    """建/读 AHSI 匹配库，返回 dict。"""
    npy_path = os.path.join(out_library_dir, "ahsi_library.npy")
    meta_path = os.path.join(out_library_dir, "ahsi_library_meta.csv")
    band_path = os.path.join(out_library_dir, "ahsi_band_spec.csv")
    if not force and os.path.isfile(npy_path) and os.path.isfile(meta_path):
        library = np.load(npy_path)
        meta = []
        with open(meta_path, "r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):
                meta.append(row)
        if library.ndim != 2 or library.shape[1] != band_wl.size or len(meta) != library.shape[0]:
            logger.warning("已有匹配库与当前波段不一致，强制重建")
        else:
            logger.info("复用已有匹配库: %s (%d 端元 × %d 波段)",
                        npy_path, library.shape[0], library.shape[1])
            return _assemble_library(library, meta, band_wl, band_fwhm)

    logger.info("构建 AHSI 匹配库 ...")
    os.makedirs(out_library_dir, exist_ok=True)
    records = scan_library(lib_root, logger)

    cache = {}
    endmembers, meta_rows, coverage_rows = [], [], []
    skipped_mismatch = 0
    filled_values = 0
    for record in records:
        grid = load_lib_wavelength_grid(lib_root, record["family"], logger, cache)
        values = load_spectrum_values(record["path"])
        if values.size != grid.size:
            skipped_mismatch += 1
            logger.warning("通道数不一致，跳过: %s (%d vs %d)",
                           record["file"], values.size, grid.size)
            continue
        resampled, coverage = gaussian_resample(grid, values[None, :], band_wl, band_fwhm)
        row = resampled[0]
        if not np.all(np.isfinite(row)):
            filled_values += int(np.count_nonzero(~np.isfinite(row)))
            row = _fill_spectrum_nan(row)
        if not np.all(np.isfinite(row)):
            skipped_mismatch += 1
            logger.warning("缺测值过多，跳过: %s", record["file"])
            continue
        endmembers.append(row)
        coverage_rows.append(coverage)
        meta_rows.append({
            "endmember_id": len(meta_rows),
            "class_id": MINERAL_EN.index(record["mineral"]) + 1,
            "mineral_en": record["mineral"],
            "mineral_cn": MINERAL_CN[record["mineral"]],
            "sample": record["sample"],
            "spectrometer": record["spectrometer"],
            "source_file": record["file"],
        })
    if skipped_mismatch:
        logger.warning("共 %d 个端元因通道数不一致被跳过", skipped_mismatch)
    if filled_values:
        logger.info("库光谱缺测值填补: %d 个（splib07a -1.23e34 哨兵/检出限空值）",
                    filled_values)
    if len(endmembers) < LIB_MIN_ENDMEMBERS:
        raise ValueError("有效端元数 %d < %d，无法建库" % (len(endmembers), LIB_MIN_ENDMEMBERS))

    library = np.vstack(endmembers).astype(np.float32)
    np.save(npy_path, library)
    write_meta_csv(meta_path, meta_rows)
    coverage_mean = np.mean(np.vstack(coverage_rows), axis=0)

    with open(band_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["band", "wavelength_um", "fwhm_nm", "used_in_matching",
                         "in_aloh_window", "lib_coverage"])
        for index, (wl, fwhm) in enumerate(zip(band_wl, band_fwhm), start=1):
            writer.writerow([index, "%.6f" % wl, "%.3f" % (fwhm * 1000.0),
                             "TRUE" if wl <= MATCH_MAX_WAVELENGTH_UM else "FALSE",
                             "TRUE" if ALOH_CONTINUUM_RANGE[0] <= wl <= ALOH_CONTINUUM_RANGE[1] else "FALSE",
                             "%.3f" % coverage_mean[index - 1]])

    logger.info("匹配库: %d 端元 × %d 波段，%s", library.shape[0], library.shape[1],
                human_size(os.path.getsize(npy_path)))
    make_library_overview(out_library_dir, band_wl, library, meta_rows, logger)
    return _assemble_library(library, meta_rows, band_wl, band_fwhm)


def _assemble_library(library, meta_rows, band_wl, band_fwhm):
    class_ids = np.asarray([int(row["class_id"]) for row in meta_rows], np.int32)
    groups = []
    for class_index in range(1, len(MINERAL_EN) + 1):
        groups.append((class_index, np.where(class_ids == class_index)[0]))
    return {
        "spectra": np.asarray(library, np.float32),
        "meta": meta_rows,
        "class_ids": class_ids,
        "groups": groups,
        "band_wl": np.asarray(band_wl, np.float64),
        "band_fwhm": np.asarray(band_fwhm, np.float64),
        "match_idx": np.where(np.asarray(band_wl) <= MATCH_MAX_WAVELENGTH_UM)[0],
    }


def write_meta_csv(path, meta_rows):
    fields = ["endmember_id", "class_id", "mineral_en", "mineral_cn",
              "sample", "spectrometer", "source_file"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in meta_rows:
            writer.writerow(row)


def make_library_overview(out_dir, band_wl, library, meta_rows, logger):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        logger.warning("matplotlib 不可用，跳过库光谱图: %s", exc)
        return
    class_ids = np.asarray([int(row["class_id"]) for row in meta_rows])
    class_names = [item[0] for item in TARGET_MINERALS]
    figure, axes = plt.subplots(1, 2, figsize=(15, 5.6))
    for class_index, class_name in enumerate(class_names, start=1):
        idx = np.where(class_ids == class_index)[0]
        if idx.size == 0:
            continue
        color = MINERAL_COLOR[class_name]
        for row in idx:
            axes[0].plot(band_wl, library[row], color=color, linewidth=0.7, alpha=0.6,
                         label=class_name if row == idx[0] else None)
            zoom = (band_wl >= 2.0) & (band_wl <= 2.5)
            axes[1].plot(band_wl[zoom], library[row, zoom], color=color, linewidth=0.8,
                         alpha=0.7)
    axes[0].set_xlabel("Wavelength (um)")
    axes[0].set_ylabel("Reflectance")
    axes[0].set_title("AHSI-matched splib07a Al-OH endmembers (n=%d)" % library.shape[0])
    axes[0].grid(alpha=0.3)
    axes[1].axvspan(2.10, 2.40, color="0.9", zorder=0)
    axes[1].set_xlim(2.0, 2.5)
    axes[1].set_xlabel("Wavelength (um)")
    axes[1].set_title("Al-OH 2.2 um feature region")
    axes[1].grid(alpha=0.3)
    handles, labels = axes[0].get_legend_handles_labels()
    axes[0].legend(handles, labels, fontsize=7, ncol=2, loc="best")
    figure.tight_layout()
    path = os.path.join(out_dir, "library_overview.png")
    figure.savefig(path, dpi=130, bbox_inches="tight")
    plt.close(figure)
    logger.info("库光谱图: %s", path)


# ╔══════════════════════════════════════════════════════════════╗
# ║  区域处理                                                     ║
# ╚══════════════════════════════════════════════════════════════╝

def discover_areas(in_root, patterns=None):
    areas = []
    if not os.path.isdir(in_root):
        raise NotADirectoryError("步骤3 成果目录不存在: %s" % in_root)
    for name in sorted(os.listdir(in_root)):
        mosaic_dir = os.path.join(in_root, name, "mosaic")
        if not os.path.isdir(mosaic_dir):
            continue
        if not os.path.isfile(os.path.join(mosaic_dir, IN_MOSAIC_NAME)):
            continue
        if patterns and not any(fnmatch.fnmatch(name, pattern) for pattern in patterns):
            continue
        areas.append(name)
    return areas


def _nearest_band(wavelengths, target):
    return int(np.argmin(np.abs(np.asarray(wavelengths) - target)))


def _stretch(array, valid):
    values = array[valid]
    if values.size == 0:
        return np.zeros_like(array)
    low, high = np.percentile(values, [2, 98])
    if high <= low:
        high = low + 1e-6
    return np.clip((array - low) / (high - low), 0.0, 1.0)


def _valid_pixels(array):
    """有效值：有限且不是 -9999 nodata。"""
    return np.isfinite(array) & (array > NODATA_THRESHOLD)


def estimate_band_missing(ds, mask_ds, candidate_idx, logger, max_rows=200):
    """抽稀行统计“清晰像元中各候选波段的有效比例”，用于剔除坏波段。"""
    width, height = ds.RasterXSize, ds.RasterYSize
    step = max(1, height // max_rows)
    valid_counts = np.zeros(candidate_idx.size, np.float64)
    clear_total = 0
    for y in range(step // 2, height, step):
        cube = ds.ReadAsArray(0, y, width, 1)
        if cube is None:
            continue
        mask = mask_ds.ReadAsArray(0, y, width, 1)
        if mask.ndim == 2:
            mask = mask[0]
        clear = mask == 1
        if not np.any(clear):
            continue
        block = cube[candidate_idx][:, 0, :]
        valid = _valid_pixels(block)
        valid_counts += valid[:, clear].sum(axis=1)
        clear_total += int(np.count_nonzero(clear))
    if clear_total == 0:
        raise RuntimeError("采样行中没有清晰像元，检查掩膜")
    valid_fraction = valid_counts / float(clear_total)
    bad = valid_fraction < (1.0 - MAX_BAND_MISSING_FRACTION)
    logger.info("坏波段检查: 采样 %d 个清晰像元；有效比例 <%.0f%% 的波段 %d 个",
                clear_total, 100.0 * (1.0 - MAX_BAND_MISSING_FRACTION),
                int(np.count_nonzero(bad)))
    return valid_fraction, clear_total


def sample_background(ds, mask_ds, match_idx, target, logger):
    """从清晰像元抽稀采样，估计背景均值与协方差逆（n_match×n_match）。"""
    width, height = ds.RasterXSize, ds.RasterYSize
    step = max(1, height // 200)
    chunks = []
    total = 0
    for y in range(step // 2, height, step):
        cube = ds.ReadAsArray(0, y, width, 1)
        if cube is None:
            continue
        mask = mask_ds.ReadAsArray(0, y, width, 1)
        if mask.ndim == 2:
            mask = mask[0]
        block = cube[match_idx].reshape(len(match_idx), -1).astype(np.float64)
        finite = _valid_pixels(block)
        missing = np.count_nonzero(~finite, axis=0)
        low = np.where(finite, block, np.inf).min(axis=0)
        high = np.where(finite, block, -np.inf).max(axis=0)
        ok = ((mask.reshape(-1) == 1) & (missing <= MAX_MISSING_BANDS_MATCH) &
              (low > REFL_MIN) & (high < REFL_MAX))
        if not np.any(ok):
            continue
        chunks.append(np.where(finite[:, ok], block[:, ok], 0.0).T)
        total += int(ok.sum())
        if total >= target:
            break
    if not chunks:
        raise RuntimeError("无法获取清晰像元样本，检查掩膜/反射率范围")
    samples = np.vstack(chunks)
    if samples.shape[0] > target:
        stride = max(1, samples.shape[0] // target)
        samples = samples[::stride][:target]
    mean = samples.mean(axis=0)
    covariance = np.cov(samples, rowvar=False)
    ridge = BG_RIDGE * np.trace(covariance) / covariance.shape[0]
    covariance = covariance + np.eye(covariance.shape[0]) * ridge
    inverse = np.linalg.pinv(covariance, rcond=1e-10)
    logger.info("背景统计: %d 像元, 条件数 1e%.1f",
                samples.shape[0], np.log10(np.linalg.cond(covariance) + 1e-30))
    return mean, inverse


def _class_reduce(values, groups, reducer, fill):
    """values: (n_end, n) → (n_class, n)，按端元分组做 min/max。"""
    n_class = len(groups)
    out = np.full((n_class, values.shape[1]), fill, np.float64)
    for position, (_, idx) in enumerate(groups):
        if idx.size:
            out[position] = reducer(values[idx], axis=0)
    return out


def _create_gtiff(tmp_path, width, height, nbands, dtype, geo, projection,
                  nodata, descriptions=None):
    creation = list(GTIFF_OPTIONS)
    if GEOTIFF_COMPRESS:
        creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS,
                         "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL,
                         "INTERLEAVE=BAND"])
        if dtype != gdal.GDT_Byte:
            creation.append("PREDICTOR=3")
    driver = gdal.GetDriverByName("GTiff")
    dataset = driver.Create(tmp_path, width, height, nbands, dtype, options=creation)
    if dataset is None:
        raise RuntimeError("无法创建 GeoTIFF: %s" % tmp_path)
    dataset.SetGeoTransform(geo)
    if projection:
        dataset.SetProjection(projection)
    for index in range(1, nbands + 1):
        band = dataset.GetRasterBand(index)
        band.SetNoDataValue(nodata)
        if descriptions and index <= len(descriptions):
            band.SetDescription(str(descriptions[index - 1]))
    return dataset


def process_area(area_name, src_dir, out_root, library, response_source, args, logger):
    started = time.time()
    mosaic_dir = os.path.join(src_dir, area_name, "mosaic")
    in_tif = os.path.join(mosaic_dir, IN_MOSAIC_NAME)
    in_mask = os.path.join(mosaic_dir, IN_MASK_NAME)
    in_wl = os.path.join(mosaic_dir, IN_WAVELENGTHS_NAME)
    out_dir = os.path.join(out_root, area_name)

    final_class = os.path.join(out_dir, OUT_CLASS_NAME)
    final_info = os.path.join(out_dir, OUT_INFO_NAME)
    if not args.overwrite and os.path.isfile(final_class) and os.path.isfile(final_info):
        logger.info("[%s] 跳过（成果已存在）: %s", area_name, out_dir)
        return {"area": area_name, "status": "跳过", "message": "成果已存在",
                "elapsed": 0.0, "classes": None,
                "valid_pixels": 0, "classified_pixels": 0}

    if not os.path.isfile(in_wl):
        raise FileNotFoundError("缺少波长文件: %s" % in_wl)
    wavelengths = np.load(in_wl).astype(np.float64)
    if wavelengths.ndim != 1:
        raise ValueError("波长文件不是 1 维: %s" % in_wl)

    ds = gdal.Open(in_tif)
    if ds is None:
        raise RuntimeError("无法打开镶嵌反射率: %s" % in_tif)
    if ds.RasterCount != wavelengths.size:
        raise ValueError("镶嵌波段数 %d != 波长数 %d" % (ds.RasterCount, wavelengths.size))
    width, height = ds.RasterXSize, ds.RasterYSize
    geo = ds.GetGeoTransform()
    projection = ds.GetProjection()
    pixel_area_m2 = abs(float(geo[1]) * float(geo[5]))   # 30 m 网格为 900，60 m 网格为 3600

    mask_ds = gdal.Open(in_mask)
    if mask_ds is None:
        raise RuntimeError("无法打开掩膜: %s" % in_mask)
    if (mask_ds.RasterXSize, mask_ds.RasterYSize) != (width, height):
        raise ValueError("掩膜尺寸与镶嵌不一致")

    library_spectra = library["spectra"]
    groups = library["groups"]
    candidate_idx = library["match_idx"]

    # 先统计坏波段（1.44 μm 等水汽带附近，步骤2 有大量 nodata），再定参与匹配的波段
    valid_fraction, _ = estimate_band_missing(ds, mask_ds, candidate_idx, logger)
    keep_mask = valid_fraction >= (1.0 - MAX_BAND_MISSING_FRACTION)
    match_idx = candidate_idx[keep_mask]
    for position, band_index in enumerate(candidate_idx):
        if keep_mask[position]:
            continue
        logger.info("  剔除坏波段 band %d (%.4f μm, 有效 %.1f%%)",
                    band_index + 1, wavelengths[band_index],
                    100.0 * valid_fraction[position])
    if match_idx.size < 100:
        raise ValueError("可用匹配波段只有 %d 个，数据质量异常" % match_idx.size)
    n_class = len(MINERAL_EN)
    n_end = library_spectra.shape[0]

    # 连续统 / 特征窗口
    win_idx = np.where((wavelengths >= ALOH_CONTINUUM_RANGE[0]) &
                       (wavelengths <= ALOH_CONTINUUM_RANGE[1]))[0]
    if win_idx.size < MIN_ALOH_BANDS:
        raise ValueError("Al-OH 窗口内只有 %d 个波段，不足 %d" % (win_idx.size, MIN_ALOH_BANDS))
    win_wl = wavelengths[win_idx]
    left_index, right_index = win_idx[0], win_idx[-1]
    win_t = ((win_wl - wavelengths[left_index]) /
             (wavelengths[right_index] - wavelengths[left_index]))
    win_dwl_nm = float(np.mean(np.diff(win_wl))) * 1000.0
    eval_mask = (win_wl >= ALOH_FEATURE_RANGE[0]) & (win_wl <= ALOH_FEATURE_RANGE[1])
    if eval_mask.sum() < 3:
        raise ValueError("Al-OH 特征评价范围内波段数不足")
    win_wl_eval = win_wl[eval_mask]

    # 库的连续统去除光谱（同一窗口 + 同一锚点索引）
    lib_win = library_spectra[:, win_idx].astype(np.float64)
    lib_left = lib_win[:, 0:1]
    lib_right = lib_win[:, -1:]
    lib_cont = lib_left + (lib_right - lib_left) * win_t[None, :]
    if np.any(lib_cont <= 0):
        logger.warning("%d 个库端元连续统为负，可能异常",
                       int((lib_cont <= 0).any(axis=1).sum()))
    lib_cr = lib_win / np.maximum(lib_cont, 1e-6)
    lib_cr_eval = lib_cr[:, eval_mask]
    lib_cr_sq = np.sum(lib_cr_eval * lib_cr_eval, axis=1)

    # 匹配用库光谱：MF 用原始反射率；SAM 用连续统去除后的 Al-OH 窗口
    lib_match = library_spectra[:, match_idx].astype(np.float64)
    lib_cr_unit = lib_cr_eval / np.maximum(np.sqrt(lib_cr_sq)[:, None], 1e-9)

    # 背景统计
    mean, inv_cov = sample_background(ds, mask_ds, match_idx, BG_SAMPLE_PIXELS, logger)
    mf_q = inv_cov @ (lib_match - mean).T                            # (151, n_end)
    mf_denom = np.sum((lib_match - mean) * mf_q.T, axis=1)           # (n_end,)

    logger.info("[%s] %d×%d, %d 波段，端元 %d，类别 %d；匹配波段 %d",
                area_name, width, height, wavelengths.size, n_end, n_class, match_idx.size)
    logger.info("[%s] Al-OH 窗口: %.3f~%.3f μm (%d 波段, 宽 %.1f nm/波段)",
                area_name, win_wl.min(), win_wl.max(), win_wl.size, win_dwl_nm)

    os.makedirs(out_dir, exist_ok=True)
    created = {}

    def make_output(name, nbands, dtype, nodata, descriptions=None):
        tmp_path = os.path.join(out_dir, name + ".tmp.tif")
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        dataset = _create_gtiff(tmp_path, width, height, nbands, dtype, geo,
                                projection, nodata, descriptions)
        created[name] = (dataset, tmp_path, os.path.join(out_dir, name))
        return dataset

    descs = ["%s %s" % (name, MINERAL_CN[name]) for name in MINERAL_EN]
    out_class = make_output(OUT_CLASS_NAME, 1, gdal.GDT_Byte, CLASS_NODATA,
                            ["Al-OH mineral class"])
    out_confidence = make_output(OUT_CONFIDENCE_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                                 ["classification confidence (%)"])
    out_depth = make_output(OUT_DEPTH_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                            ["Al-OH 2.2um absorption depth"])
    out_position = make_output(OUT_POSITION_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                               ["Al-OH absorption center (um)"])
    out_width = make_output(OUT_WIDTH_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                            ["Al-OH half-depth width (nm)"])
    out_sam = make_output(OUT_SAM_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)
    out_sff = make_output(OUT_SFF_ABUNDANCE_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)
    out_rms = make_output(OUT_SFF_RMS_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)
    out_mf = out_mtmf = None
    if ENABLE_MF and not args.no_mf:
        out_mf = make_output(OUT_MF_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)
        out_mtmf = make_output(OUT_MTMF_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)

    # 统计累加器
    class_pixels = np.zeros(n_class + 1, np.int64)
    class_depth_sum = np.zeros(n_class + 1, np.float64)
    class_conf_sum = np.zeros(n_class + 1, np.float64)
    class_angle_sum = np.zeros(n_class + 1, np.float64)
    valid_pixels = 0
    clear_pixels = 0
    classified_pixels = 0
    depth_values = []
    angle_values = []

    def write_block(dataset, band_index, full, idx, values, rows, y0):
        block = full.copy()
        block[idx] = values
        block = block.reshape(rows, width)
        dataset.GetRasterBand(band_index).WriteArray(block, 0, y0)

    try:
        for y0 in range(0, height, BLOCK_ROWS):
            rows = min(BLOCK_ROWS, height - y0)
            cube = ds.ReadAsArray(0, y0, width, rows)
            if cube is None:
                raise RuntimeError("读取镶嵌失败: y=%d h=%d" % (y0, rows))
            mask = mask_ds.ReadAsArray(0, y0, width, rows)
            if mask.ndim == 3:
                mask = mask[0]
            clear_pixels += int(np.count_nonzero(mask == 1))
            if USE_CLOUD_FILLED:
                pixel_ok = (mask == 1) | (mask == 2) | (mask == 3)
            else:
                pixel_ok = (mask == 1)
            npix = width * rows
            pixel_ok = pixel_ok.reshape(-1)

            match_block = cube[match_idx].reshape(match_idx.size, -1).astype(np.float64)
            win_block = cube[win_idx].reshape(win_idx.size, -1).astype(np.float64)
            match_finite = _valid_pixels(match_block)
            missing = np.count_nonzero(~match_finite, axis=0)
            low = np.where(match_finite, match_block, np.inf).min(axis=0)
            high = np.where(match_finite, match_block, -np.inf).max(axis=0)
            pixel_ok &= (missing <= MAX_MISSING_BANDS_MATCH)
            pixel_ok &= (low > REFL_MIN) & (high < REFL_MAX)
            pixel_ok &= _valid_pixels(win_block).all(axis=0)
            idx = np.where(pixel_ok)[0]
            if idx.size == 0:
                continue

            # 缺测波段按 0 填充（只影响 151 个匹配波段里的 1~2 个，对 SAM/MF 影响很小）
            match_filled = np.where(match_finite, match_block, 0.0)

            # 连续统去除
            win_valid = win_block[:, idx]
            left = win_valid[0][None, :]
            right = win_valid[-1][None, :]
            continuum = left + (right - left) * win_t[:, None]
            cr_ok = continuum.min(axis=0) > MIN_CONTINUUM_VALUE
            idx = idx[cr_ok]
            if idx.size == 0:
                continue
            match_valid = match_filled[:, idx]
            win_valid = win_valid[:, cr_ok]
            continuum = continuum[:, cr_ok]
            cr = win_valid / continuum

            depth_eval = 1.0 - cr[eval_mask]
            depth = np.clip(depth_eval.max(axis=0), -1.0, 1.0)
            argmax = depth_eval.argmax(axis=0)
            position = win_wl_eval[argmax]
            half = depth_eval >= (0.5 * depth[None, :])
            width_nm = half.sum(axis=0) * win_dwl_nm

            # ── SAM：连续统去除后的 Al-OH 窗口光谱角（消除反照率/大气残差影响） ──
            crn = cr[eval_mask]
            crn_norm = np.linalg.norm(crn, axis=0)
            good_norm = crn_norm > 1e-9
            crn_unit = crn / np.maximum(crn_norm[None, :], 1e-9)
            cos = lib_cr_unit @ crn_unit
            np.clip(cos, -1.0, 1.0, out=cos)
            angle_end = np.degrees(np.arccos(cos))                 # (n_end, nv)

            # ── SFF（连续统去除光谱特征拟合） ──
            cr_eval = cr[eval_mask]
            dot = cr_eval.T @ lib_cr_eval.T                        # (nv, n_end)
            scale = dot / np.maximum(lib_cr_sq[None, :], 1e-12)    # (nv, n_end)
            cr_sq = np.sum(cr_eval * cr_eval, axis=0)              # (nv,)
            resid = cr_sq[:, None] - (scale ** 2) * lib_cr_sq[None, :]
            rms = np.sqrt(np.maximum(resid, 0.0) / lib_cr_eval.shape[1])

            # ── 类别级汇总 ──
            angle_class = _class_reduce(angle_end, groups, np.min, np.inf)
            rms_class = _class_reduce(rms.T, groups, np.min, np.inf)
            abundance_class = _class_reduce(scale.T, groups, np.max, 0.0)
            abundance_class = np.clip(abundance_class, 0.0, 1.0)

            best_index = np.argmin(angle_class, axis=0)
            arange = np.arange(angle_class.shape[1])
            best_angle = angle_class[best_index, arange]
            best_rms = rms_class[best_index, arange]
            best_abundance = abundance_class[best_index, arange]

            if args.class_by == "sff":
                best_index = np.argmin(rms_class, axis=0)
                best_angle = angle_class[best_index, arange]
                best_rms = rms_class[best_index, arange]
                best_abundance = abundance_class[best_index, arange]

            usable = (depth >= args.min_depth) & good_norm
            classify = usable & (best_angle <= args.max_angle) & (best_rms <= args.max_rms)
            class_values = np.where(classify, best_index + 1, 0).astype(np.uint8)
            confidence = np.where(usable,
                                  np.clip(100.0 * (1.0 - best_angle / max(args.max_angle, 1e-6)),
                                          0.0, 100.0), 0.0)

            # ── MF / MTMF ──
            mf_class = infeas_class = None
            if out_mf is not None:
                centered = (match_valid - mean[:, None]).T          # (nv, 151)
                mf_score = centered @ mf_q / np.maximum(mf_denom[None, :], 1e-12)
                mf_class = _class_reduce(mf_score.T, groups, np.max, 0.0).T
                maha_sq = np.maximum(np.sum((centered @ inv_cov) * centered, axis=1), 0.0)
                infeas = np.sqrt(np.maximum(
                    maha_sq[:, None] - (mf_score ** 2) * mf_denom[None, :], 0.0))
                infeas_class = _class_reduce(infeas.T, groups, np.min, np.inf).T

            # ── 写块 ──
            class_block = np.full(npix, CLASS_NODATA, np.uint8)
            class_block[idx] = class_values
            out_class.GetRasterBand(1).WriteArray(class_block.reshape(rows, width), 0, y0)

            for dataset, values, nodata in (
                    (out_confidence, confidence, NODATA_OUT),
                    (out_depth, depth, NODATA_OUT)):
                buffer = np.full(npix, nodata, np.float32)
                buffer[idx] = np.where(np.isfinite(values), values, nodata)
                dataset.GetRasterBand(1).WriteArray(buffer.reshape(rows, width), 0, y0)

            # 吸收位置/宽度只在真正有吸收（深度≥阈值）的像元上输出，避免“无特征”像元
            # 在窗口边缘取值造成整幅图都是 2.28 μm 的假象
            feature_ok = usable & (depth >= args.min_depth)
            for dataset, values in ((out_position, position), (out_width, width_nm)):
                buffer = np.full(npix, NODATA_OUT, np.float32)
                buffer[idx] = np.where(feature_ok & np.isfinite(values),
                                       values, NODATA_OUT)
                dataset.GetRasterBand(1).WriteArray(buffer.reshape(rows, width), 0, y0)

            for band_index in range(n_class):
                buffer = np.full(npix, NODATA_OUT, np.float32)
                buffer[idx] = np.where(np.isfinite(angle_class[band_index]),
                                       angle_class[band_index], NODATA_OUT)
                out_sam.GetRasterBand(band_index + 1).WriteArray(
                    buffer.reshape(rows, width), 0, y0)
                buffer = np.full(npix, NODATA_OUT, np.float32)
                buffer[idx] = np.where(np.isfinite(abundance_class[band_index]),
                                       abundance_class[band_index], NODATA_OUT)
                out_sff.GetRasterBand(band_index + 1).WriteArray(
                    buffer.reshape(rows, width), 0, y0)
                buffer = np.full(npix, NODATA_OUT, np.float32)
                buffer[idx] = np.where(np.isfinite(rms_class[band_index]),
                                       rms_class[band_index], NODATA_OUT)
                out_rms.GetRasterBand(band_index + 1).WriteArray(
                    buffer.reshape(rows, width), 0, y0)
                if out_mf is not None:
                    buffer = np.full(npix, NODATA_OUT, np.float32)
                    buffer[idx] = np.where(np.isfinite(mf_class[:, band_index]),
                                           mf_class[:, band_index], NODATA_OUT)
                    out_mf.GetRasterBand(band_index + 1).WriteArray(
                        buffer.reshape(rows, width), 0, y0)
                    buffer = np.full(npix, NODATA_OUT, np.float32)
                    buffer[idx] = np.where(np.isfinite(infeas_class[:, band_index]),
                                           infeas_class[:, band_index], NODATA_OUT)
                    out_mtmf.GetRasterBand(band_index + 1).WriteArray(
                        buffer.reshape(rows, width), 0, y0)

            # ── 统计 ──
            class_pixels += np.bincount(class_values.astype(np.int64), minlength=n_class + 1)
            valid_pixels += int(idx.size)
            classified_pixels += int(np.count_nonzero(classify))
            for class_index in range(n_class + 1):
                mask_class = class_values == class_index
                if np.any(mask_class):
                    class_depth_sum[class_index] += float(depth[mask_class].sum())
                    class_conf_sum[class_index] += float(confidence[mask_class].sum())
                    class_angle_sum[class_index] += float(best_angle[mask_class].sum())
            depth_values.append(depth.astype(np.float32))
            angle_values.append(best_angle.astype(np.float32))

            if (y0 // BLOCK_ROWS) % 4 == 0 or y0 + rows >= height:
                logger.info("  [%s] 行 %d/%d（%.0f%%），有效像元累计 %d",
                            area_name, y0 + rows, height,
                            100.0 * (y0 + rows) / height, valid_pixels)
    except Exception:
        out_class = out_confidence = out_depth = out_position = out_width = None
        out_sam = out_sff = out_rms = out_mf = out_mtmf = None
        for name in list(created):
            dataset, tmp_path, final_path = created.pop(name)
            try:
                dataset.FlushCache()
            except Exception:
                pass
            dataset = None
            try:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)
            except OSError:
                pass
        raise
    finally:
        ds = None
        mask_ds = None

    # 关闭并原子替换（先清掉局部变量引用，Windows 下文件句柄才能释放）
    n_written = len(created)
    out_class = out_confidence = out_depth = out_position = out_width = None
    out_sam = out_sff = out_rms = out_mf = out_mtmf = None
    for name in list(created):
        dataset, tmp_path, final_path = created.pop(name)
        dataset.FlushCache()
        dataset = None
        os.replace(tmp_path, final_path)

    logger.info("[%s] 写出 %d 个成果文件", area_name, n_written)

    # 汇总/图例/信息
    rows = []
    for class_index in range(n_class + 1):
        count = int(class_pixels[class_index])
        if class_index == 0:
            name_en, name_cn = "unclassified", "未分类"
        else:
            name_en = MINERAL_EN[class_index - 1]
            name_cn = MINERAL_CN[name_en]
        mean_depth = (class_depth_sum[class_index] / count) if count else 0.0
        mean_conf = (class_conf_sum[class_index] / count) if count else 0.0
        mean_angle = (class_angle_sum[class_index] / count) if count else 0.0
        rows.append({
            "class_id": class_index,
            "mineral_en": name_en,
            "mineral_cn": name_cn,
            "pixels": count,
            "area_km2": count * pixel_area_m2 / 1.0e6,
            "percent_of_valid": (100.0 * count / valid_pixels) if valid_pixels else 0.0,
            "mean_depth": mean_depth,
            "mean_confidence": mean_conf,
            "mean_sam_deg": mean_angle,
        })
    write_class_summary(os.path.join(out_dir, OUT_SUMMARY_NAME), rows)
    write_legend(os.path.join(out_dir, OUT_LEGEND_NAME), library)

    depth_stats = np.concatenate(depth_values) if depth_values else np.array([])
    angle_stats = np.concatenate(angle_values) if angle_values else np.array([])
    try:
        srs = osr.SpatialReference(wkt=projection) if projection else None
        epsg = int(srs.GetAuthorityCode(None) or 0) if srs is not None else 0
    except Exception:
        epsg = 0
    info = {
        "area": area_name,
        "source": in_tif,
        "grid": {"width": int(width), "height": int(height),
                 "geotransform": [float(v) for v in geo], "epsg": epsg},
        "bands": int(wavelengths.size),
        "matching_bands": int(match_idx.size),
        "sensor_response": response_source,
        "library": {"root": os.path.abspath(args.lib_root),
                    "endmembers": int(n_end), "classes": n_class},
        "parameters": {
            "aloh_continuum_range_um": list(ALOH_CONTINUUM_RANGE),
            "aloh_feature_range_um": list(ALOH_FEATURE_RANGE),
            "min_aloh_depth": float(args.min_depth),
            "sam_max_deg": float(args.max_angle),
            "sff_max_rms": float(args.max_rms),
            "class_by": args.class_by,
            "mf_enabled": bool(out_mf is not None),
            "use_cloud_filled": bool(USE_CLOUD_FILLED),
        },
        "coverage": {
            "clear_pixels": int(clear_pixels),
            "valid_pixels": int(valid_pixels),
            "classified_pixels": int(classified_pixels),
            "classified_percent_of_valid": (100.0 * classified_pixels / valid_pixels) if valid_pixels else 0.0,
        },
        "statistics": {
            "depth_percentiles": _percentiles(depth_stats, [1, 5, 50, 95, 99]),
            "sam_percentiles": _percentiles(angle_stats, [1, 5, 50, 95, 99]),
        },
        "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(final_info + ".tmp", "w", encoding="utf-8") as handle:
        json.dump(info, handle, ensure_ascii=False, indent=2)
    os.replace(final_info + ".tmp", final_info)

    if MAKE_QC:
        try:
            make_area_qc(out_dir, area_name, in_tif, wavelengths, rows, logger)
        except Exception as exc:
            logger.warning("[%s] 质检图失败: %s", area_name, exc, exc_info=True)

    elapsed = time.time() - started
    logger.info("[%s] 完成: 有效像元 %d，分类像元 %d（%.1f%% of valid），用时 %.1f s",
                area_name, valid_pixels, classified_pixels,
                100.0 * classified_pixels / max(1, valid_pixels), elapsed)
    for row in rows[1:]:
        if row["pixels"]:
            logger.info("    %-18s %8d 像元  %8.2f km2  %6.2f%%  深度 %.3f  SAM %.1f°",
                        row["mineral_en"], row["pixels"], row["area_km2"],
                        row["percent_of_valid"], row["mean_depth"], row["mean_sam_deg"])
    return {"area": area_name, "status": "成功",
            "message": "分类像元 %d" % classified_pixels,
            "elapsed": elapsed, "classes": rows,
            "valid_pixels": int(valid_pixels),
            "classified_pixels": int(classified_pixels)}


def _percentiles(values, percents):
    if values.size == 0:
        return [None for _ in percents]
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return [None for _ in percents]
    return [float(v) for v in np.percentile(finite, percents)]


def write_class_summary(path, rows):
    fields = ["class_id", "mineral_en", "mineral_cn", "pixels", "area_km2",
              "percent_of_valid", "mean_depth", "mean_confidence", "mean_sam_deg"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            out = dict(row)
            out["area_km2"] = "%.4f" % out["area_km2"]
            out["percent_of_valid"] = "%.3f" % out["percent_of_valid"]
            out["mean_depth"] = "%.5f" % out["mean_depth"]
            out["mean_confidence"] = "%.2f" % out["mean_confidence"]
            out["mean_sam_deg"] = "%.3f" % out["mean_sam_deg"]
            writer.writerow(out)


def write_legend(path, library):
    class_ids = library["class_ids"]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["class_id", "mineral_en", "mineral_cn", "color",
                         "endmembers", "samples"])
        for class_index, name_en in enumerate(MINERAL_EN, start=1):
            idx = np.where(class_ids == class_index)[0]
            samples = ";".join(library["meta"][i]["sample"] for i in idx)
            writer.writerow([class_index, name_en, MINERAL_CN[name_en],
                             MINERAL_COLOR[name_en], int(idx.size), samples])


def make_area_qc(out_dir, area_name, mosaic_path, wavelengths, rows, logger):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    class_path = os.path.join(out_dir, OUT_CLASS_NAME)
    depth_path = os.path.join(out_dir, OUT_DEPTH_NAME)
    pos_path = os.path.join(out_dir, OUT_POSITION_NAME)
    conf_path = os.path.join(out_dir, OUT_CONFIDENCE_NAME)
    sam_path = os.path.join(out_dir, OUT_SAM_NAME)

    dataset = gdal.Open(class_path)
    if dataset is None:
        raise RuntimeError("无法打开分类成果: %s" % class_path)
    width, height = dataset.RasterXSize, dataset.RasterYSize
    geo = dataset.GetGeoTransform()
    dataset = None
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))
    bw, bh = max(1, width // scale), max(1, height // scale)

    def read_band(path, band=1, dtype=np.float64):
        data = gdal.Open(path)
        if data is None:
            return None
        if data.RasterCount == 1:
            array = data.ReadAsArray(buf_xsize=bw, buf_ysize=bh)
        else:
            band_handle = data.GetRasterBand(band)
            array = band_handle.ReadAsArray(buf_xsize=bw, buf_ysize=bh)
            band_handle = None
        data = None
        return array.astype(dtype)

    class_array = read_band(class_path, 1, np.int32)
    depth_array = read_band(depth_path)
    pos_array = read_band(pos_path)
    conf_array = read_band(conf_path)

    sam_data = gdal.Open(sam_path)
    sam_best = None
    if sam_data is not None:
        sam_best = np.full((bh, bw), np.inf, np.float64)
        for band_index in range(1, sam_data.RasterCount + 1):
            band_handle = sam_data.GetRasterBand(band_index)
            band = band_handle.ReadAsArray(buf_xsize=bw, buf_ysize=bh).astype(np.float64)
            band_handle = None
            band = np.where(band < NODATA_THRESHOLD, np.inf, band)
            sam_best = np.minimum(sam_best, band)
        sam_data = None
        sam_best[~np.isfinite(sam_best)] = np.nan

    # 真彩色
    color_bands = [_nearest_band(wavelengths, TRUE_COLOR_TARGETS[key]) + 1
                   for key in ("R", "G", "B")]
    rgb = None
    mosaic = gdal.Open(mosaic_path)
    if mosaic is not None:
        channels = []
        for band_index in color_bands:
            band_handle = mosaic.GetRasterBand(band_index)
            channels.append(band_handle.ReadAsArray(
                buf_xsize=bw, buf_ysize=bh).astype(np.float64))
            band_handle = None
        mosaic = None
        valid = np.ones(channels[0].shape, bool)
        for channel in channels:
            valid &= np.isfinite(channel) & (channel > -9000)
        rgb = np.ones(channels[0].shape + (3,), np.float64)
        for position, channel in enumerate(channels):
            rgb[..., position] = np.where(valid, _stretch(channel, valid), 1.0)

    colors = ["#dddddd"] + [MINERAL_COLOR[name] for name in MINERAL_EN]
    cmap = ListedColormap(colors)
    bounds = np.arange(-0.5, len(colors) + 0.5, 1.0)
    extent = [geo[0] / 1000.0, (geo[0] + geo[1] * width) / 1000.0,
              (geo[3] + geo[5] * height) / 1000.0, geo[3] / 1000.0]
    class_display = np.where(class_array == CLASS_NODATA, 0, class_array)

    figure, axes = plt.subplots(2, 3, figsize=(18, 10.5))
    axes = axes.ravel()

    if rgb is not None:
        axes[0].imshow(rgb, extent=extent, origin="upper", interpolation="nearest")
    axes[0].set_title("True color (R %.2f / G %.2f / B %.2f um)"
                      % tuple(wavelengths[band - 1] for band in color_bands), fontsize=10)

    axes[1].imshow(class_display, extent=extent, origin="upper", cmap=cmap,
                   norm=BoundaryNorm(bounds, cmap.N), interpolation="nearest")
    axes[1].set_title("Al-OH mineral class", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[index])
               for index in range(1, len(colors))]
    axes[1].legend(handles, MINERAL_EN, fontsize=6.5, ncol=2, loc="lower left",
                   framealpha=0.85)

    depth_display = np.where(depth_array > -9000, depth_array, np.nan)
    image = axes[2].imshow(depth_display, extent=extent, origin="upper", cmap="viridis",
                           vmin=-0.5, vmax=0.5, interpolation="nearest")
    axes[2].set_title("2.2 um absorption depth (clip -0.5~0.5)", fontsize=10)
    figure.colorbar(image, ax=axes[2], fraction=0.035, pad=0.02)

    pos_display = np.where(pos_array > -9000, pos_array, np.nan)
    image = axes[3].imshow(pos_display, extent=extent, origin="upper", cmap="turbo",
                           vmin=float(ALOH_FEATURE_RANGE[0]),
                           vmax=float(ALOH_FEATURE_RANGE[1]), interpolation="nearest")
    axes[3].set_title("Absorption center (um)", fontsize=10)
    figure.colorbar(image, ax=axes[3], fraction=0.035, pad=0.02)

    conf_display = np.where(conf_array > -9000, conf_array, np.nan)
    image = axes[4].imshow(conf_display, extent=extent, origin="upper", cmap="magma",
                           vmin=0, vmax=100, interpolation="nearest")
    axes[4].set_title("Classification confidence (%)", fontsize=10)
    figure.colorbar(image, ax=axes[4], fraction=0.035, pad=0.02)

    labels = [row["mineral_en"] for row in rows[1:]]
    values = [row["area_km2"] for row in rows[1:]]
    order = np.argsort(values)[::-1]
    axes[5].barh([labels[i] for i in order], [values[i] for i in order],
                 color=[MINERAL_COLOR[labels[i]] for i in order])
    axes[5].set_xlabel("Area (km2)")
    axes[5].set_title("Al-OH class area", fontsize=10)
    axes[5].tick_params(labelsize=7.5)
    axes[5].grid(axis="x", alpha=0.3)

    for axis in axes[:5]:
        axis.set_xlabel("Easting (km)", fontsize=8)
        axis.set_ylabel("Northing (km)", fontsize=8)
        axis.tick_params(labelsize=7)

    try:
        depth_p50 = float(np.nanpercentile(depth_display, 50))
        angle_p50 = float(np.nanpercentile(sam_best, 50))
    except Exception:
        depth_p50 = angle_p50 = float("nan")
    figure.suptitle("%s | Al-OH extraction | depth p50 %.3f | SAM p50 %.1f deg"
                    % (area_name, depth_p50, angle_p50), fontsize=12)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(out_dir, OUT_QC_NAME)
    figure.savefig(path, dpi=QC_DPI, bbox_inches="tight")
    plt.close(figure)
    logger.info("[%s] 质检图: %s", area_name, path)


# ╔══════════════════════════════════════════════════════════════╗
# ║  主流程                                                       ║
# ╚══════════════════════════════════════════════════════════════╝

def write_areas_summary(path, results):
    fields = ["area", "status", "valid_pixels", "classified_pixels",
              "dominant_mineral", "dominant_percent", "message", "elapsed"]
    existing = {}
    if os.path.isfile(path):
        try:
            with open(path, "r", encoding="utf-8-sig", newline="") as handle:
                for row in csv.DictReader(handle):
                    existing[row.get("area", "")] = row
        except OSError:
            existing = {}
    for result in results:
        classes = result.get("classes") or []
        dominant = None
        if len(classes) > 1:
            dominant = max(classes[1:], key=lambda row: row["pixels"])
        existing[result["area"]] = {
            "area": result["area"],
            "status": result["status"],
            "valid_pixels": str(result.get("valid_pixels", "")),
            "classified_pixels": str(result.get("classified_pixels", "")),
            "dominant_mineral": dominant["mineral_en"] if dominant else "",
            "dominant_percent": ("%.2f" % dominant["percent_of_valid"]) if dominant else "",
            "message": result.get("message", ""),
            "elapsed": "%.1f" % result.get("elapsed", 0.0),
        }
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for area in sorted(existing):
            writer.writerow({key: existing[area].get(key, "") for key in fields})


def main(argv=None):
    args = parse_args(argv)
    out_root = os.path.abspath(args.out_root)
    os.makedirs(out_root, exist_ok=True)
    catalog_dir = os.path.join(out_root, CATALOG_DIR_NAME)
    os.makedirs(catalog_dir, exist_ok=True)
    log_path = args.log or os.path.join(catalog_dir, "step4_batch.log")
    logger = setup_logging(log_path, args.verbose)
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)

    logger.info("=" * 72)
    logger.info("ZY1E AHSI 步骤4：Al-OH 蚀变矿物提取（splib07a + SAM/SFF/MF）")
    logger.info("步骤3 目录: %s", os.path.abspath(args.in_root))
    logger.info("输出目录  : %s", out_root)
    logger.info("光谱库    : %s", args.lib_root)
    logger.info("Al-OH 类  : %s", ", ".join(MINERAL_EN))
    logger.info("日志      : %s", log_path)
    logger.info("=" * 72)

    try:
        areas = discover_areas(args.in_root, args.area)
    except Exception as exc:
        logger.error("扫描步骤3 成果失败: %s", exc)
        return 2
    if args.limit > 0:
        areas = areas[:args.limit]
    if not areas:
        logger.error("在 %s 下未找到可处理的区域（需要 <区名>/mosaic/%s）",
                     args.in_root, IN_MOSAIC_NAME)
        return 2
    logger.info("待处理区域: %d (%s)", len(areas), ", ".join(areas))

    first_mosaic = os.path.join(args.in_root, areas[0], "mosaic")
    wavelengths = np.load(os.path.join(first_mosaic, IN_WAVELENGTHS_NAME)).astype(np.float64)
    band_fwhm, response_source = load_sensor_fwhm(TAR_DIR, wavelengths, logger)

    if args.dry_run:
        logger.info("-" * 72)
        logger.info("dry-run：不执行实际处理")
        library_dir = os.path.join(out_root, LIB_DIR_NAME)
        logger.info("匹配库目录: %s（%s）", library_dir,
                    "已存在，将复用" if os.path.isfile(os.path.join(library_dir, "ahsi_library.npy"))
                    else "尚未构建，正式运行时构建")
        for name in areas:
            logger.info("  %-10s <- %s", name, os.path.join(args.in_root, name, "mosaic"))
        logger.info("阈值: min_depth=%.3f, max_angle=%.1f°, max_rms=%.3f, class_by=%s",
                    args.min_depth, args.max_angle, args.max_rms, args.class_by)
        return 0

    library_dir = os.path.join(out_root, LIB_DIR_NAME)
    library = build_library(args.lib_root, wavelengths, band_fwhm, library_dir, logger,
                            force=args.rebuild_library)

    if args.library_only:
        logger.info("-" * 72)
        logger.info("匹配库: %d 端元 × %d 波段", library["spectra"].shape[0],
                    library["spectra"].shape[1])
        for class_index, name_en in enumerate(MINERAL_EN, start=1):
            count = int(np.count_nonzero(library["class_ids"] == class_index))
            logger.info("  %-16s %d 个端元", name_en, count)
        logger.info("库目录: %s", library_dir)
        return 0

    results = []
    started = time.time()
    for index, area in enumerate(areas, start=1):
        logger.info("-" * 72)
        logger.info("(%d/%d) %s", index, len(areas), area)
        try:
            results.append(process_area(area, os.path.abspath(args.in_root), out_root,
                                        library, response_source, args, logger))
        except Exception as exc:
            logger.exception("[%s] 处理失败", area)
            results.append({"area": area, "status": "失败", "message": str(exc),
                            "elapsed": 0.0, "classes": None,
                            "valid_pixels": 0, "classified_pixels": 0})

    elapsed = time.time() - started
    logger.info("=" * 72)
    success = [r for r in results if r["status"] == "成功"]
    skipped = [r for r in results if r["status"] == "跳过"]
    failed = [r for r in results if r["status"] == "失败"]
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %.1f s",
                len(success), len(skipped), len(failed), elapsed)
    for result in results:
        logger.info("  [%s] %s - %s", result["status"], result["area"],
                    result.get("message", ""))
    if results:
        write_areas_summary(os.path.join(catalog_dir, "areas_summary.csv"), results)
    if failed:
        return 1
    logger.info("成果目录: %s", out_root)
    return 0


if __name__ == "__main__":
    sys.exit(main())
