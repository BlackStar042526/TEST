# -*- coding: utf-8 -*-
"""
04_aloh_mineral_extraction.py
ZY1E AHSI 步骤4：云母-绿泥石等蚀变矿物提取（USGS splib07a 库 + CR-SAM / SFF / MF，v1.1.5）

输入：步骤3 每区目录下的 *_clip.tif（逐景裁剪，首次运行自动合并为 VRT），或已镶嵌的
      mosaic/{mosaic_reflectance.tif, mosaic_mask.tif, wavelengths_clean.npy}。
流程：①输入校验与波段选择 ②Al-OH 连续统窗口与背景统计 ③逐像元连续统去除、
      CR-SAM / SFF / MF 分类 ④写栅格、统计与质检图（各阶段日志打印耗时）。
输出：<OUT_ROOT>/<区>/ 下分类图、SAM 角度、SFF 丰度、吸收深度、分类统计 CSV、
      参数 JSON 与质检图（文件名见配置区 OUT_* 常量）；库缓存在 <OUT_ROOT>/_library/。
用法：python 04_aloh_mineral_extraction_v1.1.5.py --dry-run
依赖：GDAL >= 3.0、numpy；matplotlib 可选（质检图）。
v1.1.5：将目标矿物由 11 种 Al-OH 精简为 6 种（云母族 + 绿泥石族 + 蒙脱石），剔除黄玉等非目标矿物；算法与输出文件名不变。
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


# =================================================================
# Windows / conda 的 GDAL DLL 引导
# =================================================================
# Python 3.8+ 的扩展模块不再用 PATH 查找依赖 DLL，需要显式
# os.add_dll_directory。conda 环境的 gdal.dll 在 Library\bin；
# 个别环境因 conda 事务中断，sqlite3.dll 只留在 pkgs 缓存里，
# 这里一并补进 DLL 搜索路径（不影响正常环境）。
def _bootstrap_gdal_dll():
    if os.name != "nt":                             # 仅 Windows 需要手动补 DLL 路径
        return                                      # 其他系统直接返回
    prefixes = []                                   # 候选环境前缀列表
    for key in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "VIRTUAL_ENV"):  # 常见环境变量名
        value = os.environ.get(key)                  # 读取环境变量值
        if value and os.path.isdir(value):           # 值非空且目录存在
            prefixes.append(value)                   # 记为候选前缀
    prefixes.append(os.path.dirname(sys.executable))  # 当前解释器所在目录
    prefixes += [r"C:\ProgramData\anaconda3", r"D:\ProgramData\anaconda3"]  # 补充默认安装路径
    seen = set()                                     # 已处理目录集合，避免重复
    for prefix in prefixes:                          # 逐个前缀查找 DLL 目录
        if not prefix or not os.path.isdir(prefix):  # 前缀无效则跳过
            continue
        candidates = [os.path.join(prefix, "Library", "bin")]  # conda 标准 DLL 目录
        candidates += glob.glob(os.path.join(prefix, "pkgs", "libsqlite-*", "Library", "bin"))  # pkgs 缓存
        for path in candidates:                      # 逐个候选目录
            if path in seen or not os.path.isdir(path):  # 重复或不存在则跳过
                continue
            seen.add(path)                           # 记录已处理
            try:
                os.add_dll_directory(path)           # 加入 DLL 搜索路径
            except OSError:
                pass                                 # 旧版 Python 无此接口，忽略


_bootstrap_gdal_dll()

try:
    from osgeo import gdal, osr                    # GDAL 栅格与投影模块
except ImportError as exc:                         # 解释器缺少 GDAL 时
    raise SystemExit(
        "当前 Python 解释器缺少 GDAL（osgeo）：%s\n"
        "请用装好 GDAL 的 Anaconda 解释器运行本脚本，例如：\n"
        "    E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe 04_aloh_mineral_extraction.py"
        % exc)                                     # 给出安装提示后退出

gdal.UseExceptions()                               # GDAL 出错即抛 Python 异常
gdal.PushErrorHandler("CPLQuietErrorHandler")      # 屏蔽 GDAL 非致命错误输出


# =================================================================
# 配置区（按需修改）
# =================================================================

# 输入/输出：步骤3 成果目录、本步成果根目录
IN_ROOT_STEP3 = r"E:\ZhouYang\test\ZY1E_output_step3"    # 步骤3 裁剪成果根目录
OUT_ROOT_STEP4 = r"E:\ZhouYang\test\ZY1E_output_step4"   # 本步成果根目录

# 步骤3 直接产物（按研究区裁剪的逐景反射率）后缀；旧后缀自动兼容
STEP3_CLIP_SUFFIXES = ("_cal_badband_ortho_atmcorr_clip.tif",   # 新命名后缀
                       "_ortho_atmcorr_clip.tif")               # 旧命名后缀
# 波长来源：优先步骤1，其次步骤2（<场景名>/wavelengths_clean.npy）
STEP1_ROOT = r"E:\ZhouYang\test\ZY1E_output_step1"       # 步骤1 成果目录
STEP2_ROOT = r"E:\ZhouYang\test\ZY1E_output_step2"       # 步骤2 成果目录
# 裁剪成果合并后的 VRT/掩膜缓存目录（位于本步输出根目录下）
CACHE_DIR_NAME = "_cache"                                # 缓存子目录名
CACHE_VRT_NAME = "mosaic_reflectance.vrt"                # 合并 VRT 文件名
MASK_SAMPLE_BANDS = 12         # 生成掩膜时抽样检查的波段数

# 输入：USGS splib07a ASCII 光谱库目录
LIB_ROOT = r"E:\ZhouYang\test\ASCIIdata_splib07a"        # 光谱库默认根目录
LIB_CHAPTER_DIR = "ChapterM_Minerals"                    # 矿物章子目录名

# 输入：原始 tar 目录（用于读 AHSI 光谱响应 FWHM）
TAR_DIR = r"E:\ZhouYang\test\ZY1E_data"                  # 原始数据 tar 目录
RESPONSE_SUFFIXES = ("_VNIR_Spectralresponse.raw", "_SWIR_Spectralresponse.raw")  # 响应曲线文件后缀

# 参数：传感器响应（读不到 raw 时的回退值）
VNIR_FWHM_NM_FALLBACK = 8.76       # VNIR 回退 FWHM（nm）
SWIR_FWHM_NM_FALLBACK = 16.26      # SWIR 回退 FWHM（nm）
VNIR_SWIR_SPLIT_UM = 1.05      # 回退时按该波长分 VNIR/SWIR

# 参数：目标矿物清单（库记录首词必须精确等于英文名；顺序决定类别号）
# 值：(中文名, 绘图颜色)
TARGET_MINERALS = [
    ("Muscovite",       "白云母/绢云母", "#e76f51"),   # 云母族：白云母（绢云母为其细粒变种），Al-OH 2.20μm
    ("Illite",          "伊利石",       "#2a9d8f"),     # 云母族：伊利石，Al-OH 2.20μm
    ("Chlorite",        "绿泥石",       "#2a6f97"),     # 绿泥石族：泛指，Fe-OH 2.25μm
    ("Clinochlore",     "镁绿泥石",     "#61a5c2"),     # 绿泥石族：富镁，Mg-OH 2.33μm（落在匹配窗口外，识别偏弱）
    ("Clinochlore_Fe",  "铁绿泥石",     "#014f86"),     # 绿泥石族：富铁，Fe-OH 2.25μm（窗口内，识别较好）
    ("Montmorillonite", "蒙脱石",       "#9b5de5"),     # 粘土矿物，Al-OH 2.20μm
]
MINERAL_EN = [item[0] for item in TARGET_MINERALS]              # 英文名列表（决定类别号）
MINERAL_CN = {item[0]: item[1] for item in TARGET_MINERALS}     # 英文名 → 中文名
MINERAL_COLOR = {item[0]: item[2] for item in TARGET_MINERALS}  # 英文名 → 绘图颜色

# 参数：光谱库筛选
LIB_MAX_SAMPLES_PER_MINERAL = 6     # 每类最多保留几个端元（按光谱仪优先级/样本名排序）
LIB_MIN_ENDMEMBERS = 3              # 建库后至少要有这么多端元，否则报错
# 光谱仪优先顺序（数字越小越优先）：ASD 系列 2151 通道 → BECK 480 通道 → AVIRIS 224
LIB_SPECTROMETER_RANK = {
    "ASDFRa": 0, "ASDFRb": 1, "ASDHRa": 2, "ASDHRb": 3,   # ASD 高分辨率，优先级最高
    "ASDNGa": 4, "ASDNGb": 5,                             # ASD 下一代
    "BECKa": 6, "BECKb": 7, "BECKc": 8, "BECKd": 9,       # BECK 次之
    "AVIRIS": 10,                                         # AVIRIS 最后
}
# 光谱仪族 → (波长文件, 期望通道数)
LIB_GRIDS = {
    "ASD": ("splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt", 2151),    # ASD 2151 通道
    "BECK": ("splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt", 480),    # BECK 480 通道
    "AVIRIS": ("splib07a_Wavelengths_AVIRIS_1996_0.37-2.5_microns.txt", 224),  # AVIRIS 224 通道
}
LIB_WAVELENGTH_FILES = {
    "ASD": "splib07a_Wavelengths_ASD_0.35-2.5_microns_2151_ch.txt",            # ASD 波长文件
    "BECK": "splib07a_Wavelengths_BECK_Beckman_0.2-3.0_microns.txt",           # BECK 波长文件
    "AVIRIS": "splib07a_Wavelengths_AVIRIS_1996_0.37-2.5_microns.txt",         # AVIRIS 波长文件
}

# 参数：库→AHSI 重采样
LIB_MAX_WAVELENGTH_UM = 2.50   # 库光谱波长上限；AHSI 末波段 2.5011 μm 不参与匹配
# 参与匹配（SAM/SFF/MF 的背景统计与丰度）的波长上限：本数据 2.35 μm 以后
# 辐亮度只有 1~6 DN、地表反射率大于表观 TOA 反射率（物理上不可能），
# 属于噪声被 6S 强水汽吸收放大的伪值，故匹配一律只用 ≤2.30 μm 的波段。
MATCH_MAX_WAVELENGTH_UM = 2.30   # 参与 SAM/SFF/MF 的波长上限

# 参数：Al-OH 2.2 μm 吸收特征
# 连续统锚点（左右肩）：2.08 / 2.30 μm；特征评价范围：2.12~2.28 μm。
# 说明：本数据 2.35 μm 以后辐射/反射率快速趋近 0（长波 SWIR 信噪比很低），
# 若把右肩放到 2.45 μm 会把连续统压到近零、吸收深度失真，故窗口收窄到 2.30。
ALOH_CONTINUUM_RANGE = (2.08, 2.30)   # 连续统左右肩波长（μm）
ALOH_FEATURE_RANGE = (2.12, 2.28)     # 吸收特征评价范围（μm）
MIN_ALOH_BANDS = 8             # 窗口内波段数下限，不足则报错

# 参数：像元有效性
REFL_MIN = -0.05               # 反射率下限（允许轻微负值，大气校正残留）
REFL_MAX = 1.20                # 反射率上限
NODATA_THRESHOLD = -9000.0     # 小于该值视为 nodata（-9999）
MIN_CONTINUUM_VALUE = 0.002    # 连续统绝对下限（反射率量级），低于此判为无效
USE_CLOUD_FILLED = False       # False: 只用 mask==1 清晰像元
NODATA_OUT = -9999.0           # 成果栅格统一 nodata 值
# 本数据 SWIR 个别波段（如 1.44 μm 附近）有约 40% 像元是步骤2 的 nodata，
# 这类“坏波段”整体剔除，而不是把一半像元丢掉。
MAX_BAND_MISSING_FRACTION = 0.02   # 清晰像元中缺测比例超过 2% 的波段不参与匹配
MAX_MISSING_BANDS_MATCH = 2        # 单像元参与匹配的波段中允许缺测的个数（缺测按 0 填充）

# 参数：分类阈值（先用默认值试点，按直方图/质检图再调）
MIN_ALOH_DEPTH = 0.005         # 2.2 μm 吸收深度下限（连续统去除后）
SAM_MAX_DEG = 10.0             # 连续统去除光谱角（CR-SAM）上限（度）
SFF_MAX_RMS = 0.30             # SFF 连续统去除拟合残差上限（本数据噪声大，先放宽）
CLASS_BY = "sam"               # 分类依据："sam"（默认）或 "sff"（按残差最小）

# 参数：MF / MTMF（辅助产物，不参与默认分类判决）
ENABLE_MF = True                            # 总开关：是否计算 MF/MTMF
BG_SAMPLE_PIXELS = 150000      # 背景统计抽样像元数
BG_RIDGE = 1.0e-6              # 协方差求逆的岭参数（相对迹）

# 参数：分块与压缩
BLOCK_ROWS = 96                             # 每块处理的行数
GDAL_CACHE_MB = 1024                        # GDAL 栅格缓存（MB）
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]        # 大文件自动切换 BigTIFF
GEOTIFF_COMPRESS = "DEFLATE"                # 压缩算法
GEOTIFF_COMPRESS_LEVEL = 6                  # 压缩级别

# 参数：质检图
MAKE_QC = True                              # 是否生成质检图
QC_MAX_SIZE = 1500                          # 质检图最长边像素上限
QC_DPI = 150                                # 输出分辨率
TRUE_COLOR_TARGETS = {"R": 0.66, "G": 0.55, "B": 0.47}   # 真彩色合成目标波长（μm）

# 文件名约定
IN_MOSAIC_NAME = "mosaic_reflectance.tif"   # 输入：镶嵌反射率
IN_MASK_NAME = "mosaic_mask.tif"            # 输入：掩膜
IN_WAVELENGTHS_NAME = "wavelengths_clean.npy"  # 输入：波长数组
IN_INFO_NAME = "mosaic_info.json"           # 输入：镶嵌元信息

LIB_DIR_NAME = "_library"                   # 匹配库缓存目录名
CATALOG_DIR_NAME = "_catalog"               # 批处理日志/汇总目录名
OUT_CLASS_NAME = "aloh_class.tif"           # 输出：分类图
OUT_CONFIDENCE_NAME = "aloh_confidence.tif"  # 输出：置信度
OUT_DEPTH_NAME = "aloh_depth.tif"           # 输出：吸收深度
OUT_POSITION_NAME = "aloh_position_um.tif"  # 输出：吸收中心波长
OUT_WIDTH_NAME = "aloh_width_nm.tif"        # 输出：半深宽
OUT_SAM_NAME = "sam_angle_deg.tif"          # 输出：CR-SAM 角度
OUT_SFF_ABUNDANCE_NAME = "sff_abundance.tif"  # 输出：SFF 丰度
OUT_SFF_RMS_NAME = "sff_rms.tif"            # 输出：SFF 残差
OUT_MF_NAME = "mf_abundance.tif"            # 输出：MF 丰度
OUT_MTMF_NAME = "mtmf_infeasibility.tif"    # 输出：MTMF 不可行性
OUT_LEGEND_NAME = "classification_legend.csv"  # 输出：图例
OUT_SUMMARY_NAME = "aloh_summary.csv"       # 输出：分类统计
OUT_INFO_NAME = "aloh_info.json"            # 输出：运行参数与统计
OUT_QC_NAME = "preview_quality.png"         # 输出：质检图
CLASS_NODATA = 255                          # 分类图 nodata（0=未分类）

LOG_NAME = "zy1e.step4"                     # 日志器名称


# =================================================================
# 通用工具
# =================================================================

def parse_args(argv=None):
    """解析命令行参数；默认值全部来自上方配置区。"""
    parser = argparse.ArgumentParser(                      # 创建命令行解析器
        description="ZY1E AHSI 步骤4：云母-绿泥石等蚀变矿物提取（USGS splib07a + SAM/SFF/MF）")
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP3,           # 输入根目录
                        help="步骤3 成果根目录（默认: %(default)s）")
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP4,        # 输出根目录
                        help="本步成果根目录（默认: %(default)s）")
    parser.add_argument("--lib", dest="lib_root", default=LIB_ROOT,              # 光谱库目录
                        help="USGS splib07a ASCII 光谱库目录（默认: %(default)s）")
    parser.add_argument("--area", action="append", default=None,                 # 按区名筛选
                        help="只处理匹配的研究区，可多次指定，如 \"2.LSK\" 或 \"*LSK*\"")
    parser.add_argument("--limit", type=int, default=0,                          # 限制区域个数
                        help="只处理前 N 个研究区（默认 0 = 全部）")
    parser.add_argument("--overwrite", action="store_true", default=False,       # 覆盖已有成果
                        help="已有成果也重做（默认跳过）")
    parser.add_argument("--rebuild-library", action="store_true", default=False,  # 强制重建匹配库
                        help="强制重建 AHSI 匹配库")
    parser.add_argument("--library-only", action="store_true", default=False,    # 只建库
                        help="只构建/检查 AHSI 匹配库，不做区域提取")
    parser.add_argument("--dry-run", action="store_true",                        # 只打印计划
                        help="只打印处理计划，不执行")
    parser.add_argument("--no-mf", action="store_true", default=False,           # 关闭 MF
                        help="关闭 MF/MTMF 计算（更快，产物更少）")
    parser.add_argument("--max-angle", type=float, default=SAM_MAX_DEG,          # SAM 角上限
                        help="SAM 角上限（度，默认: %(default)s）")
    parser.add_argument("--min-depth", type=float, default=MIN_ALOH_DEPTH,       # 深度下限
                        help="2.2 μm 吸收深度下限（默认: %(default)s）")
    parser.add_argument("--max-rms", type=float, default=SFF_MAX_RMS,            # SFF 残差上限
                        help="SFF 拟合残差上限（默认: %(default)s）")
    parser.add_argument("--class-by", choices=("sam", "sff"), default=CLASS_BY,  # 分类依据
                        help="分类依据（默认: %(default)s）")
    parser.add_argument("--log", default=None,                                   # 日志路径
                        help="日志文件路径（默认 <out>/_catalog/step4_batch.log）")
    parser.add_argument("--verbose", action="store_true", help="输出调试日志")   # 调试级别
    return parser.parse_args(argv)                         # 解析并返回参数对象


def setup_logging(log_path, verbose=False):
    """初始化日志：控制台 INFO（--verbose 时 DEBUG）+ 文件 DEBUG，重复调用会重建。"""
    logger = logging.getLogger(LOG_NAME)                   # 取得具名日志器
    logger.setLevel(logging.DEBUG)                         # 总级别设为 DEBUG
    logger.handlers.clear()                                # 清空旧处理器
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",   # 统一日志格式
                            datefmt="%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(stream=sys.stdout)     # 控制台处理器
    console.setLevel(logging.DEBUG if verbose else logging.INFO)   # 按 --verbose 定级别
    console.setFormatter(fmt)                              # 套用统一格式
    logger.addHandler(console)                             # 挂到日志器
    if log_path:                                           # 指定了日志文件才写文件
        try:
            parent = os.path.dirname(os.path.abspath(log_path))   # 日志文件所在目录
            if parent:
                os.makedirs(parent, exist_ok=True)         # 确保目录存在
            handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")  # 追加写
            handler.setLevel(logging.DEBUG)                # 文件记录最详细级别
            handler.setFormatter(fmt)                      # 套用统一格式
            logger.addHandler(handler)                     # 挂到日志器
        except OSError as exc:                             # 日志文件无法创建
            console.handle(logging.LogRecord(              # 退化为控制台告警
                LOG_NAME, logging.WARNING, __file__, 0,
                "无法写入日志文件 %s: %s", (log_path, exc), None))
    logger.propagate = False                               # 不再向根日志器传播
    return logger                                          # 返回配置好的日志器


def human_size(num_bytes):
    """把字节数转成 B/KB/MB/GB 的可读字符串（仅用于日志显示）。"""
    size = float(num_bytes)                 # 统一转浮点，避免整除
    for unit in ("B", "KB", "MB", "GB", "TB"):   # 逐级换算单位
        if size < 1024.0 or unit == "TB":   # 小于 1024 或已到最大单位
            return "%.1f %s" % (size, unit)  # 返回格式化结果
        size /= 1024.0                      # 继续除以 1024
    return "%.1f TB" % size                 # 兜底返回 TB


def read_text_lines(path):
    """读文本文件：去掉 UTF-8 BOM，统一成 \\n 换行后按行返回。"""
    with open(path, "r", encoding="utf-8", errors="ignore") as handle:  # 忽略非法编码
        text = handle.read().lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")   # 读全文并统一换行
    return text.split("\n")                 # 按行拆分返回


# =================================================================
# 传感器波段规格（波长 + FWHM）
# =================================================================

def load_sensor_fwhm(tar_dir, band_wavelengths, logger):
    """
    用 tar 内 ZY1E_AHSI_VNIR/SWIR_Spectralresponse.raw 给每个输出波段匹配 FWHM。
    返回 (fwhm_um, 来源说明)；找不到时用回退常量。
    """
    centers_nm = np.array([], np.float64)          # 响应曲线中心波长（nm）
    fwhms_nm = np.array([], np.float64)            # 响应曲线 FWHM（nm）
    source = "回退常量"                            # 默认来源标记
    try:
        tar_list = sorted(glob.glob(os.path.join(tar_dir, "*.tar")) +    # 收集 tar 包
                          glob.glob(os.path.join(tar_dir, "*.tar.gz")) +
                          glob.glob(os.path.join(tar_dir, "*.tgz")))
        if tar_list:                               # 有 tar 包才尝试读取
            with tarfile.open(tar_list[0], "r:*") as tar:   # 打开第一个包
                names = tar.getnames()             # 包内成员名列表
                for suffix in RESPONSE_SUFFIXES:   # 依次找 VNIR/SWIR 响应
                    hits = [name for name in names if name.lower().endswith(suffix.lower())]   # 匹配响应文件名
                    if not hits:                   # 没有该响应文件
                        continue
                    handle = tar.extractfile(hits[0])   # 取成员文件对象
                    if handle is None:             # 目录项等无效成员
                        continue
                    with handle:
                        raw = handle.read()        # 一次性读入原始字节
                    text = raw.decode("utf-8", errors="ignore")   # 按文本解码
                    for line in text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n").split("\n"):   # 逐行解析
                        line = line.strip()        # 去掉首尾空白
                        if not line:               # 跳过空行
                            continue
                        parts = re.split(r"[,\s]+", line)   # 按逗号/空白拆分
                        if len(parts) < 2:         # 至少要有两列
                            continue
                        try:
                            centers_nm = np.append(centers_nm, float(parts[0]))   # 中心波长
                            fwhms_nm = np.append(fwhms_nm, float(parts[1]))       # FWHM
                        except ValueError:         # 表头等非数值行
                            continue
            if centers_nm.size:                    # 成功读到响应曲线
                source = "%s（响应曲线 %d 个中心）" % (
                    os.path.basename(tar_list[0]), centers_nm.size)   # 记录来源
                logger.info("传感器响应: %s, FWHM %.2f~%.2f nm",
                            source, float(fwhms_nm.min()), float(fwhms_nm.max()))
    except Exception as exc:                       # 读取失败不致命
        logger.warning("读取光谱响应失败，使用回退 FWHM: %s", exc, exc_info=True)

    band_wl_nm = np.asarray(band_wavelengths, np.float64) * 1000.0   # 波段中心转 nm
    if centers_nm.size:                            # 有响应曲线则逐波段匹配
        idx = np.abs(band_wl_nm[:, None] - centers_nm[None, :]).argmin(axis=1)  # 最近中心
        dist = np.abs(band_wl_nm - centers_nm[idx])    # 匹配距离
        fwhm_nm = fwhms_nm[idx].copy()             # 取对应 FWHM
        bad = dist > 0.5  # 0.5 nm 以上认为没匹配上（理论上应精确命中）
        if np.any(bad):                            # 存在未命中波段
            fwhm_nm[bad] = SWIR_FWHM_NM_FALLBACK  # 未命中波段用回退值
            logger.warning("%d 个波段未在响应曲线中找到对应中心，按 SWIR FWHM 回退",
                           int(bad.sum()))
    else:
        fwhm_nm = np.where(band_wl_nm < VNIR_SWIR_SPLIT_UM * 1000.0,   # 按 1.05 μm 分界
                           VNIR_FWHM_NM_FALLBACK, SWIR_FWHM_NM_FALLBACK)
        logger.info("传感器响应: %s（VNIR %.2f nm / SWIR %.2f nm）",
                    source, VNIR_FWHM_NM_FALLBACK, SWIR_FWHM_NM_FALLBACK)
    return fwhm_nm / 1000.0, source                # 返回 μm 与来源说明


# =================================================================
# 光谱库：扫描 / 解析 / 传感器重采样
# =================================================================

def parse_record_header(line):
    """解析 splib07a 记录头，返回 (样本名, 光谱仪, 类型) 或 None。"""
    line = line.strip().lstrip("\ufeff")        # 去空白与 BOM
    match = re.match(r"^splib07a\s+Record=\d+:\s*(.*)$", line, re.IGNORECASE)  # 匹配记录头
    if not match:                               # 不是记录头
        return None
    rest = match.group(1).strip()               # 记录头后面的内容
    parts = rest.rsplit(None, 2)                # 从右切分成 3 段
    if len(parts) != 3:                         # 段数不对则视为无效
        return None
    sample, spectrometer, rtype = parts         # 样本名 / 光谱仪 / 类型
    return sample.strip(), spectrometer.strip(), rtype.strip().upper()   # 类型统一大写


def spectrometer_family(spectrometer):
    """把 splib07a 光谱仪代号归类为 ASD / BECK / AVIRIS（对应三套波长网格）。"""
    upper = spectrometer.upper()                # 统一转大写比较
    if upper.startswith("ASDFR") or upper.startswith("ASDHR") or upper.startswith("ASDNG"):
        return "ASD"                            # ASD 系列
    if upper.startswith("BECK"):
        return "BECK"                           # BECK 系列
    if upper.startswith("AVIRIS"):
        return "AVIRIS"                         # AVIRIS 系列
    return None                                 # 未知光谱仪


def load_lib_wavelength_grid(lib_root, family, logger, cache):
    """读某光谱仪族的波长网格（μm），校验通道数并用 cache 避免重复读盘。"""
    if family in cache:                             # 已读过则直接复用
        return cache[family]
    filename = LIB_WAVELENGTH_FILES[family]         # 对应族的波长文件名
    path = os.path.join(lib_root, filename)         # 拼出完整路径
    if not os.path.isfile(path):                    # 文件必须存在
        raise FileNotFoundError("缺少光谱库波长文件: %s" % path)
    values = []                                     # 波长值列表
    for line in read_text_lines(path)[1:]:          # 跳过首行表头
        line = line.strip()                         # 去首尾空白
        if not line:                                # 跳过空行
            continue
        try:
            values.append(float(line.split()[0]))   # 取第一列波长
        except ValueError:                          # 非数值行跳过
            continue
    grid = np.asarray(values, np.float64)           # 转成 float64 数组
    expected = LIB_GRIDS[family][1]                 # 该族期望通道数
    if grid.size != expected:                       # 通道数必须一致
        raise ValueError("光谱库波长文件 %s 有 %d 个通道，预期 %d"
                         % (os.path.basename(path), grid.size, expected))
    logger.info("库波长网格 %s: %d 通道, %.4f ~ %.4f μm",
                family, grid.size, grid.min(), grid.max())
    cache[family] = grid                            # 写入缓存
    return grid                                     # 返回波长网格


def scan_library(lib_root, logger):
    """扫描矿物章，返回 Al-OH 目标端元候选记录列表。"""
    chapter = os.path.join(lib_root, LIB_CHAPTER_DIR)    # 矿物章目录
    if not os.path.isdir(chapter):                       # 目录必须存在
        raise NotADirectoryError("光谱库矿物章目录不存在: %s" % chapter)
    paths = sorted(glob.glob(os.path.join(chapter, "*.txt")))   # 全部光谱文件
    logger.info("光谱库矿物章文件数: %d", len(paths))

    skipped_type = skipped_family = parse_fail = 0       # 三类跳过计数
    candidates = []                                      # 候选端元记录
    for path in paths:                                   # 逐个文件扫描
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as handle:
                head = handle.readline()                 # 只读首行记录头
        except OSError:                                  # 文件打不开
            parse_fail += 1
            continue
        parsed = parse_record_header(head)               # 解析记录头
        if parsed is None:                               # 解析失败
            parse_fail += 1
            continue
        sample, spectrometer, rtype = parsed             # 拆出三个字段
        if rtype != "AREF":                              # 只要原始反射率
            skipped_type += 1
            continue
        token = sample.split()[0]                        # 样本名首词即矿物名
        if token not in MINERAL_CN:                      # 不在 Al-OH 清单内
            continue
        family = spectrometer_family(spectrometer)       # 归类光谱仪
        if family is None:                               # 不支持的仪器
            skipped_family += 1
            continue
        rank = LIB_SPECTROMETER_RANK.get(spectrometer, 99)   # 优先级，未知排最后
        candidates.append({
            "path": path,                                # 源文件路径
            "file": os.path.basename(path),              # 源文件名
            "sample": sample,                            # 样本名
            "mineral": token,                            # 矿物名（英文）
            "spectrometer": spectrometer,                # 光谱仪代号
            "family": family,                            # 光谱仪族
            "rank": rank,                                # 优先级
        })
    logger.info("Al-OH 候选端元(AREF): %d；跳过 非AREF %d / 不支持光谱仪 %d / 头部解析失败 %d",
                len(candidates), skipped_type, skipped_family, parse_fail)

    by_mineral = {}                                      # 按矿物分组
    for record in candidates:
        by_mineral.setdefault(record["mineral"], []).append(record)
    selected = []                                        # 最终选中的端元
    for mineral in MINERAL_EN:                           # 按清单顺序选端元
        group = by_mineral.get(mineral, [])              # 该矿物的候选
        group.sort(key=lambda item: (item["rank"], item["sample"], item["file"]))   # 排序
        chosen = group[:LIB_MAX_SAMPLES_PER_MINERAL]     # 每类最多取 N 个
        if not chosen:                                   # 一个都没有
            logger.warning("库中缺少矿物 %s 的 AREF 端元", mineral)
        else:
            selected.extend(chosen)                      # 追加到结果
        logger.info("  %-16s 候选 %2d → 选用 %d（%s）", mineral, len(group),
                    len(chosen), ",".join(sorted({item["spectrometer"] for item in chosen})) or "-")
    if len(selected) < LIB_MIN_ENDMEMBERS:               # 端元太少无法建库
        raise ValueError("Al-OH 端元数 %d < %d，请检查光谱库路径/清单"
                         % (len(selected), LIB_MIN_ENDMEMBERS))
    return selected                                      # 返回选中记录


def load_spectrum_values(path):
    """
    读 splib07a ASCII 数值列（跳过首行记录头）。

    注意：splib07a ASCII 用 -1.23e34 作为“缺测”哨兵值（多见于 BECK 光谱仪的
    水汽带等位置），这里统一转成 NaN，交给后续卷积与插值处理。
    """
    values = []                                     # 光谱数值列表
    for line in read_text_lines(path)[1:]:          # 跳过首行记录头
        line = line.strip()                         # 去首尾空白
        if not line:                                # 跳过空行
            continue
        try:
            values.append(float(line.split()[0]))   # 取第一列数值
        except (ValueError, IndexError):            # 非数值/空列跳过
            continue
    array = np.asarray(values, np.float64)          # 转 float64 数组
    if array.size:                                  # 有数据才做清洗
        invalid = ~np.isfinite(array) | (np.abs(array) > 1.0e3)   # 哨兵值判定
        array[invalid] = np.nan                     # 哨兵/异常转 NaN
    return array                                    # 返回光谱数组


def _fill_spectrum_nan(spectrum):
    """按波长顺序对 NaN 做线性插值填补（只用于库光谱的少量缺测点）。"""
    spectrum = np.array(spectrum, np.float64, copy=True)   # 复制一份避免改原数组
    valid = np.isfinite(spectrum)                   # 有效值掩膜
    if valid.sum() < 5 or valid.all():              # 有效点太少或全有效则不处理
        return spectrum
    index = np.arange(spectrum.size)                # 波段序号
    spectrum[~valid] = np.interp(index[~valid], index[valid], spectrum[valid])   # 线性插值
    return spectrum                                 # 返回填补后的光谱


def gaussian_resample(wl_lib_um, spectra, band_wl_um, band_fwhm_um):
    """
    把库光谱（N, nlib）按 AHSI 波段中心做高斯卷积，返回 (N, nband) float64。

    逐波段独立采样：本数据 AHSI 波长非单调（1.0402 → 1.0056），不能整体插值。
    每个波段只取 ±3σ 内的库采样点并归一化权重；库光谱中的 NaN（缺测哨兵）
    不参与加权；窗口内没有有效采样点时返回 NaN。
    """
    spectra = np.asarray(spectra, np.float64)       # 统一转 float64
    n_band = band_wl_um.size                        # 输出波段数
    out = np.full((spectra.shape[0], n_band), np.nan, np.float64)   # 结果先填 NaN
    coverage = np.zeros(n_band, np.float64)         # 每波段采样覆盖率
    sigma = band_fwhm_um / 2.354820045              # FWHM → 高斯 σ
    spacing = float(np.median(np.abs(np.diff(wl_lib_um))))   # 库采样间隔中位数
    for index in range(n_band):                     # 逐输出波段卷积
        center = band_wl_um[index]                  # 波段中心波长
        width = sigma[index]                        # 波段 σ
        if not np.isfinite(center) or not np.isfinite(width) or width <= 0:
            continue                                # 无效波段跳过
        mask = np.abs(wl_lib_um - center) <= 3.0 * width   # ±3σ 窗口
        count = int(mask.sum())                     # 窗口内采样点数
        expected = max(1.0, (6.0 * width) / spacing)   # 理论应有点数
        coverage[index] = min(1.0, count / expected)   # 覆盖率
        if count < 1:                               # 窗口内没有采样点
            continue
        delta = (wl_lib_um[mask] - center) / width  # 归一化偏移
        weights = np.exp(-0.5 * delta * delta)      # 高斯权重
        subset = spectra[:, mask]                   # 取窗口内光谱
        finite = np.isfinite(subset)                # 有效点掩膜
        weight_sum = finite @ weights               # 有效权重和
        weighted = np.where(finite, subset, 0.0) @ weights   # 加权和
        good = weight_sum > 0                       # 有有效权重才算
        out[good, index] = weighted[good] / weight_sum[good]  # 归一化结果
    return out, coverage                            # 返回重采样光谱与覆盖率


def build_library(lib_root, band_wl, band_fwhm, out_library_dir, logger, force=False):
    """建立（或复用）AHSI 匹配库。

    输出 _library/：
        ahsi_library.npy       (端元数, 142) 传感器匹配后的库光谱
        ahsi_library_meta.csv  端元清单（类别 / 样本 / 光谱仪 / 源文件）
        ahsi_band_spec.csv     波段表（中心波长 / FWHM / 匹配与窗口标记）
    已存在且 force=False 时直接加载，避免每次运行重复扫描光谱库。
    """
    npy_path = os.path.join(out_library_dir, "ahsi_library.npy")      # 库光谱矩阵
    meta_path = os.path.join(out_library_dir, "ahsi_library_meta.csv")  # 端元清单
    band_path = os.path.join(out_library_dir, "ahsi_band_spec.csv")   # 波段表
    if not force and os.path.isfile(npy_path) and os.path.isfile(meta_path):   # 可复用
        library = np.load(npy_path)                 # 读库光谱
        meta = []                                   # 端元元数据
        with open(meta_path, "r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle):      # 逐行读 CSV
                meta.append(row)                    # 收集端元行
        if library.ndim != 2 or library.shape[1] != band_wl.size or len(meta) != library.shape[0]:   # 校验形状一致
            logger.warning("已有匹配库与当前波段不一致，强制重建")
        else:
            logger.info("复用已有匹配库: %s (%d 端元 × %d 波段)",
                        npy_path, library.shape[0], library.shape[1])
            return _assemble_library(library, meta, band_wl, band_fwhm)   # 直接组装返回

    logger.info("构建 AHSI 匹配库 ...")
    os.makedirs(out_library_dir, exist_ok=True)     # 确保库目录存在
    records = scan_library(lib_root, logger)        # 扫描候选端元

    cache = {}                                      # 波长网格缓存
    endmembers, meta_rows, coverage_rows = [], [], []   # 端元/元数据/覆盖率
    skipped_mismatch = 0                            # 通道数不一致计数
    filled_values = 0                               # 缺测填补计数
    for record in records:                          # 逐个端元重采样
        grid = load_lib_wavelength_grid(lib_root, record["family"], logger, cache)  # 库波长
        values = load_spectrum_values(record["path"])   # 读库光谱
        if values.size != grid.size:                # 通道数必须一致
            skipped_mismatch += 1
            logger.warning("通道数不一致，跳过: %s (%d vs %d)",
                           record["file"], values.size, grid.size)
            continue
        resampled, coverage = gaussian_resample(grid, values[None, :], band_wl, band_fwhm)  # 重采样
        row = resampled[0]                          # 取出该端元结果
        if not np.all(np.isfinite(row)):            # 有缺测则插值填补
            filled_values += int(np.count_nonzero(~np.isfinite(row)))   # 统计缺测个数
            row = _fill_spectrum_nan(row)           # 线性插值填补
        if not np.all(np.isfinite(row)):            # 仍无效则跳过
            skipped_mismatch += 1
            logger.warning("缺测值过多，跳过: %s", record["file"])
            continue
        endmembers.append(row)                      # 收集端元光谱
        coverage_rows.append(coverage)              # 收集覆盖率
        meta_rows.append({
            "endmember_id": len(meta_rows),         # 端元序号
            "class_id": MINERAL_EN.index(record["mineral"]) + 1,   # 类别号从 1 起
            "mineral_en": record["mineral"],        # 矿物英文名
            "mineral_cn": MINERAL_CN[record["mineral"]],   # 矿物中文名
            "sample": record["sample"],             # 样本名
            "spectrometer": record["spectrometer"],  # 光谱仪
            "source_file": record["file"],          # 源文件名
        })
    if skipped_mismatch:                            # 有跳过则汇总提示
        logger.warning("共 %d 个端元因通道数不一致被跳过", skipped_mismatch)
    if filled_values:                               # 有填补则记录一次
        logger.info("库光谱缺测值填补: %d 个（splib07a -1.23e34 哨兵/检出限空值）",
                    filled_values)
    if len(endmembers) < LIB_MIN_ENDMEMBERS:        # 端元数不足
        raise ValueError("有效端元数 %d < %d，无法建库" % (len(endmembers), LIB_MIN_ENDMEMBERS))

    library = np.vstack(endmembers).astype(np.float32)   # 堆叠成矩阵
    np.save(npy_path, library)                      # 保存库光谱
    write_meta_csv(meta_path, meta_rows)            # 写端元清单
    coverage_mean = np.mean(np.vstack(coverage_rows), axis=0)   # 求平均覆盖率

    with open(band_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)                 # CSV 写入器
        writer.writerow(["band", "wavelength_um", "fwhm_nm", "used_in_matching",   # 表头
                         "in_aloh_window", "lib_coverage"])
        for index, (wl, fwhm) in enumerate(zip(band_wl, band_fwhm), start=1):   # 逐波段写
            writer.writerow([index, "%.6f" % wl, "%.3f" % (fwhm * 1000.0),
                             "TRUE" if wl <= MATCH_MAX_WAVELENGTH_UM else "FALSE",   # 是否参与匹配
                             "TRUE" if ALOH_CONTINUUM_RANGE[0] <= wl <= ALOH_CONTINUUM_RANGE[1] else "FALSE",   # 是否在窗口
                             "%.3f" % coverage_mean[index - 1]])

    logger.info("匹配库: %d 端元 × %d 波段，%s", library.shape[0], library.shape[1],
                human_size(os.path.getsize(npy_path)))
    make_library_overview(out_library_dir, band_wl, library, meta_rows, logger)   # 画库光谱图
    return _assemble_library(library, meta_rows, band_wl, band_fwhm)   # 组装返回


def _assemble_library(library, meta_rows, band_wl, band_fwhm):
    """把光谱矩阵与元数据组装成后续算法使用的 dict（含按矿物类分组的端元索引）。"""
    class_ids = np.asarray([int(row["class_id"]) for row in meta_rows], np.int32)  # 端元类别号
    groups = []                                     # 类别 → 端元下标
    for class_index in range(1, len(MINERAL_EN) + 1):   # 按类别号遍历
        groups.append((class_index, np.where(class_ids == class_index)[0]))
    return {
        "spectra": np.asarray(library, np.float32),   # 库光谱矩阵
        "meta": meta_rows,                          # 端元元数据
        "class_ids": class_ids,                     # 端元类别号
        "groups": groups,                           # 分组索引
        "band_wl": np.asarray(band_wl, np.float64),   # 波段中心（μm）
        "band_fwhm": np.asarray(band_fwhm, np.float64),   # 波段 FWHM（μm）
        "match_idx": np.where(np.asarray(band_wl) <= MATCH_MAX_WAVELENGTH_UM)[0],   # 参与匹配的波段
    }


def write_meta_csv(path, meta_rows):
    """写出端元清单 CSV（可直接在 Excel 中核对库的组成）。"""
    fields = ["endmember_id", "class_id", "mineral_en", "mineral_cn",
              "sample", "spectrometer", "source_file"]   # CSV 列顺序
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)   # 字典写入器
        writer.writeheader()                        # 写表头
        for row in meta_rows:                       # 逐端元写一行
            writer.writerow(row)


def make_library_overview(out_dir, band_wl, library, meta_rows, logger):
    """画库质检图：左图为全谱，右图为 2.0~2.5 μm 放大（看 Al-OH 吸收）。"""
    try:
        import matplotlib                         # 绘图库
        matplotlib.use("Agg")                     # 无界面后端
        import matplotlib.pyplot as plt
    except Exception as exc:                      # 没装 matplotlib 则跳过
        logger.warning("matplotlib 不可用，跳过库光谱图: %s", exc)
        return
    class_ids = np.asarray([int(row["class_id"]) for row in meta_rows])   # 端元类别号
    class_names = [item[0] for item in TARGET_MINERALS]   # 类别英文名
    figure, axes = plt.subplots(1, 2, figsize=(15, 5.6))   # 左全谱右放大
    for class_index, class_name in enumerate(class_names, start=1):   # 逐类绘制
        idx = np.where(class_ids == class_index)[0]   # 该类端元下标
        if idx.size == 0:                         # 该类没有端元
            continue
        color = MINERAL_COLOR[class_name]         # 该类绘图颜色
        for row in idx:                           # 逐个端元画曲线
            axes[0].plot(band_wl, library[row], color=color, linewidth=0.7, alpha=0.6,
                         label=class_name if row == idx[0] else None)   # 只标一次图例
            zoom = (band_wl >= 2.0) & (band_wl <= 2.5)   # 2.0~2.5 μm 放大范围
            axes[1].plot(band_wl[zoom], library[row, zoom], color=color, linewidth=0.8,
                         alpha=0.7)
    axes[0].set_xlabel("Wavelength (um)")         # 横轴：波长
    axes[0].set_ylabel("Reflectance")             # 纵轴：反射率
    axes[0].set_title("AHSI-matched splib07a endmembers (n=%d)" % library.shape[0])   # 左图标题
    axes[0].grid(alpha=0.3)                       # 网格线
    axes[1].axvspan(2.10, 2.40, color="0.9", zorder=0)   # 2.2 μm 特征区底纹
    axes[1].set_xlim(2.0, 2.5)                    # 放大图横轴范围
    axes[1].set_xlabel("Wavelength (um)")         # 放大图横轴
    axes[1].set_title("2.2 um feature region")   # 放大图标题
    axes[1].grid(alpha=0.3)                       # 放大图网格
    handles, labels = axes[0].get_legend_handles_labels()   # 收集图例
    axes[0].legend(handles, labels, fontsize=7, ncol=2, loc="best")   # 绘制图例
    figure.tight_layout()                         # 紧凑布局
    path = os.path.join(out_dir, "library_overview.png")   # 输出路径
    figure.savefig(path, dpi=130, bbox_inches="tight")     # 保存图片
    plt.close(figure)                             # 释放画布
    logger.info("库光谱图: %s", path)


# =================================================================
# 区域处理
# =================================================================

def _clip_suffix(filename):
    for suffix in STEP3_CLIP_SUFFIXES:               # 依次比对已知后缀
        if filename.endswith(suffix):                # 命中后缀
            return suffix
    return None                                      # 不是裁剪成果


def area_clip_files(in_root, area):
    """列出 <in_root>/<area>/ 下的逐景裁剪反射率（按文件名排序）。"""
    area_dir = os.path.join(in_root, area)           # 区域目录
    if not os.path.isdir(area_dir):                  # 目录不存在
        return []
    return [os.path.join(area_dir, name) for name in sorted(os.listdir(area_dir))
            if _clip_suffix(name) and os.path.isfile(os.path.join(area_dir, name))]   # 只留裁剪成果


def _scene_wavelengths_path(scene):
    for root in (STEP1_ROOT, STEP2_ROOT):            # 先步骤1 再步骤2
        path = os.path.join(root, scene, "wavelengths_clean.npy")   # 场景波长文件
        if os.path.isfile(path):                     # 找到即返回
            return path
    return None                                      # 两处都没有


def load_area_wavelengths(clips, band_count=None):
    """由裁剪文件名反推场景名，读取步骤1/2 的 wavelengths_clean.npy。"""
    for path in clips:                               # 逐个裁剪文件尝试
        name = os.path.basename(path)                # 文件名
        suffix = _clip_suffix(name)                  # 裁剪后缀
        scene = name[:-len(suffix)] if suffix else os.path.splitext(name)[0]   # 反推场景名
        wl_path = _scene_wavelengths_path(scene)     # 找波长文件
        if wl_path is None:                          # 找不到则试下一景
            continue
        wavelengths = np.load(wl_path).astype(np.float64)   # 读波长数组
        if wavelengths.ndim != 1:                    # 必须是一维
            raise ValueError("波长文件不是 1 维: %s" % wl_path)
        if band_count is not None and wavelengths.size != band_count:   # 波段数必须一致
            raise ValueError("裁剪影像波段数 %d != 波长数 %d（%s）",
                             band_count, wavelengths.size, wl_path)
        return wavelengths, scene                   # 返回波长与场景名
    raise FileNotFoundError(
        "未找到裁剪场景对应的 wavelengths_clean.npy（查过步骤1/2 目录）")


def build_area_vrt(clips, vrt_path):
    """用逐景裁剪成果构建虚拟镶嵌（不落地大栅格）。"""
    tmp_path = vrt_path + ".tmp.vrt"                 # 先写临时 VRT
    options = gdal.BuildVRTOptions(resolution="highest", VRTNodata=NODATA_OUT)   # 最高分辨率
    vrt = gdal.BuildVRT(tmp_path, clips, options=options)   # 构建虚拟镶嵌
    if vrt is None:                                  # 构建失败
        raise RuntimeError("构建 VRT 失败: %s" % vrt_path)
    vrt.FlushCache()                                 # 落盘
    vrt = None                                       # 释放句柄
    os.replace(tmp_path, vrt_path)                   # 原子改名


def build_area_mask(vrt_path, mask_path, logger):
    """按抽样波段的有效值与反射率范围生成 0/1 掩膜（1=清晰像元）。"""
    ds = gdal.Open(vrt_path)                         # 打开镶嵌
    if ds is None:                                   # 打开失败
        raise RuntimeError("无法打开镶嵌: %s" % vrt_path)
    width, height, nbands = ds.RasterXSize, ds.RasterYSize, ds.RasterCount   # 尺寸与波段
    count = min(MASK_SAMPLE_BANDS, nbands)           # 抽样波段数
    sample = np.unique(np.linspace(0, nbands - 1, count).astype(int))   # 均匀抽样的波段下标
    tmp_path = mask_path + ".tmp.tif"                # 掩膜临时文件
    driver = gdal.GetDriverByName("GTiff")           # GeoTIFF 驱动
    out = driver.Create(tmp_path, width, height, 1, gdal.GDT_Byte,
                        options=["TILED=YES", "COMPRESS=DEFLATE", "ZLEVEL=6"])   # 单波段 Byte
    if out is None:                                  # 创建失败
        raise RuntimeError("无法创建掩膜: %s" % mask_path)
    out.SetGeoTransform(ds.GetGeoTransform())        # 复制地理变换
    out.SetProjection(ds.GetProjection())            # 复制投影
    band = out.GetRasterBand(1)                      # 掩膜波段
    band.SetNoDataValue(0)                           # 0 为无效值
    clear_total = 0                                  # 清晰像元累计
    for y in range(0, height, BLOCK_ROWS):           # 分块处理
        rows = min(BLOCK_ROWS, height - y)           # 本块行数
        valid = np.ones((rows, width), bool)         # 有效性掩膜
        low = np.full((rows, width), np.inf)         # 各像元最小值
        high = np.full((rows, width), -np.inf)       # 各像元最大值
        for index in sample:                         # 逐个抽样波段
            array = ds.GetRasterBand(int(index) + 1).ReadAsArray(0, y, width, rows)   # 读块
            array = np.asarray(array, np.float64)    # 转浮点
            ok = _valid_pixels(array)                # 有效值判定
            valid &= ok                              # 累积有效
            low = np.where(ok, np.minimum(low, array), low)     # 更新最小值
            high = np.where(ok, np.maximum(high, array), high)  # 更新最大值
        valid &= (low > REFL_MIN) & (high < REFL_MAX)   # 反射率需在合理范围
        mask = valid.astype(np.uint8)                # 转 0/1
        band.WriteArray(mask, 0, y)                  # 写块
        clear_total += int(mask.sum())               # 累计清晰像元
    out.FlushCache()                                 # 落盘
    out = None                                       # 释放句柄
    ds = None                                        # 关闭输入数据集
    os.replace(tmp_path, mask_path)                  # 原子改名
    logger.info("生成掩膜: %s（清晰像元 %.1f%%）",
                mask_path, 100.0 * clear_total / max(1, width * height))


def resolve_area_inputs(in_root, area, cache_root, overwrite, logger, build=True):
    """解析区域输入：优先 mosaic 三件套，否则由步骤3 裁剪成果构建 VRT + 掩膜。"""
    mosaic_dir = os.path.join(in_root, area, "mosaic")   # mosaic 子目录
    mosaic_tif = os.path.join(mosaic_dir, IN_MOSAIC_NAME)   # 镶嵌反射率
    mosaic_mask = os.path.join(mosaic_dir, IN_MASK_NAME)    # 镶嵌掩膜
    mosaic_wl = os.path.join(mosaic_dir, IN_WAVELENGTHS_NAME)   # 镶嵌波长
    if (os.path.isfile(mosaic_tif) and os.path.isfile(mosaic_mask)
            and os.path.isfile(mosaic_wl)):          # 三件套齐全则直接用
        return {"refl_path": mosaic_tif, "mask_path": mosaic_mask,
                "wavelengths": np.load(mosaic_wl).astype(np.float64),
                "source": "mosaic", "scenes": []}

    clips = area_clip_files(in_root, area)           # 否则列出裁剪成果
    if not clips:                                    # 两种输入都没有
        raise FileNotFoundError(
            "区域 %s 下既没有 mosaic/ 三件套，也没有步骤3 裁剪成果" % area)
    probe = gdal.Open(clips[0])                      # 打开首个文件探波段数
    if probe is None:                                # 打开失败
        raise RuntimeError("无法打开裁剪成果: %s" % clips[0])
    band_count = probe.RasterCount                   # 波段数
    probe = None                                     # 释放句柄
    wavelengths, _ = load_area_wavelengths(clips, band_count)   # 读波长

    cache_dir = os.path.join(cache_root, CACHE_DIR_NAME, area)   # 区域缓存目录
    vrt_path = os.path.join(cache_dir, CACHE_VRT_NAME)   # VRT 路径
    mask_path = os.path.join(cache_dir, IN_MASK_NAME)    # 掩膜路径
    if not build:                                    # 只解析不构建（dry-run）
        return {"refl_path": vrt_path, "mask_path": mask_path,
                "wavelengths": wavelengths, "source": "clip", "scenes": clips}

    os.makedirs(cache_dir, exist_ok=True)            # 确保缓存目录存在
    if overwrite or not os.path.isfile(vrt_path):    # 需要时构建 VRT
        build_area_vrt(clips, vrt_path)
        logger.info("[%s] 合并 %d 景裁剪成果 -> %s", area, len(clips), vrt_path)
    if overwrite or not os.path.isfile(mask_path):   # 需要时构建掩膜
        build_area_mask(vrt_path, mask_path, logger)
    return {"refl_path": vrt_path, "mask_path": mask_path,
            "wavelengths": wavelengths, "source": "clip", "scenes": clips}


def discover_areas(in_root, patterns=None):
    areas = []                                       # 可处理区域名
    if not os.path.isdir(in_root):                   # 根目录必须存在
        raise NotADirectoryError("步骤3 成果目录不存在: %s" % in_root)
    for name in sorted(os.listdir(in_root)):         # 逐个子目录检查
        if not os.path.isdir(os.path.join(in_root, name)):   # 只处理目录
            continue
        has_mosaic = os.path.isfile(
            os.path.join(in_root, name, "mosaic", IN_MOSAIC_NAME))   # 是否有镶嵌
        if not has_mosaic and not area_clip_files(in_root, name):    # 两种输入都没有
            continue
        if patterns and not any(fnmatch.fnmatch(name, pattern) for pattern in patterns):
            continue                                 # 不符合 --area 筛选
        areas.append(name)                           # 记录该区域
    return areas                                     # 返回区域列表


def _nearest_band(wavelengths, target):
    """返回中心波长最接近 target（μm）的波段下标。"""
    return int(np.argmin(np.abs(np.asarray(wavelengths) - target)))   # 最近波段下标


def _stretch(array, valid):
    """对影像做 2%~98% 线性拉伸到 0~1，仅用于质检出图。"""
    values = array[valid]                            # 只取有效值
    if values.size == 0:                             # 没有有效值
        return np.zeros_like(array)
    low, high = np.percentile(values, [2, 98])       # 2%/98% 分位
    if high <= low:                                  # 防止除零
        high = low + 1e-6
    return np.clip((array - low) / (high - low), 0.0, 1.0)   # 拉伸并截断


def _valid_pixels(array):
    """有效值：有限且不是 -9999 nodata。"""
    return np.isfinite(array) & (array > NODATA_THRESHOLD)   # 有限且高于 nodata 阈值


def estimate_band_missing(ds, mask_ds, candidate_idx, logger, max_rows=200):
    """抽稀行统计“清晰像元中各候选波段的有效比例”，用于剔除坏波段。"""
    width, height = ds.RasterXSize, ds.RasterYSize   # 影像尺寸
    step = max(1, height // max_rows)                # 抽样行间隔
    valid_counts = np.zeros(candidate_idx.size, np.float64)   # 各波段有效像元数
    clear_total = 0                                  # 清晰像元总数
    for y in range(step // 2, height, step):         # 等间隔抽行
        cube = ds.ReadAsArray(0, y, width, 1)        # 读一行全部波段
        if cube is None:                             # 读取失败
            continue
        mask = mask_ds.ReadAsArray(0, y, width, 1)   # 读对应掩膜行
        if mask.ndim == 2:                           # 去掉波段维
            mask = mask[0]
        clear = mask == 1                            # 清晰像元
        if not np.any(clear):                        # 本行没有清晰像元
            continue
        block = cube[candidate_idx][:, 0, :]         # 取候选波段
        valid = _valid_pixels(block)                 # 有效值掩膜
        valid_counts += valid[:, clear].sum(axis=1)  # 累计有效数
        clear_total += int(np.count_nonzero(clear))  # 累计清晰像元
    if clear_total == 0:                             # 一个清晰像元都没有
        raise RuntimeError("采样行中没有清晰像元，检查掩膜")
    valid_fraction = valid_counts / float(clear_total)   # 各波段有效比例
    bad = valid_fraction < (1.0 - MAX_BAND_MISSING_FRACTION)   # 坏波段判定
    logger.info("坏波段检查: 采样 %d 个清晰像元；有效比例 <%.0f%% 的波段 %d 个",
                clear_total, 100.0 * (1.0 - MAX_BAND_MISSING_FRACTION),
                int(np.count_nonzero(bad)))
    return valid_fraction, clear_total               # 返回有效比例与样本数


def sample_background(ds, mask_ds, match_idx, target, logger):
    """从清晰像元抽稀采样，估计背景均值与协方差逆（n_match×n_match）。"""
    width, height = ds.RasterXSize, ds.RasterYSize   # 影像尺寸
    step = max(1, height // 200)                     # 抽样行间隔
    chunks = []                                      # 分块样本
    total = 0                                        # 已采样像元数
    for y in range(step // 2, height, step):         # 等间隔抽行
        cube = ds.ReadAsArray(0, y, width, 1)        # 读一行全部波段
        if cube is None:                             # 读取失败
            continue
        mask = mask_ds.ReadAsArray(0, y, width, 1)   # 对应掩膜行
        if mask.ndim == 2:                           # 去掉波段维
            mask = mask[0]
        block = cube[match_idx].reshape(len(match_idx), -1).astype(np.float64)   # 匹配波段
        finite = _valid_pixels(block)                # 有效值掩膜
        missing = np.count_nonzero(~finite, axis=0)  # 每个像元缺测波段数
        low = np.where(finite, block, np.inf).min(axis=0)    # 每像元最小值
        high = np.where(finite, block, -np.inf).max(axis=0)  # 每像元最大值
        ok = ((mask.reshape(-1) == 1) & (missing <= MAX_MISSING_BANDS_MATCH) &
              (low > REFL_MIN) & (high < REFL_MAX))  # 清晰且反射率合理
        if not np.any(ok):                           # 本行无可用像元
            continue
        chunks.append(np.where(finite[:, ok], block[:, ok], 0.0).T)   # 缺测补 0
        total += int(ok.sum())                       # 累计样本数
        if total >= target:                          # 达到目标即停止
            break
    if not chunks:                                   # 没采到任何样本
        raise RuntimeError("无法获取清晰像元样本，检查掩膜/反射率范围")
    samples = np.vstack(chunks)                      # 合并所有样本
    if samples.shape[0] > target:                    # 超量则等间隔抽稀
        stride = max(1, samples.shape[0] // target)
        samples = samples[::stride][:target]
    mean = samples.mean(axis=0)                      # 背景均值
    covariance = np.cov(samples, rowvar=False)       # 背景协方差
    ridge = BG_RIDGE * np.trace(covariance) / covariance.shape[0]   # 岭参数
    covariance = covariance + np.eye(covariance.shape[0]) * ridge   # 加岭稳定求逆
    inverse = np.linalg.pinv(covariance, rcond=1e-10)   # 伪逆
    logger.info("背景统计: %d 像元, 条件数 1e%.1f",
                samples.shape[0], np.log10(np.linalg.cond(covariance) + 1e-30))
    return mean, inverse                             # 返回均值与协方差逆


def _class_reduce(values, groups, reducer, fill):
    """values: (n_end, n) → (n_class, n)，按端元分组做 min/max。"""
    n_class = len(groups)                            # 类别数
    out = np.full((n_class, values.shape[1]), fill, np.float64)   # 先填初值
    for position, (_, idx) in enumerate(groups):     # 逐类聚合
        if idx.size:                                 # 该类有端元
            out[position] = reducer(values[idx], axis=0)   # 类内取 min/max
    return out                                       # 返回类别级结果


def _create_gtiff(tmp_path, width, height, nbands, dtype, geo, projection,
                  nodata, descriptions=None):
    """建临时 GeoTIFF（分块 + 压缩 + 波段描述），成功后由调用方原子改名。"""
    creation = list(GTIFF_OPTIONS)                   # 复制基础创建选项
    if GEOTIFF_COMPRESS:                             # 需要压缩时追加参数
        creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS,
                         "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL,
                         "INTERLEAVE=BAND"])
        if dtype != gdal.GDT_Byte:                   # 非整型加预测器
            creation.append("PREDICTOR=3")
    driver = gdal.GetDriverByName("GTiff")           # GeoTIFF 驱动
    dataset = driver.Create(tmp_path, width, height, nbands, dtype, options=creation)   # 建临时文件
    if dataset is None:                              # 创建失败
        raise RuntimeError("无法创建 GeoTIFF: %s" % tmp_path)
    dataset.SetGeoTransform(geo)                     # 写地理变换
    if projection:                                   # 有投影则写入
        dataset.SetProjection(projection)
    for index in range(1, nbands + 1):               # 逐波段设置
        band = dataset.GetRasterBand(index)          # 取波段句柄
        band.SetNoDataValue(nodata)                  # 设 nodata
        if descriptions and index <= len(descriptions):   # 有描述则设置
            band.SetDescription(str(descriptions[index - 1]))
    return dataset                                   # 返回数据集


def process_area(area_name, inputs, out_root, library, response_source, args, logger):
    started = time.time()                            # 记录开始时间
    in_tif = inputs["refl_path"]                     # 输入反射率路径
    in_mask = inputs["mask_path"]                    # 输入掩膜路径
    wavelengths = inputs["wavelengths"]              # 波长数组
    out_dir = os.path.join(out_root, area_name)      # 本区输出目录

    final_class = os.path.join(out_dir, OUT_CLASS_NAME)   # 分类图路径
    final_info = os.path.join(out_dir, OUT_INFO_NAME)     # 信息 JSON 路径
    if not args.overwrite and os.path.isfile(final_class) and os.path.isfile(final_info):
        logger.info("[%s] 跳过（成果已存在）: %s", area_name, out_dir)
        return {"area": area_name, "status": "跳过", "message": "成果已存在",
                "elapsed": 0.0, "classes": None,
                "valid_pixels": 0, "classified_pixels": 0}   # 直接返回跳过

    if wavelengths.ndim != 1:                        # 波长必须是一维
        raise ValueError("波长数组不是 1 维: %s" % in_tif)

    # ========== ① 输入校验：波段数 / 波长数 / 掩膜必须一致 ==========
    source_note = ("mosaic" if inputs["source"] == "mosaic"
                   else "裁剪成果合并 %d 景" % len(inputs["scenes"]))   # 输入来源说明
    logger.info("[%s] 开始处理：%s（%s）", area_name, in_tif, source_note)
    logger.info("[%s] 掩膜: %s；波段数 %d", area_name, in_mask, wavelengths.size)

    ds = gdal.Open(in_tif)                           # 打开镶嵌反射率
    if ds is None:                                   # 打开失败
        raise RuntimeError("无法打开镶嵌反射率: %s" % in_tif)
    if ds.RasterCount != wavelengths.size:           # 波段数必须与波长数一致
        raise ValueError("镶嵌波段数 %d != 波长数 %d" % (ds.RasterCount, wavelengths.size))
    width, height = ds.RasterXSize, ds.RasterYSize   # 影像宽高
    geo = ds.GetGeoTransform()                       # 地理变换
    projection = ds.GetProjection()                  # 投影 WKT
    pixel_area_m2 = abs(float(geo[1]) * float(geo[5]))   # 30 m 网格为 900，60 m 网格为 3600

    mask_ds = gdal.Open(in_mask)                     # 打开掩膜
    if mask_ds is None:                              # 打开失败
        raise RuntimeError("无法打开掩膜: %s" % in_mask)
    if (mask_ds.RasterXSize, mask_ds.RasterYSize) != (width, height):
        raise ValueError("掩膜尺寸与镶嵌不一致")     # 尺寸必须一致

    library_spectra = library["spectra"]             # 库光谱矩阵
    groups = library["groups"]                       # 类别分组
    candidate_idx = library["match_idx"]             # 候选匹配波段

    # ========== ② 波段选择 ==========
    # 先统计坏波段（1.44 μm 等水汽带附近，步骤2 有大量 nodata），再定参与匹配的波段
    valid_fraction, _ = estimate_band_missing(ds, mask_ds, candidate_idx, logger)   # 统计有效比例
    keep_mask = valid_fraction >= (1.0 - MAX_BAND_MISSING_FRACTION)   # 保留达标波段
    match_idx = candidate_idx[keep_mask]             # 实际参与匹配的波段
    for position, band_index in enumerate(candidate_idx):   # 记录被剔除的波段
        if keep_mask[position]:                      # 保留则跳过
            continue
        logger.info("  剔除坏波段 band %d (%.4f μm, 有效 %.1f%%)",
                    band_index + 1, wavelengths[band_index],
                    100.0 * valid_fraction[position])
    if match_idx.size < 100:                         # 可用波段过少
        raise ValueError("可用匹配波段只有 %d 个，数据质量异常" % match_idx.size)
    logger.info("[%s] 波段检查完成：候选 %d，可用 %d，用时 %.1f s",
                area_name, int(candidate_idx.size), int(match_idx.size),
                time.time() - started)
    stage_started = time.time()                      # 阶段计时起点
    n_class = len(MINERAL_EN)                        # 类别数
    n_end = library_spectra.shape[0]                 # 端元数

    # ========== ③ Al-OH 连续统窗口与库的连续统光谱 ==========
    # 连续统 / 特征窗口
    win_idx = np.where((wavelengths >= ALOH_CONTINUUM_RANGE[0]) &
                       (wavelengths <= ALOH_CONTINUUM_RANGE[1]))[0]   # 窗口内波段下标
    if win_idx.size < MIN_ALOH_BANDS:                # 窗口波段太少
        raise ValueError("Al-OH 窗口内只有 %d 个波段，不足 %d" % (win_idx.size, MIN_ALOH_BANDS))
    win_wl = wavelengths[win_idx]                    # 窗口内波长
    left_index, right_index = win_idx[0], win_idx[-1]   # 左右肩波段下标
    win_t = ((win_wl - wavelengths[left_index]) /
             (wavelengths[right_index] - wavelengths[left_index]))   # 归一化位置 0~1
    win_dwl_nm = float(np.mean(np.diff(win_wl))) * 1000.0   # 平均波段间隔（nm）
    eval_mask = (win_wl >= ALOH_FEATURE_RANGE[0]) & (win_wl <= ALOH_FEATURE_RANGE[1])   # 特征评价范围
    if eval_mask.sum() < 3:                          # 评价波段太少
        raise ValueError("Al-OH 特征评价范围内波段数不足")
    win_wl_eval = win_wl[eval_mask]                  # 评价范围波长

    # 库的连续统去除光谱（同一窗口 + 同一锚点索引）
    lib_win = library_spectra[:, win_idx].astype(np.float64)   # 库窗口光谱
    lib_left = lib_win[:, 0:1]                       # 左肩
    lib_right = lib_win[:, -1:]                      # 右肩
    lib_cont = lib_left + (lib_right - lib_left) * win_t[None, :]   # 两点直线连续统
    if np.any(lib_cont <= 0):                        # 连续统为负属异常
        logger.warning("%d 个库端元连续统为负，可能异常",
                       int((lib_cont <= 0).any(axis=1).sum()))
    lib_cr = lib_win / np.maximum(lib_cont, 1e-6)    # 库连续统去除
    lib_cr_eval = lib_cr[:, eval_mask]               # 评价范围 CR 光谱
    lib_cr_sq = np.sum(lib_cr_eval * lib_cr_eval, axis=1)   # 各端元模方

    # 匹配用库光谱：MF 用原始反射率；SAM 用连续统去除后的 Al-OH 窗口
    lib_match = library_spectra[:, match_idx].astype(np.float64)   # MF 用原始光谱
    lib_cr_unit = lib_cr_eval / np.maximum(np.sqrt(lib_cr_sq)[:, None], 1e-9)   # 单位化 CR 光谱

    # ========== ④ 背景统计与 MF 预计算 ==========
    # 背景统计
    mean, inv_cov = sample_background(ds, mask_ds, match_idx, BG_SAMPLE_PIXELS, logger)   # 背景均值/协方差逆
    mf_q = inv_cov @ (lib_match - mean).T                            # (151, n_end)
    mf_denom = np.sum((lib_match - mean) * mf_q.T, axis=1)           # (n_end,)
    logger.info("[%s] Al-OH 窗口与背景/MF 预计算完成，用时 %.1f s",
                area_name, time.time() - stage_started)
    stage_started = time.time()                      # 重置阶段计时

    logger.info("[%s] %d×%d, %d 波段，端元 %d，类别 %d；匹配波段 %d",
                area_name, width, height, wavelengths.size, n_end, n_class, match_idx.size)
    logger.info("[%s] Al-OH 窗口: %.3f~%.3f μm (%d 波段, 宽 %.1f nm/波段)",
                area_name, win_wl.min(), win_wl.max(), win_wl.size, win_dwl_nm)

    # ========== ⑤ 准备输出栅格（先写 .tmp.tif，全部成功后统一改名） ==========
    os.makedirs(out_dir, exist_ok=True)              # 确保输出目录存在
    created = {}                                     # 已创建的数据集登记表

    def make_output(name, nbands, dtype, nodata, descriptions=None):
        tmp_path = os.path.join(out_dir, name + ".tmp.tif")   # 临时文件路径
        if os.path.exists(tmp_path):                 # 清掉上次残留
            os.remove(tmp_path)
        dataset = _create_gtiff(tmp_path, width, height, nbands, dtype, geo,
                                projection, nodata, descriptions)   # 建临时 GeoTIFF
        created[name] = (dataset, tmp_path, os.path.join(out_dir, name))   # 登记
        return dataset                               # 返回数据集

    descs = ["%s %s" % (name, MINERAL_CN[name]) for name in MINERAL_EN]   # 类别波段描述
    out_class = make_output(OUT_CLASS_NAME, 1, gdal.GDT_Byte, CLASS_NODATA,
                            ["mineral class"])          # 分类图
    out_confidence = make_output(OUT_CONFIDENCE_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                                 ["classification confidence (%)"])   # 置信度
    out_depth = make_output(OUT_DEPTH_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                            ["2.2um absorption depth"])   # 吸收深度
    out_position = make_output(OUT_POSITION_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                               ["absorption center (um)"])   # 吸收中心
    out_width = make_output(OUT_WIDTH_NAME, 1, gdal.GDT_Float32, NODATA_OUT,
                            ["half-depth width (nm)"])   # 半深宽
    out_sam = make_output(OUT_SAM_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)   # SAM 角度
    out_sff = make_output(OUT_SFF_ABUNDANCE_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)   # SFF 丰度
    out_rms = make_output(OUT_SFF_RMS_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)   # SFF 残差
    out_mf = out_mtmf = None                         # MF 产物默认不建
    if ENABLE_MF and not args.no_mf:                 # 开启时才创建
        out_mf = make_output(OUT_MF_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)
        out_mtmf = make_output(OUT_MTMF_NAME, n_class, gdal.GDT_Float32, NODATA_OUT, descs)

    # 统计累加器
    class_pixels = np.zeros(n_class + 1, np.int64)   # 各类像元数
    class_depth_sum = np.zeros(n_class + 1, np.float64)   # 各类深度合计
    class_conf_sum = np.zeros(n_class + 1, np.float64)    # 各类置信度合计
    class_angle_sum = np.zeros(n_class + 1, np.float64)   # 各类角度合计
    valid_pixels = 0                                 # 有效像元累计
    clear_pixels = 0                                 # 清晰像元累计
    classified_pixels = 0                            # 分类像元累计
    depth_values = []                                # 深度样本（供分位数）
    angle_values = []                                # 角度样本（供分位数）

    def write_block(dataset, band_index, full, idx, values, rows, y0):
        block = full.copy()                          # 复制整块缓冲
        block[idx] = values                          # 只填有效像元
        block = block.reshape(rows, width)           # 还原二维
        dataset.GetRasterBand(band_index).WriteArray(block, 0, y0)   # 写盘

    logger.info("[%s] 输出栅格准备完成，用时 %.1f s", area_name, time.time() - stage_started)
    stage_started = time.time()                      # 重置阶段计时

    # ========== ⑥ 分块逐像元提取 ==========
    logger.info("[%s] 开始分块提取：块高 %d 行 × %d 波段", area_name, BLOCK_ROWS,
                wavelengths.size)
    try:
        for y0 in range(0, height, BLOCK_ROWS):      # 逐块向下扫描
            rows = min(BLOCK_ROWS, height - y0)      # 本块行数
            cube = ds.ReadAsArray(0, y0, width, rows)   # 读反射率块
            if cube is None:                         # 读取失败
                raise RuntimeError("读取镶嵌失败: y=%d h=%d" % (y0, rows))
            mask = mask_ds.ReadAsArray(0, y0, width, rows)   # 读掩膜块
            if mask.ndim == 3:                       # 去掉波段维
                mask = mask[0]
            clear_pixels += int(np.count_nonzero(mask == 1))   # 累计清晰像元
            if USE_CLOUD_FILLED:                     # 允许云填充像元
                pixel_ok = (mask == 1) | (mask == 2) | (mask == 3)
            else:
                pixel_ok = (mask == 1)               # 只用清晰像元
            npix = width * rows                      # 本块像元数
            pixel_ok = pixel_ok.reshape(-1)          # 拉平

            match_block = cube[match_idx].reshape(match_idx.size, -1).astype(np.float64)   # 匹配波段块
            win_block = cube[win_idx].reshape(win_idx.size, -1).astype(np.float64)         # 窗口波段块
            match_finite = _valid_pixels(match_block)   # 匹配波段有效掩膜
            missing = np.count_nonzero(~match_finite, axis=0)   # 每像元缺测波段数
            low = np.where(match_finite, match_block, np.inf).min(axis=0)    # 每像元最小值
            high = np.where(match_finite, match_block, -np.inf).max(axis=0)  # 每像元最大值
            pixel_ok &= (missing <= MAX_MISSING_BANDS_MATCH)   # 缺测不能过多
            pixel_ok &= (low > REFL_MIN) & (high < REFL_MAX)   # 反射率需在范围内
            pixel_ok &= _valid_pixels(win_block).all(axis=0)   # 窗口波段须全有效
            idx = np.where(pixel_ok)[0]              # 有效像元下标
            if idx.size == 0:                        # 本块无有效像元
                continue

            # 缺测波段按 0 填充（只影响 151 个匹配波段里的 1~2 个，对 SAM/MF 影响很小）
            match_filled = np.where(match_finite, match_block, 0.0)   # 缺测补 0

            # 连续统去除：用窗口左右肩两点连一条直线，逐波段按位置线性插值。
            # cr = 实测反射率 / 连续统；Al-OH 吸收越强，cr 在 2.2 μm 处越低。
            win_valid = win_block[:, idx]            # 窗口光谱（仅有效像元）
            left = win_valid[0][None, :]             # 左肩反射率
            right = win_valid[-1][None, :]           # 右肩反射率
            continuum = left + (right - left) * win_t[:, None]   # 逐像元连续统
            cr_ok = continuum.min(axis=0) > MIN_CONTINUUM_VALUE   # 连续统不能过低
            idx = idx[cr_ok]                         # 再筛一次像元
            if idx.size == 0:                        # 无剩余像元
                continue
            match_valid = match_filled[:, idx]       # 匹配光谱（仅有效像元）
            win_valid = win_valid[:, cr_ok]          # 窗口光谱同步筛选
            continuum = continuum[:, cr_ok]          # 连续统同步筛选
            cr = win_valid / continuum               # 连续统去除光谱

            # 吸收三要素：深度 = 1-cr 的最大值；位置 = 最大值对应波长；
            # 宽度 = 深度过半的波段数 × 波段间隔（半深宽，近似 FWHM）
            depth_eval = 1.0 - cr[eval_mask]         # 特征范围内吸收深度
            depth = np.clip(depth_eval.max(axis=0), -1.0, 1.0)   # 取最大深度并限幅
            argmax = depth_eval.argmax(axis=0)       # 最大深度位置
            position = win_wl_eval[argmax]           # 对应中心波长
            half = depth_eval >= (0.5 * depth[None, :])   # 超过半深的波段
            width_nm = half.sum(axis=0) * win_dwl_nm      # 半深宽（nm）

            # ── SAM：连续统去除后的 Al-OH 窗口光谱角（消除反照率/大气残差影响） ──
            crn = cr[eval_mask]                      # 评价范围 CR 光谱
            crn_norm = np.linalg.norm(crn, axis=0)   # 像元光谱模长
            good_norm = crn_norm > 1e-9              # 模长有效标记
            crn_unit = crn / np.maximum(crn_norm[None, :], 1e-9)   # 单位化像元光谱
            cos = lib_cr_unit @ crn_unit             # 与各端元的余弦相似度
            np.clip(cos, -1.0, 1.0, out=cos)         # 数值限幅
            angle_end = np.degrees(np.arccos(cos))                 # (n_end, nv)

            # ── SFF（连续统去除光谱特征拟合） ──
            # 把库光谱当作一条形状模板，用一个最小二乘比例系数 scale 去拟合像元光谱：
            #   scale = <像元, 库> / <库, 库>
            # 残差 rms 表示“形状不像”的程度；scale 可视为该类矿物的相对丰度代理。
            cr_eval = cr[eval_mask]                  # 评价范围 CR 光谱
            dot = cr_eval.T @ lib_cr_eval.T                        # (nv, n_end)
            scale = dot / np.maximum(lib_cr_sq[None, :], 1e-12)    # (nv, n_end)
            cr_sq = np.sum(cr_eval * cr_eval, axis=0)              # (nv,)
            resid = cr_sq[:, None] - (scale ** 2) * lib_cr_sq[None, :]   # 残差平方
            rms = np.sqrt(np.maximum(resid, 0.0) / lib_cr_eval.shape[1])   # 均方根残差

            # ── 类别级汇总 ──
            angle_class = _class_reduce(angle_end, groups, np.min, np.inf)   # 类内最小角
            rms_class = _class_reduce(rms.T, groups, np.min, np.inf)         # 类内最小残差
            abundance_class = _class_reduce(scale.T, groups, np.max, 0.0)    # 类内最大丰度
            abundance_class = np.clip(abundance_class, 0.0, 1.0)   # 丰度限幅 0~1

            # ── 分类判决 ──
            # 先取每类中 CR-SAM 角最小（--class-by sff 时取残差最小）的端元，
            # 再同时满足：吸收深度达标、光谱角达标、SFF 残差达标 → 归入该矿物类。
            best_index = np.argmin(angle_class, axis=0)   # 每像元最优类别
            arange = np.arange(angle_class.shape[1])      # 像元序号
            best_angle = angle_class[best_index, arange]  # 最优类别角度
            best_rms = rms_class[best_index, arange]      # 最优类别残差
            best_abundance = abundance_class[best_index, arange]   # 最优类别丰度

            if args.class_by == "sff":               # 改用 SFF 残差判类
                best_index = np.argmin(rms_class, axis=0)          # 残差最小的类别
                best_angle = angle_class[best_index, arange]       # 对应角度
                best_rms = rms_class[best_index, arange]           # 对应残差
                best_abundance = abundance_class[best_index, arange]   # 对应丰度

            usable = (depth >= args.min_depth) & good_norm   # 吸收深度达标且模长有效
            classify = usable & (best_angle <= args.max_angle) & (best_rms <= args.max_rms)   # 三条件同时满足
            class_values = np.where(classify, best_index + 1, 0).astype(np.uint8)   # 0 为未分类
            confidence = np.where(usable,
                                  np.clip(100.0 * (1.0 - best_angle / max(args.max_angle, 1e-6)),
                                          0.0, 100.0), 0.0)   # 置信度百分数

            # ── MF / MTMF（辅助产物，不参与默认判决） ──
            # MF 分数：端元与背景均值之差在背景协方差白化空间中的投影，
            #   越大越像该端元；MTMF 不可行性：白化空间里 MF 无法解释的残差长度，
            #   越大说明该像元越不可能由该类端元与背景混合而成。
            mf_class = infeas_class = None           # MF 结果先置空
            if out_mf is not None:                   # 需要 MF 时才计算
                centered = (match_valid - mean[:, None]).T          # (nv, 151)
                mf_score = centered @ mf_q / np.maximum(mf_denom[None, :], 1e-12)   # MF 分数
                mf_class = _class_reduce(mf_score.T, groups, np.max, 0.0).T   # 类内最大 MF
                maha_sq = np.maximum(np.sum((centered @ inv_cov) * centered, axis=1), 0.0)   # 马氏距离平方
                infeas = np.sqrt(np.maximum(
                    maha_sq[:, None] - (mf_score ** 2) * mf_denom[None, :], 0.0))   # 不可行性
                infeas_class = _class_reduce(infeas.T, groups, np.min, np.inf).T   # 类内最小

            # ── 写块 ──
            class_block = np.full(npix, CLASS_NODATA, np.uint8)   # 分类块先填 nodata
            class_block[idx] = class_values          # 有效像元写类别
            out_class.GetRasterBand(1).WriteArray(class_block.reshape(rows, width), 0, y0)   # 写分类图

            for dataset, values, nodata in (
                    (out_confidence, confidence, NODATA_OUT),
                    (out_depth, depth, NODATA_OUT)):
                buffer = np.full(npix, nodata, np.float32)   # 整块先填 nodata
                buffer[idx] = np.where(np.isfinite(values), values, nodata)   # 填有效值
                dataset.GetRasterBand(1).WriteArray(buffer.reshape(rows, width), 0, y0)

            # 吸收位置/宽度只在真正有吸收（深度≥阈值）的像元上输出，避免“无特征”像元
            # 在窗口边缘取值造成整幅图都是 2.28 μm 的假象
            feature_ok = usable & (depth >= args.min_depth)   # 真正有吸收的像元
            for dataset, values in ((out_position, position), (out_width, width_nm)):
                buffer = np.full(npix, NODATA_OUT, np.float32)   # 整块先填 nodata
                buffer[idx] = np.where(feature_ok & np.isfinite(values),
                                       values, NODATA_OUT)   # 只写有特征像元
                dataset.GetRasterBand(1).WriteArray(buffer.reshape(rows, width), 0, y0)

            for band_index in range(n_class):        # 逐类别写多波段成果
                buffer = np.full(npix, NODATA_OUT, np.float32)   # SAM 角度缓冲
                buffer[idx] = np.where(np.isfinite(angle_class[band_index]),
                                       angle_class[band_index], NODATA_OUT)
                out_sam.GetRasterBand(band_index + 1).WriteArray(
                    buffer.reshape(rows, width), 0, y0)
                buffer = np.full(npix, NODATA_OUT, np.float32)   # SFF 丰度缓冲
                buffer[idx] = np.where(np.isfinite(abundance_class[band_index]),
                                       abundance_class[band_index], NODATA_OUT)
                out_sff.GetRasterBand(band_index + 1).WriteArray(
                    buffer.reshape(rows, width), 0, y0)
                buffer = np.full(npix, NODATA_OUT, np.float32)   # SFF 残差缓冲
                buffer[idx] = np.where(np.isfinite(rms_class[band_index]),
                                       rms_class[band_index], NODATA_OUT)
                out_rms.GetRasterBand(band_index + 1).WriteArray(
                    buffer.reshape(rows, width), 0, y0)
                if out_mf is not None:               # MF 开启时才写
                    buffer = np.full(npix, NODATA_OUT, np.float32)   # MF 缓冲
                    buffer[idx] = np.where(np.isfinite(mf_class[:, band_index]),
                                           mf_class[:, band_index], NODATA_OUT)
                    out_mf.GetRasterBand(band_index + 1).WriteArray(
                        buffer.reshape(rows, width), 0, y0)
                    buffer = np.full(npix, NODATA_OUT, np.float32)   # MTMF 缓冲
                    buffer[idx] = np.where(np.isfinite(infeas_class[:, band_index]),
                                           infeas_class[:, band_index], NODATA_OUT)
                    out_mtmf.GetRasterBand(band_index + 1).WriteArray(
                        buffer.reshape(rows, width), 0, y0)

            # ── 统计 ──
            class_pixels += np.bincount(class_values.astype(np.int64), minlength=n_class + 1)   # 类别计数
            valid_pixels += int(idx.size)            # 累计有效像元
            classified_pixels += int(np.count_nonzero(classify))   # 累计分类像元
            for class_index in range(n_class + 1):   # 逐类累计指标
                mask_class = class_values == class_index   # 该类像元掩膜
                if np.any(mask_class):               # 该类有像元
                    class_depth_sum[class_index] += float(depth[mask_class].sum())   # 累计深度
                    class_conf_sum[class_index] += float(confidence[mask_class].sum())   # 累计置信度
                    class_angle_sum[class_index] += float(best_angle[mask_class].sum())   # 累计角度
            depth_values.append(depth.astype(np.float32))     # 收集深度样本
            angle_values.append(best_angle.astype(np.float32))   # 收集角度样本

            if (y0 // BLOCK_ROWS) % 4 == 0 or y0 + rows >= height:   # 每 4 块打印一次进度
                logger.info("  [%s] 行 %d/%d（%.0f%%），有效像元累计 %d",
                            area_name, y0 + rows, height,
                            100.0 * (y0 + rows) / height, valid_pixels)
    except Exception:                                # 任一步失败则清理临时文件
        out_class = out_confidence = out_depth = out_position = out_width = None   # 断开引用
        out_sam = out_sff = out_rms = out_mf = out_mtmf = None   # 断开引用，释放文件句柄
        for name in list(created):                   # 逐个关闭并删除临时文件
            dataset, tmp_path, final_path = created.pop(name)
            try:
                dataset.FlushCache()                 # 先落盘
            except Exception:
                pass
            dataset = None                           # 释放句柄
            try:
                if os.path.exists(tmp_path):         # 删除残留临时文件
                    os.remove(tmp_path)
            except OSError:
                pass
        raise                                        # 异常继续向上抛
    finally:
        ds = None                                    # 释放输入数据集
        mask_ds = None                               # 释放掩膜数据集

    logger.info("[%s] 分块提取完成，用时 %.1f s", area_name, time.time() - stage_started)

    # 关闭并原子替换（先清掉局部变量引用，Windows 下文件句柄才能释放）
    n_written = len(created)                         # 成果文件数
    out_class = out_confidence = out_depth = out_position = out_width = None   # 断开引用
    out_sam = out_sff = out_rms = out_mf = out_mtmf = None
    for name in list(created):                       # 逐个关闭并改名
        dataset, tmp_path, final_path = created.pop(name)
        dataset.FlushCache()                         # 落盘
        dataset = None                               # 释放句柄
        os.replace(tmp_path, final_path)             # 原子改名到正式文件

    logger.info("[%s] 写出 %d 个成果文件", area_name, n_written)

    # 汇总/图例/信息
    rows = []                                        # 逐类统计行
    for class_index in range(n_class + 1):           # 含未分类共 n_class+1 类
        count = int(class_pixels[class_index])       # 该类像元数
        if class_index == 0:                         # 0 号为未分类
            name_en, name_cn = "unclassified", "未分类"
        else:
            name_en = MINERAL_EN[class_index - 1]    # 英文名
            name_cn = MINERAL_CN[name_en]            # 中文名
        mean_depth = (class_depth_sum[class_index] / count) if count else 0.0   # 平均深度
        mean_conf = (class_conf_sum[class_index] / count) if count else 0.0     # 平均置信度
        mean_angle = (class_angle_sum[class_index] / count) if count else 0.0   # 平均角度
        rows.append({
            "class_id": class_index,                 # 类别号
            "mineral_en": name_en,                   # 矿物英文名
            "mineral_cn": name_cn,                   # 矿物中文名
            "pixels": count,                         # 像元数
            "area_km2": count * pixel_area_m2 / 1.0e6,   # 面积（km²）
            "percent_of_valid": (100.0 * count / valid_pixels) if valid_pixels else 0.0,   # 占有效像元比
            "mean_depth": mean_depth,                # 平均吸收深度
            "mean_confidence": mean_conf,            # 平均置信度
            "mean_sam_deg": mean_angle,              # 平均 SAM 角
        })
    write_class_summary(os.path.join(out_dir, OUT_SUMMARY_NAME), rows)   # 写统计 CSV
    write_legend(os.path.join(out_dir, OUT_LEGEND_NAME), library)        # 写图例 CSV

    depth_stats = np.concatenate(depth_values) if depth_values else np.array([])   # 汇总深度样本
    angle_stats = np.concatenate(angle_values) if angle_values else np.array([])   # 汇总角度样本
    try:
        srs = osr.SpatialReference(wkt=projection) if projection else None   # 解析投影
        epsg = int(srs.GetAuthorityCode(None) or 0) if srs is not None else 0   # 取 EPSG 码
    except Exception:                                # 解析失败则记 0
        epsg = 0
    info = {
        "area": area_name,                           # 区域名
        "source": in_tif,                            # 输入路径
        "grid": {"width": int(width), "height": int(height),
                 "geotransform": [float(v) for v in geo], "epsg": epsg},   # 栅格信息
        "bands": int(wavelengths.size),              # 波段数
        "matching_bands": int(match_idx.size),       # 参与匹配波段数
        "sensor_response": response_source,          # 传感器响应来源
        "library": {"root": os.path.abspath(args.lib_root),
                    "endmembers": int(n_end), "classes": n_class},   # 库信息
        "parameters": {
            "aloh_continuum_range_um": list(ALOH_CONTINUUM_RANGE),   # 连续统窗口
            "aloh_feature_range_um": list(ALOH_FEATURE_RANGE),       # 特征范围
            "min_aloh_depth": float(args.min_depth),                 # 深度阈值
            "sam_max_deg": float(args.max_angle),                    # SAM 阈值
            "sff_max_rms": float(args.max_rms),                      # SFF 阈值
            "class_by": args.class_by,                               # 分类依据
            "mf_enabled": bool(out_mf is not None),                  # 是否启用 MF
            "use_cloud_filled": bool(USE_CLOUD_FILLED),              # 是否用云填充像元
        },
        "coverage": {
            "clear_pixels": int(clear_pixels),       # 清晰像元
            "valid_pixels": int(valid_pixels),       # 有效像元
            "classified_pixels": int(classified_pixels),   # 分类像元
            "classified_percent_of_valid": (100.0 * classified_pixels / valid_pixels) if valid_pixels else 0.0,
        },
        "statistics": {
            "depth_percentiles": _percentiles(depth_stats, [1, 5, 50, 95, 99]),   # 深度分位
            "sam_percentiles": _percentiles(angle_stats, [1, 5, 50, 95, 99]),     # 角度分位
        },
        "created": datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),   # 生成时间
    }
    with open(final_info + ".tmp", "w", encoding="utf-8") as handle:
        json.dump(info, handle, ensure_ascii=False, indent=2)   # 先写临时 JSON
    os.replace(final_info + ".tmp", final_info)      # 原子改名

    if MAKE_QC:                                      # 需要质检图时
        try:
            make_area_qc(out_dir, area_name, in_tif, wavelengths, rows, logger)
        except Exception as exc:                     # 出图失败不影响主成果
            logger.warning("[%s] 质检图失败: %s", area_name, exc, exc_info=True)

    elapsed = time.time() - started                  # 总耗时
    logger.info("[%s] 完成: 有效像元 %d，分类像元 %d（%.1f%% of valid），用时 %.1f s",
                area_name, valid_pixels, classified_pixels,
                100.0 * classified_pixels / max(1, valid_pixels), elapsed)
    for row in rows[1:]:                             # 逐类打印结果
        if row["pixels"]:
            logger.info("    %-18s %8d 像元  %8.2f km2  %6.2f%%  深度 %.3f  SAM %.1f°",
                        row["mineral_en"], row["pixels"], row["area_km2"],
                        row["percent_of_valid"], row["mean_depth"], row["mean_sam_deg"])
    return {"area": area_name, "status": "成功",
            "message": "分类像元 %d" % classified_pixels,
            "elapsed": elapsed, "classes": rows,
            "valid_pixels": int(valid_pixels),
            "classified_pixels": int(classified_pixels)}   # 返回区域结果


def _percentiles(values, percents):
    """只对有限值求分位数；样本为空时返回 [None, ...]。"""
    if values.size == 0:                             # 没有样本
        return [None for _ in percents]
    finite = values[np.isfinite(values)]             # 只留有限值
    if finite.size == 0:                             # 全是 NaN
        return [None for _ in percents]
    return [float(v) for v in np.percentile(finite, percents)]   # 求分位数


def write_class_summary(path, rows):
    """写出 aloh_summary.csv：逐类像元数、面积、占比与均值指标。"""
    fields = ["class_id", "mineral_en", "mineral_cn", "pixels", "area_km2",
              "percent_of_valid", "mean_depth", "mean_confidence", "mean_sam_deg"]   # 列顺序
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)   # 字典写入器
        writer.writeheader()                         # 写表头
        for row in rows:                             # 逐类写一行
            out = dict(row)                          # 复制一行
            out["area_km2"] = "%.4f" % out["area_km2"]              # 面积保留 4 位
            out["percent_of_valid"] = "%.3f" % out["percent_of_valid"]   # 占比保留 3 位
            out["mean_depth"] = "%.5f" % out["mean_depth"]          # 深度保留 5 位
            out["mean_confidence"] = "%.2f" % out["mean_confidence"]   # 置信度 2 位
            out["mean_sam_deg"] = "%.3f" % out["mean_sam_deg"]      # 角度保留 3 位
            writer.writerow(out)


def write_legend(path, library):
    """写出 classification_legend.csv：类别号 ↔ 矿物名 ↔ 颜色 ↔ 端元数。"""
    class_ids = library["class_ids"]                 # 端元类别号
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)                  # CSV 写入器
        writer.writerow(["class_id", "mineral_en", "mineral_cn", "color",
                         "endmembers", "samples"])   # 表头
        for class_index, name_en in enumerate(MINERAL_EN, start=1):   # 逐类写
            idx = np.where(class_ids == class_index)[0]   # 该类端元下标
            samples = ";".join(library["meta"][i]["sample"] for i in idx)   # 样本名拼接
            writer.writerow([class_index, name_en, MINERAL_CN[name_en],
                             MINERAL_COLOR[name_en], int(idx.size), samples])


def make_area_qc(out_dir, area_name, mosaic_path, wavelengths, rows, logger):
    """生成 preview_quality.png：真彩色 / 分类 / 深度 / 吸收位置 / 置信度 / 面积柱状图。

    为控制内存与出图时间，先按 QC_MAX_SIZE 降采样读取各成果波段。
    """
    import matplotlib                              # 绘图库
    matplotlib.use("Agg")                          # 无界面后端
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm   # 分类色表/分段

    class_path = os.path.join(out_dir, OUT_CLASS_NAME)          # 分类图
    depth_path = os.path.join(out_dir, OUT_DEPTH_NAME)          # 深度图
    pos_path = os.path.join(out_dir, OUT_POSITION_NAME)         # 吸收中心
    conf_path = os.path.join(out_dir, OUT_CONFIDENCE_NAME)      # 置信度
    sam_path = os.path.join(out_dir, OUT_SAM_NAME)              # SAM 角度

    dataset = gdal.Open(class_path)                 # 打开分类图
    if dataset is None:                             # 打开失败
        raise RuntimeError("无法打开分类成果: %s" % class_path)
    width, height = dataset.RasterXSize, dataset.RasterYSize   # 影像尺寸
    geo = dataset.GetGeoTransform()                 # 地理变换
    dataset = None                                  # 释放句柄
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))   # 降采样倍数
    bw, bh = max(1, width // scale), max(1, height // scale)    # 降采样后尺寸

    def read_band(path, band=1, dtype=np.float64):
        data = gdal.Open(path)                      # 打开成果文件
        if data is None:                            # 打开失败返回 None
            return None
        if data.RasterCount == 1:                   # 单波段直接读
            array = data.ReadAsArray(buf_xsize=bw, buf_ysize=bh)
        else:
            band_handle = data.GetRasterBand(band)  # 多波段取指定波段
            array = band_handle.ReadAsArray(buf_xsize=bw, buf_ysize=bh)
            band_handle = None                      # 释放波段句柄
        data = None                                 # 释放数据集
        return array.astype(dtype)                  # 转目标类型返回

    class_array = read_band(class_path, 1, np.int32)   # 分类（整型）
    depth_array = read_band(depth_path)             # 吸收深度
    pos_array = read_band(pos_path)                 # 吸收中心
    conf_array = read_band(conf_path)               # 置信度

    sam_data = gdal.Open(sam_path)                  # 打开 SAM 多波段图
    sam_best = None
    if sam_data is not None:                        # 打开成功才处理
        sam_best = np.full((bh, bw), np.inf, np.float64)   # 取各类最小角
        for band_index in range(1, sam_data.RasterCount + 1):   # 逐类别波段
            band_handle = sam_data.GetRasterBand(band_index)
            band = band_handle.ReadAsArray(buf_xsize=bw, buf_ysize=bh).astype(np.float64)
            band_handle = None                      # 释放波段句柄
            band = np.where(band < NODATA_THRESHOLD, np.inf, band)   # nodata 转 inf
            sam_best = np.minimum(sam_best, band)   # 逐像元取最小
        sam_data = None                             # 释放数据集
        sam_best[~np.isfinite(sam_best)] = np.nan   # inf 转 NaN 便于出图

    # 真彩色
    color_bands = [_nearest_band(wavelengths, TRUE_COLOR_TARGETS[key]) + 1
                   for key in ("R", "G", "B")]      # 红绿蓝对应波段（1 起）
    rgb = None
    mosaic = gdal.Open(mosaic_path)                 # 打开镶嵌做底图
    if mosaic is not None:                          # 打开成功才画
        channels = []                               # 三个通道数据
        for band_index in color_bands:              # 逐通道读取
            band_handle = mosaic.GetRasterBand(band_index)
            channels.append(band_handle.ReadAsArray(
                buf_xsize=bw, buf_ysize=bh).astype(np.float64))
            band_handle = None                      # 释放波段句柄
        mosaic = None                               # 释放数据集
        valid = np.ones(channels[0].shape, bool)    # 有效像元掩膜
        for channel in channels:                    # 三个通道都有效才算
            valid &= np.isfinite(channel) & (channel > -9000)
        rgb = np.ones(channels[0].shape + (3,), np.float64)   # 白底
        for position, channel in enumerate(channels):   # 逐通道拉伸
            rgb[..., position] = np.where(valid, _stretch(channel, valid), 1.0)

    colors = ["#dddddd"] + [MINERAL_COLOR[name] for name in MINERAL_EN]   # 未分类灰 + 各类颜色
    cmap = ListedColormap(colors)                    # 离散色表
    bounds = np.arange(-0.5, len(colors) + 0.5, 1.0)   # 分类分段边界
    extent = [geo[0] / 1000.0, (geo[0] + geo[1] * width) / 1000.0,
              (geo[3] + geo[5] * height) / 1000.0, geo[3] / 1000.0]   # 投影范围（km）
    class_display = np.where(class_array == CLASS_NODATA, 0, class_array)   # nodata 归为未分类

    figure, axes = plt.subplots(2, 3, figsize=(18, 10.5))   # 2×3 子图
    axes = axes.ravel()                              # 展平成列表

    if rgb is not None:                              # 有底图才画真彩色
        axes[0].imshow(rgb, extent=extent, origin="upper", interpolation="nearest")
    axes[0].set_title("True color (R %.2f / G %.2f / B %.2f um)"
                      % tuple(wavelengths[band - 1] for band in color_bands), fontsize=10)

    axes[1].imshow(class_display, extent=extent, origin="upper", cmap=cmap,
                   norm=BoundaryNorm(bounds, cmap.N), interpolation="nearest")   # 分类图
    axes[1].set_title("mineral class", fontsize=10)
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[index])
               for index in range(1, len(colors))]   # 图例色块
    axes[1].legend(handles, MINERAL_EN, fontsize=6.5, ncol=2, loc="lower left",
                   framealpha=0.85)

    depth_display = np.where(depth_array > -9000, depth_array, np.nan)   # nodata 转 NaN
    image = axes[2].imshow(depth_display, extent=extent, origin="upper", cmap="viridis",
                           vmin=-0.5, vmax=0.5, interpolation="nearest")   # 深度图
    axes[2].set_title("2.2 um absorption depth (clip -0.5~0.5)", fontsize=10)
    figure.colorbar(image, ax=axes[2], fraction=0.035, pad=0.02)   # 色标

    pos_display = np.where(pos_array > -9000, pos_array, np.nan)   # nodata 转 NaN
    image = axes[3].imshow(pos_display, extent=extent, origin="upper", cmap="turbo",
                           vmin=float(ALOH_FEATURE_RANGE[0]),
                           vmax=float(ALOH_FEATURE_RANGE[1]), interpolation="nearest")   # 吸收中心
    axes[3].set_title("Absorption center (um)", fontsize=10)
    figure.colorbar(image, ax=axes[3], fraction=0.035, pad=0.02)

    conf_display = np.where(conf_array > -9000, conf_array, np.nan)   # nodata 转 NaN
    image = axes[4].imshow(conf_display, extent=extent, origin="upper", cmap="magma",
                           vmin=0, vmax=100, interpolation="nearest")   # 置信度
    axes[4].set_title("Classification confidence (%)", fontsize=10)
    figure.colorbar(image, ax=axes[4], fraction=0.035, pad=0.02)

    labels = [row["mineral_en"] for row in rows[1:]]   # 矿物名
    values = [row["area_km2"] for row in rows[1:]]     # 面积
    order = np.argsort(values)[::-1]                 # 面积降序
    axes[5].barh([labels[i] for i in order], [values[i] for i in order],
                 color=[MINERAL_COLOR[labels[i]] for i in order])   # 面积柱状图
    axes[5].set_xlabel("Area (km2)")
    axes[5].set_title("mineral class area", fontsize=10)
    axes[5].tick_params(labelsize=7.5)
    axes[5].grid(axis="x", alpha=0.3)

    for axis in axes[:5]:                            # 前 5 图统一坐标轴标签
        axis.set_xlabel("Easting (km)", fontsize=8)
        axis.set_ylabel("Northing (km)", fontsize=8)
        axis.tick_params(labelsize=7)

    try:
        depth_p50 = float(np.nanpercentile(depth_display, 50))   # 深度中位数
        angle_p50 = float(np.nanpercentile(sam_best, 50))        # 角度中位数
    except Exception:                                # 全 NaN 等异常
        depth_p50 = angle_p50 = float("nan")
    figure.suptitle("%s | mineral extraction | depth p50 %.3f | SAM p50 %.1f deg"
                    % (area_name, depth_p50, angle_p50), fontsize=12)   # 总标题
    figure.tight_layout(rect=(0, 0, 1, 0.96))        # 紧凑布局
    path = os.path.join(out_dir, OUT_QC_NAME)        # 质检图路径
    figure.savefig(path, dpi=QC_DPI, bbox_inches="tight")   # 保存
    plt.close(figure)                                # 释放画布
    logger.info("[%s] 质检图: %s", area_name, path)


# =================================================================
# 主流程
# =================================================================

def write_areas_summary(path, results):
    """把本次批处理每个区域的结果写进 _catalog/areas_summary.csv。"""
    fields = ["area", "status", "valid_pixels", "classified_pixels",
              "dominant_mineral", "dominant_percent", "message", "elapsed"]   # 列顺序
    existing = {}                                    # 已有汇总（按区名索引）
    if os.path.isfile(path):                         # 已有文件则先读取
        try:
            with open(path, "r", encoding="utf-8-sig", newline="") as handle:
                for row in csv.DictReader(handle):   # 逐行读入
                    existing[row.get("area", "")] = row
        except OSError:                              # 读取失败则当空表
            existing = {}
    for result in results:                           # 合并本次结果
        classes = result.get("classes") or []        # 分类统计
        dominant = None                              # 主导矿物
        if len(classes) > 1:                         # 有分类行才找主导
            dominant = max(classes[1:], key=lambda row: row["pixels"])
        existing[result["area"]] = {
            "area": result["area"],                  # 区域名
            "status": result["status"],              # 状态
            "valid_pixels": str(result.get("valid_pixels", "")),   # 有效像元
            "classified_pixels": str(result.get("classified_pixels", "")),   # 分类像元
            "dominant_mineral": dominant["mineral_en"] if dominant else "",   # 主导矿物
            "dominant_percent": ("%.2f" % dominant["percent_of_valid"]) if dominant else "",
            "message": result.get("message", ""),    # 说明
            "elapsed": "%.1f" % result.get("elapsed", 0.0),   # 耗时
        }
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)   # 字典写入器
        writer.writeheader()                         # 写表头
        for area in sorted(existing):                # 按区名排序输出
            writer.writerow({key: existing[area].get(key, "") for key in fields})


def main(argv=None):
    """主流程：解析参数 → 扫描区域 → 建/读库 → 逐区处理 → 汇总。

    返回码：0 全部成功；1 有区域失败；2 前置条件不满足（目录/建库等）。
    """
    args = parse_args(argv)                          # 解析命令行参数
    out_root = os.path.abspath(args.out_root)        # 输出根目录绝对路径
    os.makedirs(out_root, exist_ok=True)             # 确保输出目录存在
    catalog_dir = os.path.join(out_root, CATALOG_DIR_NAME)   # 日志/汇总目录
    os.makedirs(catalog_dir, exist_ok=True)          # 确保目录存在
    log_path = args.log or os.path.join(catalog_dir, "step4_batch.log")   # 日志路径
    logger = setup_logging(log_path, args.verbose)   # 初始化日志
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)    # 设置 GDAL 缓存

    logger.info("=" * 72)
    logger.info("ZY1E AHSI 步骤4：云母-绿泥石等蚀变矿物提取（splib07a + SAM/SFF/MF）")
    logger.info("步骤3 目录: %s", os.path.abspath(args.in_root))
    logger.info("输出目录  : %s", out_root)
    logger.info("光谱库    : %s", args.lib_root)
    logger.info("目标矿物类: %s", ", ".join(MINERAL_EN))
    logger.info("日志      : %s", log_path)
    logger.info("=" * 72)

    try:
        areas = discover_areas(args.in_root, args.area)   # 扫描可处理区域
    except Exception as exc:                         # 扫描失败属前置错误
        logger.error("扫描步骤3 成果失败: %s", exc)
        return 2
    if args.limit > 0:                               # --limit 限制区域数
        areas = areas[:args.limit]
    if not areas:                                    # 没有可处理区域
        logger.error("在 %s 下未找到可处理的区域（需要 <区名>/mosaic/%s 或 <区名>/*_clip.tif）",
                     args.in_root, IN_MOSAIC_NAME)
        return 2
    logger.info("待处理区域: %d (%s)", len(areas), ", ".join(areas))

    try:
        first_inputs = resolve_area_inputs(args.in_root, areas[0], out_root,
                                           False, logger, build=False)   # 试解析首个区域
    except Exception as exc:                         # 解析失败属前置错误
        logger.error("解析区域输入失败: %s", exc)
        return 2
    wavelengths = first_inputs["wavelengths"]        # 取波长数组
    band_fwhm, response_source = load_sensor_fwhm(TAR_DIR, wavelengths, logger)   # 读传感器 FWHM

    if args.dry_run:                                 # 只打印计划
        logger.info("-" * 72)
        logger.info("dry-run：不执行实际处理")
        library_dir = os.path.join(out_root, LIB_DIR_NAME)   # 匹配库目录
        logger.info("匹配库目录: %s（%s）", library_dir,
                    "已存在，将复用" if os.path.isfile(os.path.join(library_dir, "ahsi_library.npy"))
                    else "尚未构建，正式运行时构建")
        for name in areas:                           # 逐区域打印输入
            try:
                inputs = resolve_area_inputs(args.in_root, name, out_root,
                                             False, logger, build=False)
            except Exception as exc:                 # 输入不可用
                logger.error("  %-10s 输入不可用: %s", name, exc)
                continue
            if inputs["source"] == "mosaic":         # 已有镶嵌
                logger.info("  %-10s <- %s", name, inputs["refl_path"])
            else:                                    # 需合并裁剪成果
                logger.info("  %-10s <- %d 景裁剪成果（将合并为 %s）",
                            name, len(inputs["scenes"]), inputs["refl_path"])
        logger.info("阈值: min_depth=%.3f, max_angle=%.1f°, max_rms=%.3f, class_by=%s",
                    args.min_depth, args.max_angle, args.max_rms, args.class_by)
        return 0                                     # dry-run 正常结束

    library_dir = os.path.join(out_root, LIB_DIR_NAME)   # 匹配库目录
    logger.info("构建/加载 AHSI 匹配库: %s", library_dir)
    try:
        library = build_library(args.lib_root, wavelengths, band_fwhm, library_dir,
                                logger, force=args.rebuild_library)   # 建库或复用
    except Exception:                                # 建库失败属前置错误
        logger.exception("构建 AHSI 匹配库失败（--lib=%s）", args.lib_root)
        logger.error("排查建议：确认光谱库目录存在、波长数(%d)与库网格匹配；日志: %s",
                     wavelengths.size, log_path)
        return 2

    if args.library_only:                            # 只建库模式
        logger.info("-" * 72)
        logger.info("匹配库: %d 端元 × %d 波段", library["spectra"].shape[0],
                    library["spectra"].shape[1])
        for class_index, name_en in enumerate(MINERAL_EN, start=1):   # 逐类统计端元数
            count = int(np.count_nonzero(library["class_ids"] == class_index))
            logger.info("  %-16s %d 个端元", name_en, count)
        logger.info("库目录: %s", library_dir)
        return 0                                     # 建库完成即退出

    results = []                                     # 各区处理结果
    started = time.time()                            # 批处理开始时间
    for index, area in enumerate(areas, start=1):    # 逐区域处理
        logger.info("-" * 72)
        logger.info("(%d/%d) %s", index, len(areas), area)
        try:
            inputs = resolve_area_inputs(args.in_root, area, out_root,
                                         args.overwrite, logger)   # 解析输入
            results.append(process_area(area, inputs, out_root,
                                        library, response_source, args, logger))   # 分类提取
        except Exception as exc:                     # 单区失败不影响其他区
            logger.exception("[%s] 处理失败", area)
            results.append({"area": area, "status": "失败", "message": str(exc),
                            "elapsed": 0.0, "classes": None,
                            "valid_pixels": 0, "classified_pixels": 0})

    elapsed = time.time() - started                  # 总耗时
    logger.info("=" * 72)
    success = [r for r in results if r["status"] == "成功"]   # 成功列表
    skipped = [r for r in results if r["status"] == "跳过"]   # 跳过列表
    failed = [r for r in results if r["status"] == "失败"]    # 失败列表
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %.1f s",
                len(success), len(skipped), len(failed), elapsed)
    for result in results:                           # 逐条打印结果
        logger.info("  [%s] %s - %s", result["status"], result["area"],
                    result.get("message", ""))
    if results:                                      # 有结果才写汇总
        write_areas_summary(os.path.join(catalog_dir, "areas_summary.csv"), results)
    if failed:                                       # 有失败则返回 1
        logger.error("失败区域: %s", ", ".join(r["area"] for r in failed))
        logger.error("排查请查看日志: %s", log_path)
        return 1
    logger.info("成果目录: %s", out_root)
    logger.info("日志文件: %s", log_path)
    return 0                                         # 全部成功


if __name__ == "__main__":
    sys.exit(main())                                 # 运行主流程并按返回码退出
