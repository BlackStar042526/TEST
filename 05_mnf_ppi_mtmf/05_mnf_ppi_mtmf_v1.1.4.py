# -*- coding: utf-8 -*-
"""
05_mnf_ppi_mtmf.py
ZY1E AHSI 步骤4.5：MNF 去噪降维 + PPI 端元初选 + 包络线去除 + MTMF 填图（v1.1.4）

输入：步骤3 每区裁剪成果（首次运行合并为 VRT）或 mosaic 三件套；可选步骤4 的
      ahsi_library.npy。
流程：①输入校验 ②采样与 MNF ③PPI/n 维端元 ④库投影 ⑤准备输出 ⑥分块 MF/MTMF
      ⑦汇总、对比与质检（各阶段日志打印耗时）。
输出：<OUT_ROOT>/<区>/ 下 MNF 分量、端元表/散点图、MTMF 丰度与不可行性、分类图、
      统计 CSV、参数 JSON、质检图与 SAM/MTMF 对比表。
用法：python 05_mnf_ppi_mtmf_v1.1.4.py --dry-run
依赖：GDAL、numpy、scipy（ConvexHull）、matplotlib（质检图）。

本版仅补充逐行中文备注，算法与输出不变。
"""

import argparse
import csv
import glob
import json
import logging
import os
import sys
import time

import numpy as np


def _bootstrap_gdal_dll():
    if os.name != "nt":                                      # 仅 Windows 需要补 DLL 目录
        return                                               # 其他系统直接返回
    prefixes = []                                            # 收集候选安装前缀
    for key in ("CONDA_PREFIX", "CONDA_DEFAULT_ENV", "VIRTUAL_ENV"):  # 常见环境变量
        value = os.environ.get(key)                          # 读取变量值
        if value and os.path.isdir(value):                    # 非空且目录存在才采用
            prefixes.append(value)                           # 记录该前缀
    prefixes.append(os.path.dirname(sys.executable))         # 追加解释器所在目录
    prefixes += [r"C:\ProgramData\anaconda3", r"D:\ProgramData\anaconda3"]  # 追加兜底路径
    seen = set()                                             # 已注册目录去重
    for prefix in prefixes:                                  # 逐个前缀尝试
        if not prefix or not os.path.isdir(prefix):          # 跳过空值或无效目录
            continue
        candidates = [os.path.join(prefix, "Library", "bin")]  # conda 标准 DLL 目录
        candidates += glob.glob(os.path.join(prefix, "pkgs", "libsqlite-*", "Library", "bin"))  # pkgs 缓存中的 sqlite DLL
        for path in candidates:                              # 逐个候选目录
            if path in seen or not os.path.isdir(path):      # 已处理或不存在则跳过
                continue
            seen.add(path)                                   # 标记已处理
            try:
                os.add_dll_directory(path)                   # 注册为 DLL 搜索目录
            except OSError:
                pass                                         # 注册失败忽略，留给导入报错


_bootstrap_gdal_dll()                                        # 导入 GDAL 前先补 DLL 路径

try:
    from osgeo import gdal, osr                              # 延迟导入，确保先注册 DLL
except ImportError as exc:
    raise SystemExit("缺少 GDAL（osgeo）：%s\n请用装好 GDAL 的 Anaconda 解释器" % exc)  # 明确提示换解释器

gdal.UseExceptions()                                         # GDAL 出错时抛异常
gdal.PushErrorHandler("CPLQuietErrorHandler")                # 屏蔽非致命错误输出


# =================================================================
# 配置区（按需修改）
# =================================================================

IN_ROOT_STEP3 = r"E:\ZhouYang\test\ZY1E_output_step3"       # 步骤3 逐区裁剪成果根目录
STEP4_ROOT = r"E:\ZhouYang\test\ZY1E_output_step4"          # 步骤4 输出根目录（含光谱库）
OUT_ROOT = r"E:\ZhouYang\test\ZY1E_output_step4_mnf"        # 本步默认输出根目录
LIB_NPY = os.path.join(STEP4_ROOT, "_library", "ahsi_library.npy")        # 光谱库矩阵文件
LIB_META = os.path.join(STEP4_ROOT, "_library", "ahsi_library_meta.csv")  # 光谱库元数据

MASK_NAME = "mosaic_mask.tif"                               # mosaic 掩膜文件名
REFL_NAME = "mosaic_reflectance.tif"                        # mosaic 反射率文件名
WL_NAME = "wavelengths_clean.npy"                           # mosaic 波长文件名

# 步骤3 直接产物（逐景裁剪反射率）后缀；旧后缀自动兼容
STEP3_CLIP_SUFFIXES = ("_cal_badband_ortho_atmcorr_clip.tif",   # 新版含坏波段标记
                       "_ortho_atmcorr_clip.tif")               # 旧版后缀
# 波长来源：优先步骤1，其次步骤2（<场景名>/wavelengths_clean.npy）
STEP1_ROOT = r"E:\ZhouYang\test\ZY1E_output_step1"          # 步骤1 输出目录
STEP2_ROOT = r"E:\ZhouYang\test\ZY1E_output_step2"          # 步骤2 输出目录
CACHE_DIR_NAME = "_cache"                                   # VRT 缓存目录名
CACHE_VRT_NAME = "mosaic_reflectance.vrt"                    # 缓存 VRT 文件名
MASK_SAMPLE_BANDS = 12                                      # 建掩膜时抽样波段数

MATCH_MAX_WAVELENGTH_UM = 2.30     # 2.35 μm 以后信号到噪声底，不参与计算
# MNF/PPI/MTMF 只做 2.0~2.30 μm 的 SWIR 子集：可见光-近红外被植被/阴影主导，
# 用全谱会让“不可行性”普遍偏大，淹没 Al-OH 2.2 μm 的微弱特征。
MNF_BAND_MIN_UM = 2.00                                      # SWIR 子集起始波长（μm）
MNF_BAND_MAX_UM = 2.30                                      # SWIR 子集结束波长（μm）
REFL_MIN, REFL_MAX = -0.05, 1.20                            # 有效反射率区间
NODATA_THRESHOLD = -9000.0                                  # 原始无效值哨兵
NODATA_OUT = -9999.0                                        # 输出无效值
BLOCK_ROWS = 96                                             # 分块处理的行数

MNF_VARIANCE_KEEP = 0.99           # MNF 保留方差比例
MNF_MIN_COMPONENTS = 5                                      # MNF 最少分量数
MNF_MAX_COMPONENTS = 20                                     # MNF 最多分量数

STATS_SAMPLE_PIXELS = 150000       # 背景/噪声统计抽样像元
PPI_SKEWERS = 2000                 # PPI 随机投影方向数
PPI_CANDIDATES = 300               # PPI 候选纯像元数
PPI_HULL_DIMS = 5                  # 凸包用的 MNF 维数
N_IMAGE_ENDMEMBERS = 8             # 影像端元个数上限

MTMF_MIN_ABUNDANCE = 0.02           # MTMF 分类的 MF 丰度下限（本景信号弱，先取噪声级阈值）
MTMF_MAX_INFEASIBILITY = 3.0        # MTMF 分类的不可行性上限

MAKE_QC = True                                              # 是否生成质检图
QC_MAX_SIZE = 1400                                          # 质检图最长边像素数
QC_DPI = 140                                                # 质检图保存 DPI

MINERALS = [
    ("Kaolinite", "高岭石", "#5b7db1"), ("Dickite", "地开石", "#e3a857"),  # 高岭石/地开石
    ("Halloysite", "埃洛石", "#8a5a44"), ("Illite", "伊利石", "#2a9d8f"),  # 埃洛石/伊利石
    ("Muscovite", "白云母/绢云母", "#e76f51"), ("Montmorillonite", "蒙脱石", "#9b5de5"),  # 白云母/蒙脱石
    ("Beidellite", "贝得石", "#f4a261"), ("Pyrophyllite", "叶蜡石", "#4d908e"),  # 贝得石/叶蜡石
    ("Alunite", "明矾石", "#f94144"), ("Diaspore", "硬水铝石", "#577590"),  # 明矾石/硬水铝石
    ("Topaz", "黄玉", "#43aa8b"),                                       # 黄玉
]
MINERAL_EN = [m[0] for m in MINERALS]                       # 矿物英文名列表
MINERAL_COLOR = {m[0]: m[2] for m in MINERALS}              # 英文名→颜色映射
CLASS_NODATA = 255                                          # 分类图无数据值

LOG_NAME = "zy1e.step4_5"                                   # 日志器名称


# =================================================================
# 工具
# =================================================================

def parse_args(argv=None):
    """解析命令行参数；默认值全部来自上方配置区。"""
    parser = argparse.ArgumentParser(
        description="ZY1E AHSI 步骤4.5：MNF + PPI/nD 端元 + 包络线去除 + MTMF")     # 创建参数解析器
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP3)          # 步骤3 输入根目录
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT)             # 输出根目录
    parser.add_argument("--area", action="append", default=None)                # 只处理匹配区域，可多次
    parser.add_argument("--overwrite", action="store_true", default=False)      # 覆盖已有成果
    parser.add_argument("--dry-run", action="store_true")                       # 只打印计划不落盘
    parser.add_argument("--log", default=None)                                  # 自定义日志路径
    parser.add_argument("--verbose", action="store_true")                       # 控制台输出 DEBUG
    return parser.parse_args(argv)                                              # 解析并返回参数


def setup_logging(path, verbose=False):
    """初始化日志：控制台 INFO（--verbose 时 DEBUG）+ 文件 DEBUG，重复调用会重建。"""
    logger = logging.getLogger(LOG_NAME)                       # 取得具名日志器
    logger.setLevel(logging.DEBUG)                             # 记录器总级别
    logger.handlers.clear()                                    # 清空旧处理器，防止重复输出
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")  # 统一格式
    console = logging.StreamHandler(sys.stdout)                # 控制台处理器
    console.setLevel(logging.DEBUG if verbose else logging.INFO)  # 按 --verbose 决定级别
    console.setFormatter(fmt)                                  # 套用统一格式
    logger.addHandler(console)                                 # 挂到日志器
    if path:                                                   # 指定了日志文件才写文件
        try:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)  # 确保日志目录存在
            handler = logging.FileHandler(path, mode="a", encoding="utf-8")     # 追加写，保留历史
            handler.setLevel(logging.DEBUG)                    # 文件记录最详细级别
            handler.setFormatter(fmt)                          # 套用统一格式
            logger.addHandler(handler)                         # 挂到日志器
        except OSError:
            pass                                               # 日志文件不可写时静默忽略
    logger.propagate = False                                   # 不向根日志器传播
    return logger                                              # 返回配置好的日志器


def valid_pixels(array):
    """返回“有限且大于 nodata 哨兵(-9000)”的布尔掩膜。"""
    return np.isfinite(array) & (array > -9000.0)              # 有限值且高于哨兵


def read_bands(ds, band_indices, y0, rows):
    """读取指定波段的行块，返回 (波段数, 像元数) 的 float64 数组。"""
    cube = ds.ReadAsArray(0, y0, ds.RasterXSize, rows)         # 一次读入整行块
    if cube is None:
        return None                                            # 读取失败返回 None
    return cube[np.asarray(band_indices)].reshape(len(band_indices), -1).astype(np.float64)  # 选波段并拉平成二维


def _clip_suffix(filename):
    for suffix in STEP3_CLIP_SUFFIXES:                         # 逐个候选后缀比对
        if filename.endswith(suffix):
            return suffix                                      # 命中即返回该后缀
    return None                                                # 都不匹配返回 None


def area_clip_files(in_root, area):
    """列出 <in_root>/<area>/ 下的逐景裁剪反射率（按文件名排序）。"""
    area_dir = os.path.join(in_root, area)                     # 区域目录
    if not os.path.isdir(area_dir):
        return []                                              # 目录不存在返回空表
    return [os.path.join(area_dir, name) for name in sorted(os.listdir(area_dir))
            if _clip_suffix(name) and os.path.isfile(os.path.join(area_dir, name))]  # 只保留裁剪成果


def _scene_wavelengths_path(scene):
    for root in (STEP1_ROOT, STEP2_ROOT):                      # 步骤1 优先，其次步骤2
        path = os.path.join(root, scene, "wavelengths_clean.npy")  # 候选波长文件
        if os.path.isfile(path):
            return path                                        # 找到即返回
    return None                                                # 均不存在返回 None


def load_area_wavelengths(clips, band_count=None):
    """由裁剪文件名反推场景名，读取步骤1/2 的 wavelengths_clean.npy。"""
    for path in clips:                                         # 逐个裁剪文件尝试
        name = os.path.basename(path)                          # 取文件名
        suffix = _clip_suffix(name)                            # 识别裁剪后缀
        scene = name[:-len(suffix)] if suffix else os.path.splitext(name)[0]  # 反推场景名
        wl_path = _scene_wavelengths_path(scene)               # 查找波长文件
        if wl_path is None:
            continue                                           # 该景无波长文件，换下一个
        wavelengths = np.load(wl_path).astype(np.float64)      # 读入波长并转 float64
        if wavelengths.ndim != 1:
            raise ValueError("波长文件不是 1 维: %s" % wl_path)  # 维度不符直接报错
        if band_count is not None and wavelengths.size != band_count:
            raise ValueError("裁剪影像波段数 %d != 波长数 %d（%s）",  # 波段数必须与波长数一致
                             band_count, wavelengths.size, wl_path)
        return wavelengths, scene                             # 返回波长与场景名
    raise FileNotFoundError(
        "未找到裁剪场景对应的 wavelengths_clean.npy（查过步骤1/2 目录）")  # 全部失败则报错


def build_area_vrt(clips, vrt_path):
    """用逐景裁剪成果构建虚拟镶嵌（不落地大栅格）。"""
    tmp_path = vrt_path + ".tmp.vrt"                                     # 先写临时 VRT
    options = gdal.BuildVRTOptions(resolution="highest", VRTNodata=NODATA_OUT)  # 最高分辨率、统一无效值
    vrt = gdal.BuildVRT(tmp_path, clips, options=options)                # 由裁剪成果构建虚拟镶嵌
    if vrt is None:
        raise RuntimeError("构建 VRT 失败: %s" % vrt_path)                # 构建失败报错
    vrt.FlushCache()                                                     # 强制落盘
    vrt = None                                                           # 关闭写句柄
    os.replace(tmp_path, vrt_path)                                       # 原子替换为正式文件


def build_area_mask(vrt_path, mask_path, logger):
    """按抽样波段的有效值与反射率范围生成 0/1 掩膜（1=清晰像元）。"""
    ds = gdal.Open(vrt_path)                                             # 打开镶嵌 VRT
    if ds is None:
        raise RuntimeError("无法打开镶嵌: %s" % vrt_path)                # 打不开直接报错
    width, height, nbands = ds.RasterXSize, ds.RasterYSize, ds.RasterCount  # 记录宽高与波段数
    count = min(MASK_SAMPLE_BANDS, nbands)                               # 实际抽样波段数
    sample = np.unique(np.linspace(0, nbands - 1, count).astype(int))    # 均匀抽样的波段下标
    tmp_path = mask_path + ".tmp.tif"                                    # 临时掩膜路径
    driver = gdal.GetDriverByName("GTiff")                               # GeoTIFF 驱动
    out = driver.Create(tmp_path, width, height, 1, gdal.GDT_Byte,       # 建单波段字节型掩膜
                        options=["TILED=YES", "COMPRESS=DEFLATE", "ZLEVEL=6"])  # 分块 DEFLATE 压缩
    if out is None:
        raise RuntimeError("无法创建掩膜: %s" % mask_path)                # 创建失败报错
    out.SetGeoTransform(ds.GetGeoTransform())                            # 继承仿射变换
    out.SetProjection(ds.GetProjection())                                # 继承投影
    band = out.GetRasterBand(1)                                          # 取掩膜波段
    band.SetNoDataValue(0)                                               # 0 视为无效
    clear_total = 0                                                      # 清晰像元累计
    for y in range(0, height, BLOCK_ROWS):                               # 按块逐行处理
        rows = min(BLOCK_ROWS, height - y)                               # 本块实际行数
        valid = np.ones((rows, width), bool)                             # 初始全部假定有效
        low = np.full((rows, width), np.inf)                             # 各像元最小值
        high = np.full((rows, width), -np.inf)                           # 各像元最大值
        for index in sample:                                             # 逐个抽样波段
            array = ds.GetRasterBand(int(index) + 1).ReadAsArray(0, y, width, rows)  # 读一行块
            array = np.asarray(array, np.float64)                        # 统一为 float64
            ok = valid_pixels(array)                                     # 本波段有效掩膜
            valid &= ok                                                  # 累积所有波段都有效
            low = np.where(ok, np.minimum(low, array), low)              # 更新最小值
            high = np.where(ok, np.maximum(high, array), high)           # 更新最大值
        valid &= (low > REFL_MIN) & (high < REFL_MAX)                    # 反射率须在合理区间
        mask = valid.astype(np.uint8)                                    # 转 0/1 字节掩膜
        band.WriteArray(mask, 0, y)                                      # 写入本块
        clear_total += int(mask.sum())                                   # 累加清晰像元
    out.FlushCache()                                                     # 强制落盘
    out = None                                                           # 关闭写句柄
    ds = None                                                            # 关闭源数据集
    os.replace(tmp_path, mask_path)                                      # 原子替换为正式掩膜
    logger.info("生成掩膜: %s（清晰像元 %.1f%%）",                          # 记录掩膜与清晰率
                mask_path, 100.0 * clear_total / max(1, width * height))


def resolve_area_inputs(in_root, area, cache_root, overwrite, logger, build=True):
    """解析区域输入：优先 mosaic 三件套，否则由步骤3 裁剪成果构建 VRT + 掩膜。"""
    mosaic_dir = os.path.join(in_root, area, "mosaic")            # 区域 mosaic 目录
    refl_path = os.path.join(mosaic_dir, REFL_NAME)               # 反射率文件
    mask_path = os.path.join(mosaic_dir, MASK_NAME)               # 掩膜文件
    wl_path = os.path.join(mosaic_dir, WL_NAME)                   # 波长文件
    if (os.path.isfile(refl_path) and os.path.isfile(mask_path)
            and os.path.isfile(wl_path)):                         # 三件套齐全则直接用
        return {"refl_path": refl_path, "mask_path": mask_path,   # 反射率与掩膜路径
                "wavelengths": np.load(wl_path).astype(np.float64),  # 读入波长
                "source": "mosaic", "scenes": []}                 # 标记来源为 mosaic

    clips = area_clip_files(in_root, area)                        # 否则找逐景裁剪成果
    if not clips:
        raise FileNotFoundError(
            "区域 %s 下既没有 mosaic/ 三件套，也没有步骤3 裁剪成果" % area)  # 两者皆无则报错
    probe = gdal.Open(clips[0])                                   # 打开首个裁剪成果探波段数
    if probe is None:
        raise RuntimeError("无法打开裁剪成果: %s" % clips[0])       # 打不开直接报错
    band_count = probe.RasterCount                                # 记录波段数
    probe = None                                                  # 释放探针数据集
    wavelengths, _ = load_area_wavelengths(clips, band_count)     # 读取对应波长

    cache_dir = os.path.join(cache_root, CACHE_DIR_NAME, area)    # 区域缓存目录
    vrt_path = os.path.join(cache_dir, CACHE_VRT_NAME)            # 缓存 VRT 路径
    cache_mask = os.path.join(cache_dir, MASK_NAME)               # 缓存掩膜路径
    if not build:
        return {"refl_path": vrt_path, "mask_path": cache_mask,   # dry-run 只返回路径
                "wavelengths": wavelengths, "source": "clip", "scenes": clips}

    os.makedirs(cache_dir, exist_ok=True)                         # 确保缓存目录存在
    if overwrite or not os.path.isfile(vrt_path):                 # VRT 缺失或要求覆盖
        build_area_vrt(clips, vrt_path)                           # 构建 VRT
        logger.info("[%s] 合并 %d 景裁剪成果 -> %s", area, len(clips), vrt_path)  # 记录合并
    if overwrite or not os.path.isfile(cache_mask):               # 掩膜缺失或要求覆盖
        build_area_mask(vrt_path, cache_mask, logger)             # 生成掩膜
    return {"refl_path": vrt_path, "mask_path": cache_mask,       # 返回缓存路径
            "wavelengths": wavelengths, "source": "clip", "scenes": clips}


def discover_areas(in_root, patterns=None):
    import fnmatch                                                # 局部导入用于通配匹配
    areas = []                                                    # 收集符合条件的区域
    for name in sorted(os.listdir(in_root)):                      # 排序遍历输入目录
        if not os.path.isdir(os.path.join(in_root, name)):        # 只保留子目录
            continue
        has_mosaic = os.path.isfile(os.path.join(in_root, name, "mosaic", REFL_NAME))  # 是否有 mosaic
        if not has_mosaic and not area_clip_files(in_root, name):
            continue                                              # 两者都没有则跳过
        if patterns and not any(fnmatch.fnmatch(name, p) for p in patterns):
            continue                                              # 不匹配 --area 过滤
        areas.append(name)                                        # 记录该区域
    return areas                                                  # 返回区域列表


# =================================================================
# 采样：背景协方差 / 噪声协方差
# =================================================================

def sample_stats(ds, mask_ds, band_idx, target=STATS_SAMPLE_PIXELS):
    """抽稀行采样有效像元，同时用相邻像元差估计噪声协方差。"""
    width, height = ds.RasterXSize, ds.RasterYSize             # 影像宽高
    step = max(1, height // 200)                               # 最多抽 200 行，计算行间隔
    pixels, diffs = [], []                                     # 有效像元与相邻差样本
    total = 0                                                  # 已采样本数
    for y in range(step // 2, height, step):                   # 按行间隔逐行采样
        block = read_bands(ds, band_idx, y, 1)                 # 读当前行
        if block is None:
            continue                                           # 读取失败跳过
        mask = mask_ds.ReadAsArray(0, y, width, 1)             # 读同一行掩膜
        if mask.ndim == 2:
            mask = mask[0]                                     # 压成一维
        ok = (mask == 1) & valid_pixels(block).all(axis=0)     # 掩膜有效且全部波段有效
        ok &= (block.min(axis=0) > REFL_MIN) & (block.max(axis=0) < REFL_MAX)  # 反射率在区间内
        if not np.any(ok):
            continue                                           # 本行无有效像元
        pixels.append(block[:, ok].T)                          # 收集有效像元光谱
        total += int(ok.sum())                                 # 更新已采数量
        # 相邻像元差：用于 MNF 噪声估计
        pair_ok = ok[:-1] & ok[1:]                             # 相邻两像元都有效
        if np.any(pair_ok):
            diff = block[:, 1:] - block[:, :-1]                # 相邻像元差
            diffs.append(diff[:, pair_ok].T)                   # 收集有效的差值
        if total >= target:
            break                                              # 达到目标数即停止
    if not pixels:
        raise RuntimeError("采样不到有效像元")                   # 一个都没有则报错
    samples = np.vstack(pixels)                                # 纵向堆叠成样本矩阵
    if samples.shape[0] > target:
        stride = max(1, samples.shape[0] // target)            # 超出目标时抽稀
        samples = samples[::stride][:target]
    if diffs:
        noise_diffs = np.vstack(diffs)                         # 堆叠噪声差样本
        if noise_diffs.shape[0] > target:
            stride = max(1, noise_diffs.shape[0] // target)    # 同样抽稀
            noise_diffs = noise_diffs[::stride][:target]
    else:
        noise_diffs = samples[: min(5000, samples.shape[0])] * 0.0  # 无差值时用零占位
    return samples, noise_diffs                                # 返回样本与噪声差


def mnf_transform(samples, noise_diffs, variance_keep=MNF_VARIANCE_KEEP):
    """返回 (mu, M, eigenvalues, k, var_ratio)：T = (X - mu) @ M。"""
    nb = samples.shape[1]                                      # 波段数
    mu = samples.mean(axis=0)                                  # 各波段均值，用于中心化
    cov_noise = 0.5 * np.cov(noise_diffs, rowvar=False)         # 相邻差估噪声协方差（0.5 系数）
    if not np.all(np.isfinite(cov_noise)):
        cov_noise = np.eye(nb) * 1e-6                          # 协方差异常时用极小对角阵兜底
    ridge = 1e-6 * np.trace(cov_noise) / nb + 1e-12            # 岭项，保证可分解
    cov_noise = cov_noise + np.eye(nb) * ridge                 # 加岭项稳定协方差
    try:
        lower = np.linalg.cholesky(cov_noise)                  # Cholesky 白化分解
    except np.linalg.LinAlgError:
        lower = np.eye(nb)                                     # 非正定时退化为单位阵
    lower_inv = np.linalg.inv(lower)                           # 白化矩阵的逆
    whitened = (samples - mu) @ lower_inv.T                    # 去均值并白化
    cov_white = np.cov(whitened, rowvar=False)                 # 白化后的协方差
    eigvals, eigvecs = np.linalg.eigh(cov_white)               # 对称矩阵特征分解
    order = np.argsort(eigvals)[::-1]                          # 特征值从大到小排序
    eigvals = np.maximum(eigvals[order], 0.0)                  # 重排并裁掉负值
    eigvecs = eigvecs[:, order]                                # 特征向量同步重排
    total = float(eigvals.sum()) or 1.0                        # 特征值总和，避免除零
    cumulative = np.cumsum(eigvals) / total                    # 累计方差贡献
    k = int(np.searchsorted(cumulative, variance_keep) + 1)    # 达到目标方差所需分量数
    k = max(MNF_MIN_COMPONENTS, min(MNF_MAX_COMPONENTS, k, eigvecs.shape[1]))  # 限制在上下限内
    matrix = lower_inv.T @ eigvecs[:, :k]                      # 白化逆 × 特征向量 = MNF 变换矩阵
    return mu, matrix, eigvals[:k], k, float(cumulative[k - 1])  # 返回均值/矩阵/特征值/维数/方差比


# =================================================================
# PPI + n 维端元初选
# =================================================================

def ppi_extract(scores, n_skewers=PPI_SKEWERS, rng=None):
    rng = rng or np.random.default_rng(42)                     # 默认固定随机种子，保证可复现
    n, k = scores.shape                                        # 样本数、MNF 维数
    counts = np.zeros(n, np.int64)                             # 每个像元被选为极值的次数
    for _ in range(n_skewers):                                 # 逐个随机投影方向
        direction = rng.standard_normal(k)                     # 高斯随机方向
        norm = np.linalg.norm(direction)                       # 归一化因子
        if norm == 0:
            continue                                           # 零向量跳过
        direction /= norm                                      # 单位化
        projection = scores @ direction                        # 全部样本投影到该方向
        counts[int(np.argmax(projection))] += 1                # 投影最大端计数
        counts[int(np.argmin(projection))] += 1                # 投影最小端计数
    return counts                                              # 返回 PPI 计数


def select_endmembers(scores, counts, candidates=PPI_CANDIDATES, hull_dims=PPI_HULL_DIMS,
                      max_endmembers=N_IMAGE_ENDMEMBERS):
    """PPI 候选 → MNF 前 hull_dims 维凸包顶点 → 端元索引。

    先按 PPI 计数从高到低取 candidates 个候选，再在低维空间做凸包；
    返回的样本下标即影像端元（纯像元），用于后续 MTMF。
    """
    order = np.argsort(counts)[::-1]                           # PPI 计数从高到低排序
    order = order[counts[order] > 0][:candidates]              # 只取计数>0 的前 candidates 个
    if order.size < 2:
        order = np.argsort(counts)[::-1][:max(2, max_endmembers)]  # 候选太少时放宽
    points = scores[order][:, : min(hull_dims, scores.shape[1])]  # 低维空间坐标
    chosen = []                                                # 端元下标列表
    try:
        from scipy.spatial import ConvexHull                   # 局部导入 scipy 凸包
        hull = ConvexHull(points)                              # 计算凸包
        vertices = np.unique(hull.vertices)                    # 凸包顶点（去重）
        chosen = list(order[vertices])                         # 顶点映射回样本下标
        chosen.sort(key=lambda idx: -counts[idx])              # 按 PPI 计数降序
        chosen = chosen[:max_endmembers]                       # 限制端元个数
    except Exception:
        # 回退：最远点贪婪
        chosen = [int(order[0])]                               # 以计数最高者为起点
        while len(chosen) < max_endmembers:                    # 逐个补充端元
            distances = np.min(
                np.linalg.norm(scores[order][:, None, :] - scores[chosen][None, :, :], axis=2),
                axis=1)                                        # 到已选端元的最小距离
            pick = int(order[int(np.argmax(distances))])       # 取最远点
            if pick in chosen:
                break                                          # 已选过则退出
            chosen.append(pick)                                # 加入端元列表
    return np.asarray(chosen, np.int64)                        # 返回端元下标数组


# =================================================================
# 包络线去除（凸包）
# =================================================================

def convex_hull_continuum(wavelength, spectrum):
    """对单条光谱做上凸包（连续统），返回与 spectrum 同长度的连续统。"""
    x = np.asarray(wavelength, np.float64)                    # 波长转 float64
    y = np.asarray(spectrum, np.float64)                      # 反射率转 float64
    order = np.argsort(x)                                     # 按波长升序排序
    xs, ys = x[order], y[order]                               # 排序后的序列
    n = ys.size                                               # 波段数
    if n < 3:
        return np.ones_like(y)                                # 波段太少，连续统取 1
    hull = []                                                 # 上凸包顶点下标栈
    for i in range(n):                                        # 逐点构造凸包
        while len(hull) >= 2:                                 # 至少两个点才能判断
            x1, y1 = xs[hull[-2]], ys[hull[-2]]               # 栈内倒数第二点
            x2, y2 = xs[hull[-1]], ys[hull[-1]]               # 栈顶
            if (y2 - y1) * (xs[i] - x1) >= (ys[i] - y1) * (x2 - x1):  # 叉积判断是否凸
                hull.pop()                                    # 非凸则弹出栈顶
            else:
                break
        hull.append(i)                                        # 当前点入栈
    continuum_sorted = np.full(n, np.nan)                     # 连续统（排序后）
    for a, b in zip(hull[:-1], hull[1:]):                     # 逐条凸包边插值
        seg = slice(a, b + 1)                                 # 该边覆盖的下标
        continuum_sorted[seg] = ys[a] + (ys[b] - ys[a]) * (xs[seg] - xs[a]) / (xs[b] - xs[a])  # 线性插值
    continuum = np.full(n, np.nan)                            # 连续统（原始顺序）
    continuum[order] = continuum_sorted                       # 还原原始波长顺序
    return continuum                                          # 返回连续统


def continuum_depth(wavelength, spectrum, feature_range=(2.12, 2.28)):
    """包络线去除后的吸收深度/位置。"""
    x = np.asarray(wavelength, np.float64)                    # 波长转 float64
    continuum = convex_hull_continuum(x, spectrum)            # 计算上凸包连续统
    with np.errstate(divide="ignore", invalid="ignore"):
        cr = np.asarray(spectrum, np.float64) / continuum     # 包络线去除后的比值
    mask = (x >= feature_range[0]) & (x <= feature_range[1]) & np.isfinite(cr)  # 特征区间内有效点
    if not np.any(mask):
        return np.nan, np.nan                                 # 无有效点返回 NaN
    depth = 1.0 - cr[mask]                                    # 吸收深度
    index = int(np.argmax(depth))                             # 最深处的下标
    return float(depth[index]), float(x[mask][index])         # 返回深度与对应波长


# =================================================================
# 主处理
# =================================================================

def create_gtiff(path, width, height, nbands, dtype, geo, projection, nodata, descs=None):
    """建临时 GeoTIFF（分块 + DEFLATE 压缩），成功后由调用方原子改名。"""
    creation = ["BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "ZLEVEL=6", "INTERLEAVE=BAND"]  # 创建选项
    if dtype != gdal.GDT_Byte:
        creation.append("PREDICTOR=3")                        # 浮点型启用浮点预测器
    driver = gdal.GetDriverByName("GTiff")                    # GeoTIFF 驱动
    ds = driver.Create(path, width, height, nbands, dtype, options=creation)  # 创建数据集
    if ds is None:
        raise RuntimeError("无法创建: %s" % path)              # 创建失败报错
    ds.SetGeoTransform(geo)                                   # 写入仿射变换
    if projection:
        ds.SetProjection(projection)                          # 有投影则写入
    for i in range(1, nbands + 1):                            # 逐波段设置
        band = ds.GetRasterBand(i)                            # 取第 i 波段
        band.SetNoDataValue(nodata)                           # 设置无效值
        if descs and i <= len(descs):
            band.SetDescription(str(descs[i - 1]))            # 写入波段描述
    return ds                                                 # 返回打开的数据集


def load_library():
    """读步骤4 的 AHSI 匹配库；缺失时返回 (None, None) 并以影像端元继续。"""
    if not (os.path.isfile(LIB_NPY) and os.path.isfile(LIB_META)):
        return None, None                                     # 库文件缺失则不使用库端元
    library = np.load(LIB_NPY).astype(np.float64)             # 读光谱库矩阵
    meta = list(csv.DictReader(open(LIB_META, encoding="utf-8-sig")))  # 读元数据表
    class_ids = np.asarray([int(row["class_id"]) for row in meta], np.int32)  # 每行类别号
    return library, class_ids                                 # 返回库矩阵与类别号


def process_area(area, inputs, out_root, args, logger):
    started = time.time()                                        # 区域处理计时起点
    refl_path = inputs["refl_path"]                              # 反射率影像路径
    mask_path = inputs["mask_path"]                              # 掩膜路径
    wl = inputs["wavelengths"]                                   # 波长数组
    out_dir = os.path.join(out_root, area)                       # 区域输出目录
    if not args.overwrite and os.path.isfile(os.path.join(out_dir, "mtmf_class.tif")):
        logger.info("[%s] 跳过（成果已存在）", area)              # 已有成果且不覆盖
        return {"area": area, "status": "跳过", "message": "成果已存在"}  # 直接返回跳过

    source_note = ("mosaic" if inputs["source"] == "mosaic"
                   else "裁剪成果合并 %d 景" % len(inputs["scenes"]))  # 输入来源说明
    logger.info("[%s] 开始处理：%s（%s）", area, refl_path, source_note)  # 记录开始
    logger.info("[%s] 掩膜: %s；波段数 %d", area, mask_path, wl.size)     # 记录掩膜与波段数

    # ========== ① 输入校验：打开影像/掩膜，核对波段数 ==========
    ds = gdal.Open(refl_path)                                    # 打开反射率影像
    mask_ds = gdal.Open(mask_path)                               # 打开掩膜
    if ds is None or mask_ds is None:
        raise RuntimeError("无法打开输入: %s" % refl_path)        # 任一打不开即报错
    width, height = ds.RasterXSize, ds.RasterYSize               # 影像宽高
    geo = ds.GetGeoTransform()                                   # 仿射变换
    projection = ds.GetProjection()                              # 投影字符串
    if ds.RasterCount != wl.size:
        raise ValueError("波段数 %d != 波长数 %d" % (ds.RasterCount, wl.size))  # 波段数须与波长数一致

    band_idx = np.where((wl >= MNF_BAND_MIN_UM) & (wl <= MNF_BAND_MAX_UM))[0]  # 选出 SWIR 子集波段
    logger.info("[%s] %d×%d, SWIR 子集 %d 波段 (%.2f~%.2f μm)", area, width, height,
                band_idx.size, MNF_BAND_MIN_UM, MNF_BAND_MAX_UM)  # 记录子集信息

    # ========== ② 采样 + MNF ==========
    # 只用 2.0~2.30 μm 的 SWIR 子集；相邻像元差近似噪声，先白化再 PCA，
    # 得到变换矩阵 matrix（k × 波段），把原始反射率投影到 MNF 空间。
    t_stage = time.time()                                        # 阶段计时起点
    samples, noise_diffs = sample_stats(ds, mask_ds, band_idx)   # 抽样背景与噪声差
    logger.info("[%s] 采样 %d 像元, 噪声差样本 %d（用时 %.1f s）", area,
                samples.shape[0], noise_diffs.shape[0], time.time() - t_stage)  # 记录采样规模
    t_stage = time.time()                                        # 重置阶段计时
    mu, matrix, eigvals, k, var_ratio = mnf_transform(samples, noise_diffs)  # 计算 MNF 变换
    logger.info("[%s] MNF 保留 %d 分量, 累计方差 %.4f（用时 %.1f s）",
                area, k, var_ratio, time.time() - t_stage)       # 记录 MNF 结果
    scores_sample = (samples - mu) @ matrix                      # 样本投影到 MNF 空间

    # ========== ③ PPI + 影像端元 ==========
    # 在 MNF 空间投随机“针”，两端极值像元作为纯像元候选，再取凸包顶点为端元。
    t_stage = time.time()                                        # 阶段计时起点
    counts = ppi_extract(scores_sample)                          # 计算 PPI 计数
    endmember_idx = select_endmembers(scores_sample, counts)     # 选影像端元
    endmember_scores = scores_sample[endmember_idx]              # 端元的 MNF 坐标
    endmember_spectra = samples[endmember_idx]          # 影像端元（match bands）
    logger.info("[%s] PPI 端元 %d 个（PPI 计数 %s，用时 %.1f s）", area,
                endmember_idx.size, counts[endmember_idx].tolist(),
                time.time() - t_stage)                           # 记录端元数量与计数
    for idx in endmember_idx:                                    # 逐个端元估吸收特征
        depth, pos = continuum_depth(wl[band_idx], samples[idx])  # 包络线去除求深度/位置
        logger.info("    端元 PPI=%4d  λc=%.3f μm  depth=%.3f", counts[idx], pos, depth)  # 打印特征

    # ========== ④ 光谱库端元投影到同一 MNF 空间 ==========
    library, class_ids = load_library()                          # 读光谱库（可能为空）
    lib_scores = None                                            # 库端元 MNF 坐标
    if library is not None:
        lib_match = library[:, band_idx]                         # 只取 SWIR 子集波段
        lib_scores = (lib_match - mu) @ matrix                   # 投到同一 MNF 空间
        logger.info("[%s] 库端元 %d 个投影到 MNF", area, lib_scores.shape[0])  # 记录库端元数
    else:
        logger.warning("[%s] 未找到步骤4 光谱库（%s），跳过库端元 MTMF",
                       area, LIB_NPY)                            # 无库则仅用影像端元

    # ========== ⑤ 准备输出栅格（tmp 文件，全部成功后统一改名） ==========
    os.makedirs(out_dir, exist_ok=True)                          # 确保区域输出目录存在
    created = {}                                                 # 记录待改名的临时成果
    descs = ["%s %s" % (a, b) for a, b, _ in MINERALS]           # 各类矿物波段描述

    def make(name, nbands, dtype, nodata, descriptions=None):
        tmp = os.path.join(out_dir, name + ".tmp.tif")           # 临时文件名
        if os.path.exists(tmp):
            os.remove(tmp)                                       # 先删残留临时文件
        dataset = create_gtiff(tmp, width, height, nbands, dtype, geo, projection,
                               nodata, descriptions)             # 建临时 GeoTIFF
        created[name] = (dataset, tmp, os.path.join(out_dir, name))  # 登记待改名信息
        return dataset                                           # 返回写句柄

    out_mnf = make("mnf_components.tif", k, gdal.GDT_Float32, NODATA_OUT,
                   ["MNF %d" % (i + 1) for i in range(k)])        # MNF 分量输出
    out_abund = None                                             # 库丰度输出占位
    out_infeas = None                                            # 库不可行性输出占位
    if library is not None:
        out_abund = make("mtmf_library_abundance.tif", len(MINERALS), gdal.GDT_Float32,
                         NODATA_OUT, descs)                      # 每类矿物 MF 丰度
        out_infeas = make("mtmf_library_infeasibility.tif", len(MINERALS), gdal.GDT_Float32,
                          NODATA_OUT, descs)                     # 每类矿物不可行性
    out_ppi_abund = make("mtmf_ppi_abundance.tif", endmember_idx.size, gdal.GDT_Float32,
                         NODATA_OUT, ["PPI endmember %d" % (i + 1) for i in range(endmember_idx.size)])  # 影像端元丰度
    out_class = make("mtmf_class.tif", 1, gdal.GDT_Byte, CLASS_NODATA,
                     ["MTMF class"])                             # 分类图

    # MNF 空间背景统计（用全部采样分数）
    mean_scores = scores_sample.mean(axis=0)                     # 背景均值
    cov_scores = np.cov(scores_sample, rowvar=False)             # 背景协方差
    ridge = 1e-8 * np.trace(cov_scores) / k + 1e-12              # 岭项，避免奇异
    inv_cov = np.linalg.pinv(cov_scores + np.eye(k) * ridge, rcond=1e-10)  # 协方差伪逆

    # MF 预计算（匹配滤波的闭式解）：
    #   q     = 背景协方差逆 × (端元 - 背景均值)
    #   denom = (端元 - 背景均值)·q
    # 像元的 MF 分数 = (像元 - 背景均值)·q / denom，可理解为该端元的丰度估计。
    def mf_terms(end_scores):
        q = (end_scores - mean_scores) @ inv_cov            # (n_end, k)
        denom = np.sum((end_scores - mean_scores) * q, axis=1)  # 分母，逐端元
        return q.T, denom                                    # (k, n_end), (n_end,)

    lib_q = lib_denom = None                                     # 库端元 MF 参数占位
    if lib_scores is not None:
        lib_q, lib_denom = mf_terms(lib_scores)                  # 预计算库端元参数
    ppi_q, ppi_denom = mf_terms(endmember_scores)                # 预计算影像端元参数

    class_counts = np.zeros(len(MINERALS), np.int64)             # 各类像元计数
    class_sum_mf = np.zeros(len(MINERALS), np.float64)           # 各类 MF 累加
    valid_total = 0                                              # 有效像元累计
    classified_total = 0                                         # 已分类像元累计

    # ========== ⑥ 分块计算 MF / MTMF（只处理 mask==1 的清晰像元） ==========
    logger.info("[%s] 开始分块处理：%d×%d，块高 %d 行", area, width, height, BLOCK_ROWS)  # 记录分块参数
    total_blocks = (height + BLOCK_ROWS - 1) // BLOCK_ROWS      # 总块数
    for block_index, y0 in enumerate(range(0, height, BLOCK_ROWS), start=1):  # 逐块遍历
        rows = min(BLOCK_ROWS, height - y0)                      # 本块实际行数
        npix = width * rows                                      # 本块像元数
        block = read_bands(ds, band_idx, y0, rows)               # 读本块 SWIR 数据
        if block is None:
            raise RuntimeError("读取失败 y=%d" % y0)              # 读取失败直接报错
        mask = mask_ds.ReadAsArray(0, y0, width, rows)           # 读本块掩膜
        if mask.ndim == 3:
            mask = mask[0]                                       # 去掉波段维
        ok = (mask.reshape(-1) == 1) & valid_pixels(block).all(axis=0)  # 掩膜有效且各波段有效
        ok &= (block.min(axis=0) > REFL_MIN) & (block.max(axis=0) < REFL_MAX)  # 反射率在区间内
        idx = np.where(ok)[0]                                    # 有效像元下标
        if idx.size == 0:
            continue                                             # 本块无有效像元
        scores = (block[:, idx].T - mu) @ matrix           # (nv, k)
        valid_total += int(idx.size)                             # 累加有效像元

        mnf_block = np.full((k, npix), NODATA_OUT, np.float32)   # 整块先填无效值
        mnf_block[:, idx] = scores.T                             # 有效像元写入 MNF 分数
        for b in range(k):
            out_mnf.GetRasterBand(b + 1).WriteArray(
                mnf_block[b].reshape(rows, width), 0, y0)        # 逐分量写盘

        if lib_scores is not None:
            # mf：每类端元的匹配滤波丰度；maha：到背景均值的马氏距离平方；
            # infeas：扣掉端元分量后剩余的马氏距离（越小越“可行”）。
            mf = (scores - mean_scores) @ lib_q / np.maximum(lib_denom[None, :], 1e-12)  # 库端元 MF 丰度
            maha = np.sum(((scores - mean_scores) @ inv_cov) * (scores - mean_scores), axis=1)  # 马氏距离平方
            infeas = np.sqrt(np.maximum(maha[:, None] - (mf ** 2) * lib_denom[None, :], 0.0))  # 剩余马氏距离
            mf_class = np.full((len(MINERALS), idx.size), np.nan)  # 各类丰度缓冲
            inf_class = np.full((len(MINERALS), idx.size), np.nan)  # 各类不可行性缓冲
            for ci in range(len(MINERALS)):                      # 逐类别聚合
                sel = np.where(class_ids == ci + 1)[0]           # 该类对应的库端元
                if sel.size:
                    mf_class[ci] = mf[:, sel].max(axis=1)        # 多元端元取最大丰度
                    inf_class[ci] = infeas[:, sel].min(axis=1)   # 取最小不可行性
            for ci in range(len(MINERALS)):                      # 逐类别写盘
                buf = np.full(npix, NODATA_OUT, np.float32)      # 整块先填无效值
                buf[idx] = np.where(np.isfinite(mf_class[ci]), mf_class[ci], NODATA_OUT)  # 有效像元填丰度
                out_abund.GetRasterBand(ci + 1).WriteArray(buf.reshape(rows, width), 0, y0)  # 写丰度波段
                buf = np.full(npix, NODATA_OUT, np.float32)      # 缓冲区复用
                buf[idx] = np.where(np.isfinite(inf_class[ci]), inf_class[ci], NODATA_OUT)  # 填不可行性
                out_infeas.GetRasterBand(ci + 1).WriteArray(buf.reshape(rows, width), 0, y0)  # 写不可行性波段
            # 分类判决：取 MF 最大的类别，且丰度 ≥ 下限、不可行性 ≤ 上限才归入，
            # 否则记 0（未分类）；无数据像元记 255。
            best = np.argmax(np.where(np.isfinite(mf_class), mf_class, -np.inf), axis=0)  # 每像元最优类别
            best_mf = mf_class[best, np.arange(idx.size)]        # 最优类的丰度
            best_inf = inf_class[best, np.arange(idx.size)]      # 最优类的不可行性
            classify = (best_mf >= MTMF_MIN_ABUNDANCE) & (best_inf <= MTMF_MAX_INFEASIBILITY)  # 双阈值判决
            class_values = np.where(classify, best + 1, 0).astype(np.uint8)  # 类别号 1..N，未分类 0
            buf = np.full(npix, CLASS_NODATA, np.uint8)          # 整块先填无数据
            buf[idx] = class_values                              # 写入分类结果
            out_class.GetRasterBand(1).WriteArray(buf.reshape(rows, width), 0, y0)  # 写分类图
            classified_total += int(np.count_nonzero(classify))  # 累加已分类像元
            class_counts += np.bincount(class_values.astype(np.int64),
                                        minlength=len(MINERALS) + 1)[1:]  # 累加各类计数
            if classify.any():
                for ci in range(len(MINERALS)):                  # 逐类累加 MF 和
                    sel = classify & (class_values == ci + 1)    # 该类像元掩膜
                    if np.any(sel):
                        class_sum_mf[ci] += float(best_mf[sel].sum())  # 累加丰度

        mf_ppi = (scores - mean_scores) @ ppi_q / np.maximum(ppi_denom[None, :], 1e-12)  # 影像端元 MF
        for ei in range(endmember_idx.size):                     # 逐端元写丰度
            buf = np.full(npix, NODATA_OUT, np.float32)          # 整块先填无效值
            buf[idx] = np.clip(mf_ppi[:, ei], 0.0, 1.0)          # 丰度裁剪到 0~1
            out_ppi_abund.GetRasterBand(ei + 1).WriteArray(buf.reshape(rows, width), 0, y0)  # 写盘
        if block_index % 20 == 0 or block_index == total_blocks:  # 每 20 块或最后一块
            logger.info("[%s] 进度 %d/%d 块（%.0f%%）", area, block_index,
                        total_blocks, 100.0 * block_index / total_blocks)  # 打印进度

    # 关闭写盘（清掉局部引用，Windows 下句柄才能释放）
    n_written = len(created)                                     # 成果文件数量
    out_mnf = out_abund = out_infeas = out_ppi_abund = out_class = None  # 释放局部写句柄
    for key in list(created):                                    # 逐个处理临时成果
        dataset, tmp, final = created.pop(key)                   # 取出数据集/临时/正式路径
        dataset.FlushCache()                                     # 强制落盘
        dataset = None                                           # 关闭句柄
        os.replace(tmp, final)                                   # 原子改名为正式文件
    ds = None                                                    # 关闭输入影像
    mask_ds = None                                               # 关闭掩膜
    logger.info("[%s] 写出 %d 个成果文件", area, n_written)        # 记录成果数

    # ========== ⑦ 汇总：统计表 / 元数据 / 与 SAM 对比 / 质检图 ==========
    # 端元 CSV
    with open(os.path.join(out_dir, "mnf_endmembers.csv"), "w", newline="",
              encoding="utf-8") as handle:                       # 端元表
        writer = csv.writer(handle)                              # CSV 写器
        writer.writerow(["endmember_id", "ppi_count"] + ["wl_%.4f" % w for w in wl[band_idx]])  # 表头含波长
        for ei, idx in enumerate(endmember_idx):                 # 逐端元写一行
            writer.writerow([ei + 1, int(counts[idx])] + ["%.5f" % v for v in samples[idx]])  # 编号/计数/光谱

    # 汇总/信息
    rows_summary = []                                            # 各类别汇总行
    for ci, (name_en, name_cn, _) in enumerate(MINERALS):        # 逐类别统计
        count = int(class_counts[ci])                            # 该类像元数
        rows_summary.append({
            "class_id": ci + 1, "mineral_en": name_en, "mineral_cn": name_cn,  # 类别与名称
            "pixels": count,                                     # 像元数
            "area_km2": count * abs(geo[1] * geo[5]) / 1e6,      # 按像元面积折算平方公里
            "mean_mf": (class_sum_mf[ci] / count) if count else 0.0,  # 平均丰度
        })
    with open(os.path.join(out_dir, "mtmf_summary.csv"), "w", newline="",
              encoding="utf-8") as handle:                       # 汇总表
        writer = csv.DictWriter(handle, fieldnames=list(rows_summary[0].keys()))  # 字典写器
        writer.writeheader()                                     # 写表头
        for row in rows_summary:                                 # 逐行格式化
            out = dict(row)                                      # 复制一行
            out["area_km2"] = "%.4f" % row["area_km2"]           # 面积保留 4 位
            out["mean_mf"] = "%.4f" % row["mean_mf"]             # 丰度保留 4 位
            writer.writerow(out)                                 # 写入一行

    info = {
        "area": area,                                            # 区域名
        "bands_used": int(band_idx.size),                        # 参与计算的波段数
        "mnf": {"components": int(k),                            # MNF 分量数
                "eigenvalues": [float(v) for v in eigvals],      # 保留的特征值
                "noise_ridge": float(ridge)},                    # 噪声岭项
        "ppi": {"skewers": PPI_SKEWERS,                          # 投影方向数
                "endmembers": int(endmember_idx.size),           # 端元个数
                "endmember_ppi_counts": [int(counts[i]) for i in endmember_idx]},  # 各端元计数
        "mtmf": {"min_abundance": MTMF_MIN_ABUNDANCE,            # 丰度下限
                 "max_infeasibility": MTMF_MAX_INFEASIBILITY},   # 不可行性上限
        "coverage": {"valid_pixels": int(valid_total),           # 有效像元数
                     "classified_pixels": int(classified_total),  # 已分类像元数
                     "classified_percent": (100.0 * classified_total / valid_total) if valid_total else 0.0},  # 分类占比
        "classes": rows_summary,                                 # 各类别汇总
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),           # 生成时间
    }
    with open(os.path.join(out_dir, "mtmf_info.json"), "w", encoding="utf-8") as handle:
        json.dump(info, handle, ensure_ascii=False, indent=2)    # 写参数 JSON

    # 与 SAM 分类对比
    try:
        compare_sam(out_dir, area, logger)                       # 生成 SAM 对比表
    except Exception as exc:
        logger.warning("[%s] SAM 对比失败（不影响 MTMF 成果）: %s", area, exc,
                       exc_info=True)                            # 对比失败仅告警

    if MAKE_QC:
        try:
            make_qc(out_dir, area, samples, scores_sample, counts, endmember_idx,
                    wl[band_idx], logger)                        # 生成质检图
        except Exception as exc:
            logger.warning("[%s] 质检图失败: %s", area, exc, exc_info=True)  # 绘图失败仅告警

    elapsed = time.time() - started                              # 本区总耗时
    logger.info("[%s] 完成: 有效 %d, MTMF 分类 %d (%.2f%%), 用时 %.1f s",
                area, valid_total, classified_total,
                100.0 * classified_total / max(1, valid_total), elapsed)  # 记录完成统计
    for row in rows_summary:                                     # 逐类别打印
        if row["pixels"]:
            logger.info("    %-18s %8d 像元  %8.2f km2  MF %.3f",
                        row["mineral_en"], row["pixels"], row["area_km2"], row["mean_mf"])  # 打印明细
    return {"area": area, "status": "成功", "message": "分类 %d" % classified_total,
            "elapsed": elapsed}                                  # 返回成功状态


def compare_sam(out_dir, area, logger):
    """读取步骤4 的 aloh_class.tif，与本步 mtmf_class.tif 逐像元对比。

    输出 compare_sam_vs_mtmf.csv（混淆矩阵、总体一致率），用于交叉验证两种方法。
    """
    sam_path = os.path.join(STEP4_ROOT, area, "aloh_class.tif")   # 步骤4 的 SAM 分类
    mtmf_path = os.path.join(out_dir, "mtmf_class.tif")           # 本步的 MTMF 分类
    if not (os.path.isfile(sam_path) and os.path.isfile(mtmf_path)):
        return                                                    # 缺任一文件则跳过
    ds_sam = gdal.Open(sam_path)                                  # 打开 SAM 分类
    ds_mtmf = gdal.Open(mtmf_path)                                # 打开 MTMF 分类
    sam = ds_sam.ReadAsArray()                                    # 读入 SAM 数组
    mtmf = ds_mtmf.ReadAsArray()                                  # 读入 MTMF 数组
    ds_sam = None                                                 # 释放 SAM 数据集
    ds_mtmf = None                                                # 释放 MTMF 数据集
    valid = (sam != CLASS_NODATA) & (mtmf != CLASS_NODATA)        # 两边都有数据
    sam_classified = valid & (sam > 0) & (sam != CLASS_NODATA)    # SAM 已分类像元
    mtmf_classified = valid & (mtmf > 0) & (mtmf != CLASS_NODATA)  # MTMF 已分类像元
    both = sam_classified & mtmf_classified                       # 双方都分类
    agreement = float(np.mean(sam[valid] == mtmf[valid])) if np.any(valid) else 0.0  # 总体一致率
    joint = float(np.mean(sam[both] == mtmf[both])) if np.any(both) else 0.0  # 都分类时一致率
    with open(os.path.join(out_dir, "compare_sam_vs_mtmf.csv"), "w", newline="",
              encoding="utf-8") as handle:                        # 对比结果表
        writer = csv.writer(handle)                               # CSV 写器
        writer.writerow(["metric", "value"])                      # 表头
        writer.writerow(["valid_pixels", int(valid.sum())])       # 有效像元数
        writer.writerow(["overall_agreement", "%.4f" % agreement])  # 总体一致率
        writer.writerow(["both_classified", int(both.sum())])     # 双方都分类数
        writer.writerow(["agreement_when_both_classified", "%.4f" % joint])  # 双方分类时一致率
        writer.writerow(["sam_classified", int(sam_classified.sum())])  # SAM 分类数
        writer.writerow(["mtmf_classified", int(mtmf_classified.sum())])  # MTMF 分类数
    logger.info("[%s] SAM vs MTMF: 总体一致率 %.1f%%，双方都分类时一致率 %.1f%%",
                area, 100 * agreement, 100 * joint)               # 记录一致率


def make_qc(out_dir, area, samples, scores, counts, endmember_idx, wl_band, logger):
    """生成 preview_mtMF.png：MNF 分量、PPI 散点、端元光谱、MTMF 分类与 SAM 对比。"""
    import matplotlib                                      # 延迟导入，避免拖慢主流程
    matplotlib.use("Agg")                                  # 无界面后端
    import matplotlib.pyplot as plt                        # 绘图接口
    from matplotlib.colors import ListedColormap, BoundaryNorm  # 离散色带工具

    def read_preview(name, band=1, dtype=np.float64):
        path = os.path.join(out_dir, name)                 # 预览图路径
        ds = gdal.Open(path)                               # 打开成果栅格
        width, height = ds.RasterXSize, ds.RasterYSize     # 原始宽高
        scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))  # 抽稀倍率
        bw, bh = max(1, width // scale), max(1, height // scale)  # 缩略尺寸
        if ds.RasterCount == 1:
            array = ds.ReadAsArray(buf_xsize=bw, buf_ysize=bh)  # 单波段直接读
        else:
            band_handle = ds.GetRasterBand(band)           # 取指定波段
            array = band_handle.ReadAsArray(buf_xsize=bw, buf_ysize=bh)  # 降采样读取
            band_handle = None                             # 释放波段句柄
        geo = ds.GetGeoTransform()                         # 仿射变换
        extent = [geo[0] / 1000.0, (geo[0] + geo[1] * width) / 1000.0,
                  (geo[3] + geo[5] * height) / 1000.0, geo[3] / 1000.0]  # 转为公里范围
        ds = None                                          # 关闭数据集
        return array, extent                               # 返回数组与范围

    mnf1, extent = read_preview("mnf_components.tif", 1)   # 第 1 MNF 分量
    mnf2, _ = read_preview("mnf_components.tif", 2)        # 第 2 MNF 分量
    class_map, _ = read_preview("mtmf_class.tif")          # 分类图
    colors = ["#dddddd"] + [m[2] for m in MINERALS]        # 背景灰 + 各矿物色
    cmap = ListedColormap(colors)                          # 离散色表
    bounds = np.arange(-0.5, len(colors) + 0.5, 1.0)       # 各类边界

    def stretch_display(array):
        data = np.where(array < -9000, np.nan, array.astype(np.float64))  # 无效值转 NaN
        values = data[np.isfinite(data)]                   # 有效像素
        if values.size:
            low, high = np.percentile(values, [2, 98])     # 2%~98% 拉伸
        else:
            low, high = 0.0, 1.0                           # 无数据时用默认范围
        if high <= low:
            high = low + 1e-6                              # 避免范围为零
        return data, low, high                             # 返回数据与拉伸范围

    figure = plt.figure(figsize=(18, 10))                  # 2×3 画布
    ax1 = figure.add_subplot(2, 3, 1)                      # 第 1 子图
    data1, low1, high1 = stretch_display(mnf1)             # 拉伸 MNF1
    ax1.imshow(data1, extent=extent, origin="upper", cmap="gray", vmin=low1, vmax=high1,
               interpolation="nearest")                    # 显示 MNF1
    ax1.set_title("MNF component 1")                       # 子图标题
    ax2 = figure.add_subplot(2, 3, 2)                      # 第 2 子图
    data2, low2, high2 = stretch_display(mnf2)             # 拉伸 MNF2
    ax2.imshow(data2, extent=extent, origin="upper", cmap="gray", vmin=low2, vmax=high2,
               interpolation="nearest")                    # 显示 MNF2
    ax2.set_title("MNF component 2")                       # 子图标题
    ax3 = figure.add_subplot(2, 3, 3)                      # 第 3 子图
    ax3.imshow(np.where(class_map == CLASS_NODATA, 0, class_map), extent=extent,
               origin="upper", cmap=cmap, norm=BoundaryNorm(bounds, cmap.N),
               interpolation="nearest")                    # 显示分类图
    ax3.set_title("MTMF class")                            # 子图标题
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i]) for i in range(1, len(colors))]  # 图例色块
    ax3.legend(handles, MINERAL_EN, fontsize=6, ncol=2, loc="lower left", framealpha=0.85)  # 矿物图例

    sample_idx = np.random.default_rng(7).choice(scores.shape[0],
                                                 size=min(20000, scores.shape[0]),
                                                 replace=False)  # 抽样散点，最多 2 万
    ax4 = figure.add_subplot(2, 3, 4, projection="3d")     # 三维散点子图
    ax4.scatter(scores[sample_idx, 0], scores[sample_idx, 1], scores[sample_idx, 2],
                s=2, c="0.7", alpha=0.4, label="pixels")    # 样本点
    ax4.scatter(scores[endmember_idx, 0], scores[endmember_idx, 1],
                scores[endmember_idx, 2], s=60, c="red", marker="^", label="endmembers")  # 端元点
    ax4.set_xlabel("MNF1")                                 # X 轴
    ax4.set_ylabel("MNF2")                                 # Y 轴
    ax4.set_zlabel("MNF3")                                 # Z 轴
    ax4.set_title("n-D scatter (MNF1-3) + PPI endmembers")  # 子图标题

    ax5 = figure.add_subplot(2, 3, 5)                      # 端元光谱子图
    for ei, idx in enumerate(endmember_idx):               # 逐端元画曲线
        spectrum = samples[idx]                            # 端元光谱
        depth, pos = continuum_depth(wl_band, spectrum)    # 吸收特征
        ax5.plot(wl_band, spectrum, lw=0.8,
                 label="EM%d (PPI=%d, λc=%.2f)" % (ei + 1, counts[idx], pos))  # 标注编号与特征
    ax5.set_xlabel("Wavelength (um)")                      # X 轴
    ax5.set_ylabel("Reflectance")                          # Y 轴
    ax5.set_title("Image endmember spectra (2.0-2.3 um)")  # 子图标题
    ax5.legend(fontsize=6)                                 # 图例
    ax5.grid(alpha=0.3)                                    # 淡网格

    ax6 = figure.add_subplot(2, 3, 6)                      # 包络线子图
    best_em = endmember_idx[0]                             # 取首个端元示例
    continuum = convex_hull_continuum(wl_band, samples[best_em])  # 计算连续统
    ax6.plot(wl_band, samples[best_em], "k-", lw=1.0, label="EM1 spectrum")  # 原始光谱
    ax6.plot(wl_band, continuum, "r--", lw=1.0, label="convex-hull continuum")  # 连续统
    ax6.set_xlabel("Wavelength (um)")                      # X 轴
    ax6.set_title("Envelope (continuum) removal")          # 子图标题
    ax6.legend(fontsize=7)                                 # 图例
    ax6.grid(alpha=0.3)                                    # 淡网格

    for axis in (ax1, ax2, ax3):                           # 前三个子图统一坐标标注
        axis.set_xlabel("Easting (km)", fontsize=8)        # X 轴
        axis.set_ylabel("Northing (km)", fontsize=8)       # Y 轴
        axis.tick_params(labelsize=7)                      # 刻度字号

    figure.suptitle("%s | MNF + PPI + MTMF" % area, fontsize=13)  # 总标题
    figure.tight_layout(rect=(0, 0, 1, 0.96))              # 自动布局，留出总标题
    path = os.path.join(out_dir, "preview_mtMF.png")       # 质检图路径
    figure.savefig(path, dpi=QC_DPI, bbox_inches="tight")  # 保存图片
    plt.close(figure)                                      # 关闭画布释放内存
    logger.info("[%s] 质检图: %s", area, path)              # 记录输出


def main(argv=None):
    """主流程：解析参数 → 扫描区域 → 逐区处理 → 汇总。

    返回码：0 全部成功/跳过；1 有区域失败；2 输入目录或区域不满足条件。
    """
    args = parse_args(argv)                                      # 解析命令行
    out_root = os.path.abspath(args.out_root)                    # 输出根目录绝对路径
    os.makedirs(out_root, exist_ok=True)                         # 确保输出根目录存在
    log_path = args.log or os.path.join(out_root, "_catalog", "step4_5_batch.log")  # 默认日志路径
    logger = setup_logging(log_path, args.verbose)               # 初始化日志
    gdal.SetCacheMax(1024 * 1024 * 1024)                         # GDAL 缓存上限 1 GB

    logger.info("=" * 72)
    logger.info("ZY1E AHSI 步骤4.5：MNF + PPI/nD + 包络线 + MTMF")
    logger.info("步骤3 目录: %s", os.path.abspath(args.in_root))  # 打印输入目录
    logger.info("输出目录  : %s", out_root)                       # 打印输出目录
    logger.info("=" * 72)

    try:
        areas = discover_areas(os.path.abspath(args.in_root), args.area)  # 扫描待处理区域
    except Exception as exc:
        logger.error("扫描失败: %s", exc)                          # 记录扫描失败
        return 2                                                 # 返回错误码
    if not areas:
        logger.error("没有可处理区域")                            # 无区域可做
        return 2                                                 # 返回错误码

    if args.dry_run:
        for area in areas:                                       # 逐区打印计划
            try:
                inputs = resolve_area_inputs(os.path.abspath(args.in_root), area,
                                             out_root, False, logger, build=False)  # 只解析不建 VRT
            except Exception as exc:
                logger.error("  %-10s 输入不可用: %s", area, exc)  # 输入缺失仅报错
                continue
            if inputs["source"] == "mosaic":
                logger.info("  %-10s <- %s", area, inputs["refl_path"])  # mosaic 直接可用
            else:
                logger.info("  %-10s <- %d 景裁剪成果（将合并为 %s）",
                            area, len(inputs["scenes"]), inputs["refl_path"])  # 需合并
        logger.info("dry-run 结束")                               # dry-run 收尾
        return 0                                                 # 正常退出

    results = []                                                 # 收集每区结果
    batch_started = time.time()                                  # 批处理计时起点
    for index, area in enumerate(areas, start=1):                # 逐区处理
        logger.info("-" * 72)
        logger.info("(%d/%d) 处理: %s", index, len(areas), area)  # 当前进度
        try:
            inputs = resolve_area_inputs(os.path.abspath(args.in_root), area,
                                         out_root, args.overwrite, logger)  # 解析并构建输入
            results.append(process_area(area, inputs, out_root, args, logger))  # 处理该区
        except Exception as exc:
            logger.exception("[%s] 处理失败", area)                # 记录失败堆栈
            results.append({"area": area, "status": "失败", "message": str(exc),
                            "elapsed": 0.0})                     # 记录失败结果

    logger.info("=" * 72)
    success = [r for r in results if r.get("status") == "成功"]    # 成功区域
    skipped = [r for r in results if r.get("status") == "跳过"]    # 跳过区域
    failed = [r for r in results if r.get("status") == "失败"]     # 失败区域
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %.1f s",
                len(success), len(skipped), len(failed), time.time() - batch_started)  # 汇总
    for result in results:                                       # 逐区打印结果
        logger.info("  [%s] %s - %s", result.get("status"), result["area"],
                    result.get("message", ""))                   # 状态与说明
    if failed:
        logger.error("失败区域: %s", ", ".join(r["area"] for r in failed))  # 列出失败区域
        logger.error("排查请查看日志: %s", log_path)               # 提示日志位置
        return 1                                                 # 有失败返回 1
    logger.info("完成，成果: %s", out_root)                       # 全部完成
    logger.info("日志文件: %s", log_path)                         # 日志位置
    return 0                                                     # 成功返回 0


if __name__ == "__main__":
    sys.exit(main())                                             # 以返回码退出进程
