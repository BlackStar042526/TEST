# -*- coding: utf-8 -*-
"""
05_mnf_ppi_mtmf.py
ZY1E AHSI 步骤4.5：MNF 去噪降维 + PPI 端元初选 + 包络线去除 + MTMF 填图（v1.1.3）

输入：步骤3 每区裁剪成果（首次运行合并为 VRT）或 mosaic 三件套；可选步骤4 的
      ahsi_library.npy。
流程：①输入校验 ②采样与 MNF ③PPI/n 维端元 ④库投影 ⑤准备输出 ⑥分块 MF/MTMF
      ⑦汇总、对比与质检（各阶段日志打印耗时）。
输出：<OUT_ROOT>/<区>/ 下 MNF 分量、端元表/散点图、MTMF 丰度与不可行性、分类图、
      统计 CSV、参数 JSON、质检图与 SAM/MTMF 对比表。
用法：python 05_mnf_ppi_mtmf_v1.1.3.py --dry-run
依赖：GDAL、numpy、scipy（ConvexHull）、matplotlib（质检图）。
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
    raise SystemExit("缺少 GDAL（osgeo）：%s\n请用装好 GDAL 的 Anaconda 解释器" % exc)

gdal.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")


# =================================================================
# 配置区（按需修改）
# =================================================================

IN_ROOT_STEP3 = r"E:\ZhouYang\test\ZY1E_output_step3"
STEP4_ROOT = r"E:\ZhouYang\test\ZY1E_output_step4"
OUT_ROOT = r"E:\ZhouYang\test\ZY1E_output_step4_mnf"
LIB_NPY = os.path.join(STEP4_ROOT, "_library", "ahsi_library.npy")
LIB_META = os.path.join(STEP4_ROOT, "_library", "ahsi_library_meta.csv")

MASK_NAME = "mosaic_mask.tif"
REFL_NAME = "mosaic_reflectance.tif"
WL_NAME = "wavelengths_clean.npy"

# 步骤3 直接产物（逐景裁剪反射率）后缀；旧后缀自动兼容
STEP3_CLIP_SUFFIXES = ("_cal_badband_ortho_atmcorr_clip.tif",
                       "_ortho_atmcorr_clip.tif")
# 波长来源：优先步骤1，其次步骤2（<场景名>/wavelengths_clean.npy）
STEP1_ROOT = r"E:\ZhouYang\test\ZY1E_output_step1"
STEP2_ROOT = r"E:\ZhouYang\test\ZY1E_output_step2"
CACHE_DIR_NAME = "_cache"
CACHE_VRT_NAME = "mosaic_reflectance.vrt"
MASK_SAMPLE_BANDS = 12

MATCH_MAX_WAVELENGTH_UM = 2.30     # 2.35 μm 以后信号到噪声底，不参与计算
# MNF/PPI/MTMF 只做 2.0~2.30 μm 的 SWIR 子集：可见光-近红外被植被/阴影主导，
# 用全谱会让“不可行性”普遍偏大，淹没 Al-OH 2.2 μm 的微弱特征。
MNF_BAND_MIN_UM = 2.00
MNF_BAND_MAX_UM = 2.30
REFL_MIN, REFL_MAX = -0.05, 1.20
NODATA_THRESHOLD = -9000.0
NODATA_OUT = -9999.0
BLOCK_ROWS = 96

MNF_VARIANCE_KEEP = 0.99           # MNF 保留方差比例
MNF_MIN_COMPONENTS = 5
MNF_MAX_COMPONENTS = 20

STATS_SAMPLE_PIXELS = 150000       # 背景/噪声统计抽样像元
PPI_SKEWERS = 2000                 # PPI 随机投影方向数
PPI_CANDIDATES = 300               # PPI 候选纯像元数
PPI_HULL_DIMS = 5                  # 凸包用的 MNF 维数
N_IMAGE_ENDMEMBERS = 8             # 影像端元个数上限

MTMF_MIN_ABUNDANCE = 0.02           # MTMF 分类的 MF 丰度下限（本景信号弱，先取噪声级阈值）
MTMF_MAX_INFEASIBILITY = 3.0        # MTMF 分类的不可行性上限

MAKE_QC = True
QC_MAX_SIZE = 1400
QC_DPI = 140

MINERALS = [
    ("Kaolinite", "高岭石", "#5b7db1"), ("Dickite", "地开石", "#e3a857"),
    ("Halloysite", "埃洛石", "#8a5a44"), ("Illite", "伊利石", "#2a9d8f"),
    ("Muscovite", "白云母/绢云母", "#e76f51"), ("Montmorillonite", "蒙脱石", "#9b5de5"),
    ("Beidellite", "贝得石", "#f4a261"), ("Pyrophyllite", "叶蜡石", "#4d908e"),
    ("Alunite", "明矾石", "#f94144"), ("Diaspore", "硬水铝石", "#577590"),
    ("Topaz", "黄玉", "#43aa8b"),
]
MINERAL_EN = [m[0] for m in MINERALS]
MINERAL_COLOR = {m[0]: m[2] for m in MINERALS}
CLASS_NODATA = 255

LOG_NAME = "zy1e.step4_5"


# =================================================================
# 工具
# =================================================================

def parse_args(argv=None):
    """解析命令行参数；默认值全部来自上方配置区。"""
    parser = argparse.ArgumentParser(
        description="ZY1E AHSI 步骤4.5：MNF + PPI/nD 端元 + 包络线去除 + MTMF")
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP3)
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT)
    parser.add_argument("--area", action="append", default=None)
    parser.add_argument("--overwrite", action="store_true", default=False)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log", default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)


def setup_logging(path, verbose=False):
    """初始化日志：控制台 INFO（--verbose 时 DEBUG）+ 文件 DEBUG，重复调用会重建。"""
    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", "%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(fmt)
    logger.addHandler(console)
    if path:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
            handler = logging.FileHandler(path, mode="a", encoding="utf-8")
            handler.setLevel(logging.DEBUG)
            handler.setFormatter(fmt)
            logger.addHandler(handler)
        except OSError:
            pass
    logger.propagate = False
    return logger


def valid_pixels(array):
    """返回“有限且大于 nodata 哨兵(-9000)”的布尔掩膜。"""
    return np.isfinite(array) & (array > -9000.0)


def read_bands(ds, band_indices, y0, rows):
    """读取指定波段的行块，返回 (波段数, 像元数) 的 float64 数组。"""
    cube = ds.ReadAsArray(0, y0, ds.RasterXSize, rows)
    if cube is None:
        return None
    return cube[np.asarray(band_indices)].reshape(len(band_indices), -1).astype(np.float64)


def _clip_suffix(filename):
    for suffix in STEP3_CLIP_SUFFIXES:
        if filename.endswith(suffix):
            return suffix
    return None


def area_clip_files(in_root, area):
    """列出 <in_root>/<area>/ 下的逐景裁剪反射率（按文件名排序）。"""
    area_dir = os.path.join(in_root, area)
    if not os.path.isdir(area_dir):
        return []
    return [os.path.join(area_dir, name) for name in sorted(os.listdir(area_dir))
            if _clip_suffix(name) and os.path.isfile(os.path.join(area_dir, name))]


def _scene_wavelengths_path(scene):
    for root in (STEP1_ROOT, STEP2_ROOT):
        path = os.path.join(root, scene, "wavelengths_clean.npy")
        if os.path.isfile(path):
            return path
    return None


def load_area_wavelengths(clips, band_count=None):
    """由裁剪文件名反推场景名，读取步骤1/2 的 wavelengths_clean.npy。"""
    for path in clips:
        name = os.path.basename(path)
        suffix = _clip_suffix(name)
        scene = name[:-len(suffix)] if suffix else os.path.splitext(name)[0]
        wl_path = _scene_wavelengths_path(scene)
        if wl_path is None:
            continue
        wavelengths = np.load(wl_path).astype(np.float64)
        if wavelengths.ndim != 1:
            raise ValueError("波长文件不是 1 维: %s" % wl_path)
        if band_count is not None and wavelengths.size != band_count:
            raise ValueError("裁剪影像波段数 %d != 波长数 %d（%s）",
                             band_count, wavelengths.size, wl_path)
        return wavelengths, scene
    raise FileNotFoundError(
        "未找到裁剪场景对应的 wavelengths_clean.npy（查过步骤1/2 目录）")


def build_area_vrt(clips, vrt_path):
    """用逐景裁剪成果构建虚拟镶嵌（不落地大栅格）。"""
    tmp_path = vrt_path + ".tmp.vrt"
    options = gdal.BuildVRTOptions(resolution="highest", VRTNodata=NODATA_OUT)
    vrt = gdal.BuildVRT(tmp_path, clips, options=options)
    if vrt is None:
        raise RuntimeError("构建 VRT 失败: %s" % vrt_path)
    vrt.FlushCache()
    vrt = None
    os.replace(tmp_path, vrt_path)


def build_area_mask(vrt_path, mask_path, logger):
    """按抽样波段的有效值与反射率范围生成 0/1 掩膜（1=清晰像元）。"""
    ds = gdal.Open(vrt_path)
    if ds is None:
        raise RuntimeError("无法打开镶嵌: %s" % vrt_path)
    width, height, nbands = ds.RasterXSize, ds.RasterYSize, ds.RasterCount
    count = min(MASK_SAMPLE_BANDS, nbands)
    sample = np.unique(np.linspace(0, nbands - 1, count).astype(int))
    tmp_path = mask_path + ".tmp.tif"
    driver = gdal.GetDriverByName("GTiff")
    out = driver.Create(tmp_path, width, height, 1, gdal.GDT_Byte,
                        options=["TILED=YES", "COMPRESS=DEFLATE", "ZLEVEL=6"])
    if out is None:
        raise RuntimeError("无法创建掩膜: %s" % mask_path)
    out.SetGeoTransform(ds.GetGeoTransform())
    out.SetProjection(ds.GetProjection())
    band = out.GetRasterBand(1)
    band.SetNoDataValue(0)
    clear_total = 0
    for y in range(0, height, BLOCK_ROWS):
        rows = min(BLOCK_ROWS, height - y)
        valid = np.ones((rows, width), bool)
        low = np.full((rows, width), np.inf)
        high = np.full((rows, width), -np.inf)
        for index in sample:
            array = ds.GetRasterBand(int(index) + 1).ReadAsArray(0, y, width, rows)
            array = np.asarray(array, np.float64)
            ok = valid_pixels(array)
            valid &= ok
            low = np.where(ok, np.minimum(low, array), low)
            high = np.where(ok, np.maximum(high, array), high)
        valid &= (low > REFL_MIN) & (high < REFL_MAX)
        mask = valid.astype(np.uint8)
        band.WriteArray(mask, 0, y)
        clear_total += int(mask.sum())
    out.FlushCache()
    out = None
    ds = None
    os.replace(tmp_path, mask_path)
    logger.info("生成掩膜: %s（清晰像元 %.1f%%）",
                mask_path, 100.0 * clear_total / max(1, width * height))


def resolve_area_inputs(in_root, area, cache_root, overwrite, logger, build=True):
    """解析区域输入：优先 mosaic 三件套，否则由步骤3 裁剪成果构建 VRT + 掩膜。"""
    mosaic_dir = os.path.join(in_root, area, "mosaic")
    refl_path = os.path.join(mosaic_dir, REFL_NAME)
    mask_path = os.path.join(mosaic_dir, MASK_NAME)
    wl_path = os.path.join(mosaic_dir, WL_NAME)
    if (os.path.isfile(refl_path) and os.path.isfile(mask_path)
            and os.path.isfile(wl_path)):
        return {"refl_path": refl_path, "mask_path": mask_path,
                "wavelengths": np.load(wl_path).astype(np.float64),
                "source": "mosaic", "scenes": []}

    clips = area_clip_files(in_root, area)
    if not clips:
        raise FileNotFoundError(
            "区域 %s 下既没有 mosaic/ 三件套，也没有步骤3 裁剪成果" % area)
    probe = gdal.Open(clips[0])
    if probe is None:
        raise RuntimeError("无法打开裁剪成果: %s" % clips[0])
    band_count = probe.RasterCount
    probe = None
    wavelengths, _ = load_area_wavelengths(clips, band_count)

    cache_dir = os.path.join(cache_root, CACHE_DIR_NAME, area)
    vrt_path = os.path.join(cache_dir, CACHE_VRT_NAME)
    cache_mask = os.path.join(cache_dir, MASK_NAME)
    if not build:
        return {"refl_path": vrt_path, "mask_path": cache_mask,
                "wavelengths": wavelengths, "source": "clip", "scenes": clips}

    os.makedirs(cache_dir, exist_ok=True)
    if overwrite or not os.path.isfile(vrt_path):
        build_area_vrt(clips, vrt_path)
        logger.info("[%s] 合并 %d 景裁剪成果 -> %s", area, len(clips), vrt_path)
    if overwrite or not os.path.isfile(cache_mask):
        build_area_mask(vrt_path, cache_mask, logger)
    return {"refl_path": vrt_path, "mask_path": cache_mask,
            "wavelengths": wavelengths, "source": "clip", "scenes": clips}


def discover_areas(in_root, patterns=None):
    import fnmatch
    areas = []
    for name in sorted(os.listdir(in_root)):
        if not os.path.isdir(os.path.join(in_root, name)):
            continue
        has_mosaic = os.path.isfile(os.path.join(in_root, name, "mosaic", REFL_NAME))
        if not has_mosaic and not area_clip_files(in_root, name):
            continue
        if patterns and not any(fnmatch.fnmatch(name, p) for p in patterns):
            continue
        areas.append(name)
    return areas


# =================================================================
# 采样：背景协方差 / 噪声协方差
# =================================================================

def sample_stats(ds, mask_ds, band_idx, target=STATS_SAMPLE_PIXELS):
    """抽稀行采样有效像元，同时用相邻像元差估计噪声协方差。"""
    width, height = ds.RasterXSize, ds.RasterYSize
    step = max(1, height // 200)
    pixels, diffs = [], []
    total = 0
    for y in range(step // 2, height, step):
        block = read_bands(ds, band_idx, y, 1)
        if block is None:
            continue
        mask = mask_ds.ReadAsArray(0, y, width, 1)
        if mask.ndim == 2:
            mask = mask[0]
        ok = (mask == 1) & valid_pixels(block).all(axis=0)
        ok &= (block.min(axis=0) > REFL_MIN) & (block.max(axis=0) < REFL_MAX)
        if not np.any(ok):
            continue
        pixels.append(block[:, ok].T)
        total += int(ok.sum())
        # 相邻像元差：用于 MNF 噪声估计
        pair_ok = ok[:-1] & ok[1:]
        if np.any(pair_ok):
            diff = block[:, 1:] - block[:, :-1]
            diffs.append(diff[:, pair_ok].T)
        if total >= target:
            break
    if not pixels:
        raise RuntimeError("采样不到有效像元")
    samples = np.vstack(pixels)
    if samples.shape[0] > target:
        stride = max(1, samples.shape[0] // target)
        samples = samples[::stride][:target]
    if diffs:
        noise_diffs = np.vstack(diffs)
        if noise_diffs.shape[0] > target:
            stride = max(1, noise_diffs.shape[0] // target)
            noise_diffs = noise_diffs[::stride][:target]
    else:
        noise_diffs = samples[: min(5000, samples.shape[0])] * 0.0
    return samples, noise_diffs


def mnf_transform(samples, noise_diffs, variance_keep=MNF_VARIANCE_KEEP):
    """返回 (mu, M, eigenvalues, k, var_ratio)：T = (X - mu) @ M。"""
    nb = samples.shape[1]
    mu = samples.mean(axis=0)
    cov_noise = 0.5 * np.cov(noise_diffs, rowvar=False)
    if not np.all(np.isfinite(cov_noise)):
        cov_noise = np.eye(nb) * 1e-6
    ridge = 1e-6 * np.trace(cov_noise) / nb + 1e-12
    cov_noise = cov_noise + np.eye(nb) * ridge
    try:
        lower = np.linalg.cholesky(cov_noise)
    except np.linalg.LinAlgError:
        lower = np.eye(nb)
    lower_inv = np.linalg.inv(lower)
    whitened = (samples - mu) @ lower_inv.T
    cov_white = np.cov(whitened, rowvar=False)
    eigvals, eigvecs = np.linalg.eigh(cov_white)
    order = np.argsort(eigvals)[::-1]
    eigvals = np.maximum(eigvals[order], 0.0)
    eigvecs = eigvecs[:, order]
    total = float(eigvals.sum()) or 1.0
    cumulative = np.cumsum(eigvals) / total
    k = int(np.searchsorted(cumulative, variance_keep) + 1)
    k = max(MNF_MIN_COMPONENTS, min(MNF_MAX_COMPONENTS, k, eigvecs.shape[1]))
    matrix = lower_inv.T @ eigvecs[:, :k]
    return mu, matrix, eigvals[:k], k, float(cumulative[k - 1])


# =================================================================
# PPI + n 维端元初选
# =================================================================

def ppi_extract(scores, n_skewers=PPI_SKEWERS, rng=None):
    rng = rng or np.random.default_rng(42)
    n, k = scores.shape
    counts = np.zeros(n, np.int64)
    for _ in range(n_skewers):
        direction = rng.standard_normal(k)
        norm = np.linalg.norm(direction)
        if norm == 0:
            continue
        direction /= norm
        projection = scores @ direction
        counts[int(np.argmax(projection))] += 1
        counts[int(np.argmin(projection))] += 1
    return counts


def select_endmembers(scores, counts, candidates=PPI_CANDIDATES, hull_dims=PPI_HULL_DIMS,
                      max_endmembers=N_IMAGE_ENDMEMBERS):
    """PPI 候选 → MNF 前 hull_dims 维凸包顶点 → 端元索引。

    先按 PPI 计数从高到低取 candidates 个候选，再在低维空间做凸包；
    返回的样本下标即影像端元（纯像元），用于后续 MTMF。
    """
    order = np.argsort(counts)[::-1]
    order = order[counts[order] > 0][:candidates]
    if order.size < 2:
        order = np.argsort(counts)[::-1][:max(2, max_endmembers)]
    points = scores[order][:, : min(hull_dims, scores.shape[1])]
    chosen = []
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(points)
        vertices = np.unique(hull.vertices)
        chosen = list(order[vertices])
        chosen.sort(key=lambda idx: -counts[idx])
        chosen = chosen[:max_endmembers]
    except Exception:
        # 回退：最远点贪婪
        chosen = [int(order[0])]
        while len(chosen) < max_endmembers:
            distances = np.min(
                np.linalg.norm(scores[order][:, None, :] - scores[chosen][None, :, :], axis=2),
                axis=1)
            pick = int(order[int(np.argmax(distances))])
            if pick in chosen:
                break
            chosen.append(pick)
    return np.asarray(chosen, np.int64)


# =================================================================
# 包络线去除（凸包）
# =================================================================

def convex_hull_continuum(wavelength, spectrum):
    """对单条光谱做上凸包（连续统），返回与 spectrum 同长度的连续统。"""
    x = np.asarray(wavelength, np.float64)
    y = np.asarray(spectrum, np.float64)
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    n = ys.size
    if n < 3:
        return np.ones_like(y)
    hull = []
    for i in range(n):
        while len(hull) >= 2:
            x1, y1 = xs[hull[-2]], ys[hull[-2]]
            x2, y2 = xs[hull[-1]], ys[hull[-1]]
            if (y2 - y1) * (xs[i] - x1) >= (ys[i] - y1) * (x2 - x1):
                hull.pop()
            else:
                break
        hull.append(i)
    continuum_sorted = np.full(n, np.nan)
    for a, b in zip(hull[:-1], hull[1:]):
        seg = slice(a, b + 1)
        continuum_sorted[seg] = ys[a] + (ys[b] - ys[a]) * (xs[seg] - xs[a]) / (xs[b] - xs[a])
    continuum = np.full(n, np.nan)
    continuum[order] = continuum_sorted
    return continuum


def continuum_depth(wavelength, spectrum, feature_range=(2.12, 2.28)):
    """包络线去除后的吸收深度/位置。"""
    x = np.asarray(wavelength, np.float64)
    continuum = convex_hull_continuum(x, spectrum)
    with np.errstate(divide="ignore", invalid="ignore"):
        cr = np.asarray(spectrum, np.float64) / continuum
    mask = (x >= feature_range[0]) & (x <= feature_range[1]) & np.isfinite(cr)
    if not np.any(mask):
        return np.nan, np.nan
    depth = 1.0 - cr[mask]
    index = int(np.argmax(depth))
    return float(depth[index]), float(x[mask][index])


# =================================================================
# 主处理
# =================================================================

def create_gtiff(path, width, height, nbands, dtype, geo, projection, nodata, descs=None):
    """建临时 GeoTIFF（分块 + DEFLATE 压缩），成功后由调用方原子改名。"""
    creation = ["BIGTIFF=IF_SAFER", "COMPRESS=DEFLATE", "ZLEVEL=6", "INTERLEAVE=BAND"]
    if dtype != gdal.GDT_Byte:
        creation.append("PREDICTOR=3")
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(path, width, height, nbands, dtype, options=creation)
    if ds is None:
        raise RuntimeError("无法创建: %s" % path)
    ds.SetGeoTransform(geo)
    if projection:
        ds.SetProjection(projection)
    for i in range(1, nbands + 1):
        band = ds.GetRasterBand(i)
        band.SetNoDataValue(nodata)
        if descs and i <= len(descs):
            band.SetDescription(str(descs[i - 1]))
    return ds


def load_library():
    """读步骤4 的 AHSI 匹配库；缺失时返回 (None, None) 并以影像端元继续。"""
    if not (os.path.isfile(LIB_NPY) and os.path.isfile(LIB_META)):
        return None, None
    library = np.load(LIB_NPY).astype(np.float64)
    meta = list(csv.DictReader(open(LIB_META, encoding="utf-8-sig")))
    class_ids = np.asarray([int(row["class_id"]) for row in meta], np.int32)
    return library, class_ids


def process_area(area, inputs, out_root, args, logger):
    started = time.time()
    refl_path = inputs["refl_path"]
    mask_path = inputs["mask_path"]
    wl = inputs["wavelengths"]
    out_dir = os.path.join(out_root, area)
    if not args.overwrite and os.path.isfile(os.path.join(out_dir, "mtmf_class.tif")):
        logger.info("[%s] 跳过（成果已存在）", area)
        return {"area": area, "status": "跳过", "message": "成果已存在"}

    source_note = ("mosaic" if inputs["source"] == "mosaic"
                   else "裁剪成果合并 %d 景" % len(inputs["scenes"]))
    logger.info("[%s] 开始处理：%s（%s）", area, refl_path, source_note)
    logger.info("[%s] 掩膜: %s；波段数 %d", area, mask_path, wl.size)

    # ========== ① 输入校验：打开影像/掩膜，核对波段数 ==========
    ds = gdal.Open(refl_path)
    mask_ds = gdal.Open(mask_path)
    if ds is None or mask_ds is None:
        raise RuntimeError("无法打开输入: %s" % refl_path)
    width, height = ds.RasterXSize, ds.RasterYSize
    geo = ds.GetGeoTransform()
    projection = ds.GetProjection()
    if ds.RasterCount != wl.size:
        raise ValueError("波段数 %d != 波长数 %d" % (ds.RasterCount, wl.size))

    band_idx = np.where((wl >= MNF_BAND_MIN_UM) & (wl <= MNF_BAND_MAX_UM))[0]
    logger.info("[%s] %d×%d, SWIR 子集 %d 波段 (%.2f~%.2f μm)", area, width, height,
                band_idx.size, MNF_BAND_MIN_UM, MNF_BAND_MAX_UM)

    # ========== ② 采样 + MNF ==========
    # 只用 2.0~2.30 μm 的 SWIR 子集；相邻像元差近似噪声，先白化再 PCA，
    # 得到变换矩阵 matrix（k × 波段），把原始反射率投影到 MNF 空间。
    t_stage = time.time()
    samples, noise_diffs = sample_stats(ds, mask_ds, band_idx)
    logger.info("[%s] 采样 %d 像元, 噪声差样本 %d（用时 %.1f s）", area,
                samples.shape[0], noise_diffs.shape[0], time.time() - t_stage)
    t_stage = time.time()
    mu, matrix, eigvals, k, var_ratio = mnf_transform(samples, noise_diffs)
    logger.info("[%s] MNF 保留 %d 分量, 累计方差 %.4f（用时 %.1f s）",
                area, k, var_ratio, time.time() - t_stage)
    scores_sample = (samples - mu) @ matrix

    # ========== ③ PPI + 影像端元 ==========
    # 在 MNF 空间投随机“针”，两端极值像元作为纯像元候选，再取凸包顶点为端元。
    t_stage = time.time()
    counts = ppi_extract(scores_sample)
    endmember_idx = select_endmembers(scores_sample, counts)
    endmember_scores = scores_sample[endmember_idx]
    endmember_spectra = samples[endmember_idx]          # 影像端元（match bands）
    logger.info("[%s] PPI 端元 %d 个（PPI 计数 %s，用时 %.1f s）", area,
                endmember_idx.size, counts[endmember_idx].tolist(),
                time.time() - t_stage)
    for idx in endmember_idx:
        depth, pos = continuum_depth(wl[band_idx], samples[idx])
        logger.info("    端元 PPI=%4d  λc=%.3f μm  depth=%.3f", counts[idx], pos, depth)

    # ========== ④ 光谱库端元投影到同一 MNF 空间 ==========
    library, class_ids = load_library()
    lib_scores = None
    if library is not None:
        lib_match = library[:, band_idx]
        lib_scores = (lib_match - mu) @ matrix
        logger.info("[%s] 库端元 %d 个投影到 MNF", area, lib_scores.shape[0])
    else:
        logger.warning("[%s] 未找到步骤4 光谱库（%s），跳过库端元 MTMF",
                       area, LIB_NPY)

    # ========== ⑤ 准备输出栅格（tmp 文件，全部成功后统一改名） ==========
    os.makedirs(out_dir, exist_ok=True)
    created = {}
    descs = ["%s %s" % (a, b) for a, b, _ in MINERALS]

    def make(name, nbands, dtype, nodata, descriptions=None):
        tmp = os.path.join(out_dir, name + ".tmp.tif")
        if os.path.exists(tmp):
            os.remove(tmp)
        dataset = create_gtiff(tmp, width, height, nbands, dtype, geo, projection,
                               nodata, descriptions)
        created[name] = (dataset, tmp, os.path.join(out_dir, name))
        return dataset

    out_mnf = make("mnf_components.tif", k, gdal.GDT_Float32, NODATA_OUT,
                   ["MNF %d" % (i + 1) for i in range(k)])
    out_abund = None
    out_infeas = None
    if library is not None:
        out_abund = make("mtmf_library_abundance.tif", len(MINERALS), gdal.GDT_Float32,
                         NODATA_OUT, descs)
        out_infeas = make("mtmf_library_infeasibility.tif", len(MINERALS), gdal.GDT_Float32,
                          NODATA_OUT, descs)
    out_ppi_abund = make("mtmf_ppi_abundance.tif", endmember_idx.size, gdal.GDT_Float32,
                         NODATA_OUT, ["PPI endmember %d" % (i + 1) for i in range(endmember_idx.size)])
    out_class = make("mtmf_class.tif", 1, gdal.GDT_Byte, CLASS_NODATA,
                     ["MTMF class"])

    # MNF 空间背景统计（用全部采样分数）
    mean_scores = scores_sample.mean(axis=0)
    cov_scores = np.cov(scores_sample, rowvar=False)
    ridge = 1e-8 * np.trace(cov_scores) / k + 1e-12
    inv_cov = np.linalg.pinv(cov_scores + np.eye(k) * ridge, rcond=1e-10)

    # MF 预计算（匹配滤波的闭式解）：
    #   q     = 背景协方差逆 × (端元 - 背景均值)
    #   denom = (端元 - 背景均值)·q
    # 像元的 MF 分数 = (像元 - 背景均值)·q / denom，可理解为该端元的丰度估计。
    def mf_terms(end_scores):
        q = (end_scores - mean_scores) @ inv_cov            # (n_end, k)
        denom = np.sum((end_scores - mean_scores) * q, axis=1)
        return q.T, denom                                    # (k, n_end), (n_end,)

    lib_q = lib_denom = None
    if lib_scores is not None:
        lib_q, lib_denom = mf_terms(lib_scores)
    ppi_q, ppi_denom = mf_terms(endmember_scores)

    class_counts = np.zeros(len(MINERALS), np.int64)
    class_sum_mf = np.zeros(len(MINERALS), np.float64)
    valid_total = 0
    classified_total = 0

    # ========== ⑥ 分块计算 MF / MTMF（只处理 mask==1 的清晰像元） ==========
    logger.info("[%s] 开始分块处理：%d×%d，块高 %d 行", area, width, height, BLOCK_ROWS)
    total_blocks = (height + BLOCK_ROWS - 1) // BLOCK_ROWS
    for block_index, y0 in enumerate(range(0, height, BLOCK_ROWS), start=1):
        rows = min(BLOCK_ROWS, height - y0)
        npix = width * rows
        block = read_bands(ds, band_idx, y0, rows)
        if block is None:
            raise RuntimeError("读取失败 y=%d" % y0)
        mask = mask_ds.ReadAsArray(0, y0, width, rows)
        if mask.ndim == 3:
            mask = mask[0]
        ok = (mask.reshape(-1) == 1) & valid_pixels(block).all(axis=0)
        ok &= (block.min(axis=0) > REFL_MIN) & (block.max(axis=0) < REFL_MAX)
        idx = np.where(ok)[0]
        if idx.size == 0:
            continue
        scores = (block[:, idx].T - mu) @ matrix           # (nv, k)
        valid_total += int(idx.size)

        mnf_block = np.full((k, npix), NODATA_OUT, np.float32)
        mnf_block[:, idx] = scores.T
        for b in range(k):
            out_mnf.GetRasterBand(b + 1).WriteArray(
                mnf_block[b].reshape(rows, width), 0, y0)

        if lib_scores is not None:
            # mf：每类端元的匹配滤波丰度；maha：到背景均值的马氏距离平方；
            # infeas：扣掉端元分量后剩余的马氏距离（越小越“可行”）。
            mf = (scores - mean_scores) @ lib_q / np.maximum(lib_denom[None, :], 1e-12)
            maha = np.sum(((scores - mean_scores) @ inv_cov) * (scores - mean_scores), axis=1)
            infeas = np.sqrt(np.maximum(maha[:, None] - (mf ** 2) * lib_denom[None, :], 0.0))
            mf_class = np.full((len(MINERALS), idx.size), np.nan)
            inf_class = np.full((len(MINERALS), idx.size), np.nan)
            for ci in range(len(MINERALS)):
                sel = np.where(class_ids == ci + 1)[0]
                if sel.size:
                    mf_class[ci] = mf[:, sel].max(axis=1)
                    inf_class[ci] = infeas[:, sel].min(axis=1)
            for ci in range(len(MINERALS)):
                buf = np.full(npix, NODATA_OUT, np.float32)
                buf[idx] = np.where(np.isfinite(mf_class[ci]), mf_class[ci], NODATA_OUT)
                out_abund.GetRasterBand(ci + 1).WriteArray(buf.reshape(rows, width), 0, y0)
                buf = np.full(npix, NODATA_OUT, np.float32)
                buf[idx] = np.where(np.isfinite(inf_class[ci]), inf_class[ci], NODATA_OUT)
                out_infeas.GetRasterBand(ci + 1).WriteArray(buf.reshape(rows, width), 0, y0)
            # 分类判决：取 MF 最大的类别，且丰度 ≥ 下限、不可行性 ≤ 上限才归入，
            # 否则记 0（未分类）；无数据像元记 255。
            best = np.argmax(np.where(np.isfinite(mf_class), mf_class, -np.inf), axis=0)
            best_mf = mf_class[best, np.arange(idx.size)]
            best_inf = inf_class[best, np.arange(idx.size)]
            classify = (best_mf >= MTMF_MIN_ABUNDANCE) & (best_inf <= MTMF_MAX_INFEASIBILITY)
            class_values = np.where(classify, best + 1, 0).astype(np.uint8)
            buf = np.full(npix, CLASS_NODATA, np.uint8)
            buf[idx] = class_values
            out_class.GetRasterBand(1).WriteArray(buf.reshape(rows, width), 0, y0)
            classified_total += int(np.count_nonzero(classify))
            class_counts += np.bincount(class_values.astype(np.int64),
                                        minlength=len(MINERALS) + 1)[1:]
            if classify.any():
                for ci in range(len(MINERALS)):
                    sel = classify & (class_values == ci + 1)
                    if np.any(sel):
                        class_sum_mf[ci] += float(best_mf[sel].sum())

        mf_ppi = (scores - mean_scores) @ ppi_q / np.maximum(ppi_denom[None, :], 1e-12)
        for ei in range(endmember_idx.size):
            buf = np.full(npix, NODATA_OUT, np.float32)
            buf[idx] = np.clip(mf_ppi[:, ei], 0.0, 1.0)
            out_ppi_abund.GetRasterBand(ei + 1).WriteArray(buf.reshape(rows, width), 0, y0)
        if block_index % 20 == 0 or block_index == total_blocks:
            logger.info("[%s] 进度 %d/%d 块（%.0f%%）", area, block_index,
                        total_blocks, 100.0 * block_index / total_blocks)

    # 关闭写盘（清掉局部引用，Windows 下句柄才能释放）
    n_written = len(created)
    out_mnf = out_abund = out_infeas = out_ppi_abund = out_class = None
    for key in list(created):
        dataset, tmp, final = created.pop(key)
        dataset.FlushCache()
        dataset = None
        os.replace(tmp, final)
    ds = None
    mask_ds = None
    logger.info("[%s] 写出 %d 个成果文件", area, n_written)

    # ========== ⑦ 汇总：统计表 / 元数据 / 与 SAM 对比 / 质检图 ==========
    # 端元 CSV
    with open(os.path.join(out_dir, "mnf_endmembers.csv"), "w", newline="",
              encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["endmember_id", "ppi_count"] + ["wl_%.4f" % w for w in wl[band_idx]])
        for ei, idx in enumerate(endmember_idx):
            writer.writerow([ei + 1, int(counts[idx])] + ["%.5f" % v for v in samples[idx]])

    # 汇总/信息
    rows_summary = []
    for ci, (name_en, name_cn, _) in enumerate(MINERALS):
        count = int(class_counts[ci])
        rows_summary.append({
            "class_id": ci + 1, "mineral_en": name_en, "mineral_cn": name_cn,
            "pixels": count,
            "area_km2": count * abs(geo[1] * geo[5]) / 1e6,
            "mean_mf": (class_sum_mf[ci] / count) if count else 0.0,
        })
    with open(os.path.join(out_dir, "mtmf_summary.csv"), "w", newline="",
              encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows_summary[0].keys()))
        writer.writeheader()
        for row in rows_summary:
            out = dict(row)
            out["area_km2"] = "%.4f" % row["area_km2"]
            out["mean_mf"] = "%.4f" % row["mean_mf"]
            writer.writerow(out)

    info = {
        "area": area,
        "bands_used": int(band_idx.size),
        "mnf": {"components": int(k),
                "eigenvalues": [float(v) for v in eigvals],
                "noise_ridge": float(ridge)},
        "ppi": {"skewers": PPI_SKEWERS,
                "endmembers": int(endmember_idx.size),
                "endmember_ppi_counts": [int(counts[i]) for i in endmember_idx]},
        "mtmf": {"min_abundance": MTMF_MIN_ABUNDANCE,
                 "max_infeasibility": MTMF_MAX_INFEASIBILITY},
        "coverage": {"valid_pixels": int(valid_total),
                     "classified_pixels": int(classified_total),
                     "classified_percent": (100.0 * classified_total / valid_total) if valid_total else 0.0},
        "classes": rows_summary,
        "created": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    with open(os.path.join(out_dir, "mtmf_info.json"), "w", encoding="utf-8") as handle:
        json.dump(info, handle, ensure_ascii=False, indent=2)

    # 与 SAM 分类对比
    try:
        compare_sam(out_dir, area, logger)
    except Exception as exc:
        logger.warning("[%s] SAM 对比失败（不影响 MTMF 成果）: %s", area, exc,
                       exc_info=True)

    if MAKE_QC:
        try:
            make_qc(out_dir, area, samples, scores_sample, counts, endmember_idx,
                    wl[band_idx], logger)
        except Exception as exc:
            logger.warning("[%s] 质检图失败: %s", area, exc, exc_info=True)

    elapsed = time.time() - started
    logger.info("[%s] 完成: 有效 %d, MTMF 分类 %d (%.2f%%), 用时 %.1f s",
                area, valid_total, classified_total,
                100.0 * classified_total / max(1, valid_total), elapsed)
    for row in rows_summary:
        if row["pixels"]:
            logger.info("    %-18s %8d 像元  %8.2f km2  MF %.3f",
                        row["mineral_en"], row["pixels"], row["area_km2"], row["mean_mf"])
    return {"area": area, "status": "成功", "message": "分类 %d" % classified_total,
            "elapsed": elapsed}


def compare_sam(out_dir, area, logger):
    """读取步骤4 的 aloh_class.tif，与本步 mtmf_class.tif 逐像元对比。

    输出 compare_sam_vs_mtmf.csv（混淆矩阵、总体一致率），用于交叉验证两种方法。
    """
    sam_path = os.path.join(STEP4_ROOT, area, "aloh_class.tif")
    mtmf_path = os.path.join(out_dir, "mtmf_class.tif")
    if not (os.path.isfile(sam_path) and os.path.isfile(mtmf_path)):
        return
    ds_sam = gdal.Open(sam_path)
    ds_mtmf = gdal.Open(mtmf_path)
    sam = ds_sam.ReadAsArray()
    mtmf = ds_mtmf.ReadAsArray()
    ds_sam = None
    ds_mtmf = None
    valid = (sam != CLASS_NODATA) & (mtmf != CLASS_NODATA)
    sam_classified = valid & (sam > 0) & (sam != CLASS_NODATA)
    mtmf_classified = valid & (mtmf > 0) & (mtmf != CLASS_NODATA)
    both = sam_classified & mtmf_classified
    agreement = float(np.mean(sam[valid] == mtmf[valid])) if np.any(valid) else 0.0
    joint = float(np.mean(sam[both] == mtmf[both])) if np.any(both) else 0.0
    with open(os.path.join(out_dir, "compare_sam_vs_mtmf.csv"), "w", newline="",
              encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["metric", "value"])
        writer.writerow(["valid_pixels", int(valid.sum())])
        writer.writerow(["overall_agreement", "%.4f" % agreement])
        writer.writerow(["both_classified", int(both.sum())])
        writer.writerow(["agreement_when_both_classified", "%.4f" % joint])
        writer.writerow(["sam_classified", int(sam_classified.sum())])
        writer.writerow(["mtmf_classified", int(mtmf_classified.sum())])
    logger.info("[%s] SAM vs MTMF: 总体一致率 %.1f%%，双方都分类时一致率 %.1f%%",
                area, 100 * agreement, 100 * joint)


def make_qc(out_dir, area, samples, scores, counts, endmember_idx, wl_band, logger):
    """生成 preview_mtMF.png：MNF 分量、PPI 散点、端元光谱、MTMF 分类与 SAM 对比。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    def read_preview(name, band=1, dtype=np.float64):
        path = os.path.join(out_dir, name)
        ds = gdal.Open(path)
        width, height = ds.RasterXSize, ds.RasterYSize
        scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))
        bw, bh = max(1, width // scale), max(1, height // scale)
        if ds.RasterCount == 1:
            array = ds.ReadAsArray(buf_xsize=bw, buf_ysize=bh)
        else:
            band_handle = ds.GetRasterBand(band)
            array = band_handle.ReadAsArray(buf_xsize=bw, buf_ysize=bh)
            band_handle = None
        geo = ds.GetGeoTransform()
        extent = [geo[0] / 1000.0, (geo[0] + geo[1] * width) / 1000.0,
                  (geo[3] + geo[5] * height) / 1000.0, geo[3] / 1000.0]
        ds = None
        return array, extent

    mnf1, extent = read_preview("mnf_components.tif", 1)
    mnf2, _ = read_preview("mnf_components.tif", 2)
    class_map, _ = read_preview("mtmf_class.tif")
    colors = ["#dddddd"] + [m[2] for m in MINERALS]
    cmap = ListedColormap(colors)
    bounds = np.arange(-0.5, len(colors) + 0.5, 1.0)

    def stretch_display(array):
        data = np.where(array < -9000, np.nan, array.astype(np.float64))
        values = data[np.isfinite(data)]
        if values.size:
            low, high = np.percentile(values, [2, 98])
        else:
            low, high = 0.0, 1.0
        if high <= low:
            high = low + 1e-6
        return data, low, high

    figure = plt.figure(figsize=(18, 10))
    ax1 = figure.add_subplot(2, 3, 1)
    data1, low1, high1 = stretch_display(mnf1)
    ax1.imshow(data1, extent=extent, origin="upper", cmap="gray", vmin=low1, vmax=high1,
               interpolation="nearest")
    ax1.set_title("MNF component 1")
    ax2 = figure.add_subplot(2, 3, 2)
    data2, low2, high2 = stretch_display(mnf2)
    ax2.imshow(data2, extent=extent, origin="upper", cmap="gray", vmin=low2, vmax=high2,
               interpolation="nearest")
    ax2.set_title("MNF component 2")
    ax3 = figure.add_subplot(2, 3, 3)
    ax3.imshow(np.where(class_map == CLASS_NODATA, 0, class_map), extent=extent,
               origin="upper", cmap=cmap, norm=BoundaryNorm(bounds, cmap.N),
               interpolation="nearest")
    ax3.set_title("MTMF class")
    handles = [plt.Rectangle((0, 0), 1, 1, color=colors[i]) for i in range(1, len(colors))]
    ax3.legend(handles, MINERAL_EN, fontsize=6, ncol=2, loc="lower left", framealpha=0.85)

    sample_idx = np.random.default_rng(7).choice(scores.shape[0],
                                                 size=min(20000, scores.shape[0]),
                                                 replace=False)
    ax4 = figure.add_subplot(2, 3, 4, projection="3d")
    ax4.scatter(scores[sample_idx, 0], scores[sample_idx, 1], scores[sample_idx, 2],
                s=2, c="0.7", alpha=0.4, label="pixels")
    ax4.scatter(scores[endmember_idx, 0], scores[endmember_idx, 1],
                scores[endmember_idx, 2], s=60, c="red", marker="^", label="endmembers")
    ax4.set_xlabel("MNF1")
    ax4.set_ylabel("MNF2")
    ax4.set_zlabel("MNF3")
    ax4.set_title("n-D scatter (MNF1-3) + PPI endmembers")

    ax5 = figure.add_subplot(2, 3, 5)
    for ei, idx in enumerate(endmember_idx):
        spectrum = samples[idx]
        depth, pos = continuum_depth(wl_band, spectrum)
        ax5.plot(wl_band, spectrum, lw=0.8,
                 label="EM%d (PPI=%d, λc=%.2f)" % (ei + 1, counts[idx], pos))
    ax5.set_xlabel("Wavelength (um)")
    ax5.set_ylabel("Reflectance")
    ax5.set_title("Image endmember spectra (2.0-2.3 um)")
    ax5.legend(fontsize=6)
    ax5.grid(alpha=0.3)

    ax6 = figure.add_subplot(2, 3, 6)
    best_em = endmember_idx[0]
    continuum = convex_hull_continuum(wl_band, samples[best_em])
    ax6.plot(wl_band, samples[best_em], "k-", lw=1.0, label="EM1 spectrum")
    ax6.plot(wl_band, continuum, "r--", lw=1.0, label="convex-hull continuum")
    ax6.set_xlabel("Wavelength (um)")
    ax6.set_title("Envelope (continuum) removal")
    ax6.legend(fontsize=7)
    ax6.grid(alpha=0.3)

    for axis in (ax1, ax2, ax3):
        axis.set_xlabel("Easting (km)", fontsize=8)
        axis.set_ylabel("Northing (km)", fontsize=8)
        axis.tick_params(labelsize=7)

    figure.suptitle("%s | MNF + PPI + MTMF" % area, fontsize=13)
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    path = os.path.join(out_dir, "preview_mtMF.png")
    figure.savefig(path, dpi=QC_DPI, bbox_inches="tight")
    plt.close(figure)
    logger.info("[%s] 质检图: %s", area, path)


def main(argv=None):
    """主流程：解析参数 → 扫描区域 → 逐区处理 → 汇总。

    返回码：0 全部成功/跳过；1 有区域失败；2 输入目录或区域不满足条件。
    """
    args = parse_args(argv)
    out_root = os.path.abspath(args.out_root)
    os.makedirs(out_root, exist_ok=True)
    log_path = args.log or os.path.join(out_root, "_catalog", "step4_5_batch.log")
    logger = setup_logging(log_path, args.verbose)
    gdal.SetCacheMax(1024 * 1024 * 1024)

    logger.info("=" * 72)
    logger.info("ZY1E AHSI 步骤4.5：MNF + PPI/nD + 包络线 + MTMF")
    logger.info("步骤3 目录: %s", os.path.abspath(args.in_root))
    logger.info("输出目录  : %s", out_root)
    logger.info("=" * 72)

    try:
        areas = discover_areas(os.path.abspath(args.in_root), args.area)
    except Exception as exc:
        logger.error("扫描失败: %s", exc)
        return 2
    if not areas:
        logger.error("没有可处理区域")
        return 2

    if args.dry_run:
        for area in areas:
            try:
                inputs = resolve_area_inputs(os.path.abspath(args.in_root), area,
                                             out_root, False, logger, build=False)
            except Exception as exc:
                logger.error("  %-10s 输入不可用: %s", area, exc)
                continue
            if inputs["source"] == "mosaic":
                logger.info("  %-10s <- %s", area, inputs["refl_path"])
            else:
                logger.info("  %-10s <- %d 景裁剪成果（将合并为 %s）",
                            area, len(inputs["scenes"]), inputs["refl_path"])
        logger.info("dry-run 结束")
        return 0

    results = []
    batch_started = time.time()
    for index, area in enumerate(areas, start=1):
        logger.info("-" * 72)
        logger.info("(%d/%d) 处理: %s", index, len(areas), area)
        try:
            inputs = resolve_area_inputs(os.path.abspath(args.in_root), area,
                                         out_root, args.overwrite, logger)
            results.append(process_area(area, inputs, out_root, args, logger))
        except Exception as exc:
            logger.exception("[%s] 处理失败", area)
            results.append({"area": area, "status": "失败", "message": str(exc),
                            "elapsed": 0.0})

    logger.info("=" * 72)
    success = [r for r in results if r.get("status") == "成功"]
    skipped = [r for r in results if r.get("status") == "跳过"]
    failed = [r for r in results if r.get("status") == "失败"]
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %.1f s",
                len(success), len(skipped), len(failed), time.time() - batch_started)
    for result in results:
        logger.info("  [%s] %s - %s", result.get("status"), result["area"],
                    result.get("message", ""))
    if failed:
        logger.error("失败区域: %s", ", ".join(r["area"] for r in failed))
        logger.error("排查请查看日志: %s", log_path)
        return 1
    logger.info("完成，成果: %s", out_root)
    logger.info("日志文件: %s", log_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
