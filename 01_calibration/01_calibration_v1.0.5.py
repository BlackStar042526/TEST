# -*- coding: utf-8 -*-
"""
01_calibration_badbands.py
ZY1E AHSI 步骤1：辐射定标 + 坏波段剔除 + RPC 正射 + 质检（v1.0.5）

流程：读包内 VNIR/SWIR → 解析光谱/定标/足迹/RPC → 逐波段定标 L=Gain*DN+Bias
→ 按剔除表跳过 24 个坏波段 → RPC 正射到 EPSG:32650/30 m → 统计与质检。

本版仅补充逐行中文备注，算法与输出不变。
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

gdal.UseExceptions()                           # 出错即抛异常，不再静默返回
gdal.PushErrorHandler("CPLQuietErrorHandler")  # 屏蔽 GDAL 非致命错误输出

# =================================================================
# 配置区（按需修改）
# =================================================================
SRC_DIR = r"E:\ZhouYang\test\ZY1E_AHSI"           # 默认源目录：压缩包所在
OUT_ROOT = r"E:\ZhouYang\test\ZY1E_output_step1"  # 默认输出根目录
TAR_PATTERNS = ("*.tar", "*.tar.gz", "*.tgz")     # 扫描压缩包的通配模式

# 波段参数
N_VNIR = 76  # VNIR 波段数
N_SWIR = 90  # SWIR 波段数

# 坏波段剔除表（1-based）
REMOVE_BAND_GROUPS = [
    ("SWIR 首段低信噪", (77, 78, 79), 1006, 1039),                     # 组名、波段号、波长范围
    ("1.4 μm 水汽吸收", (98, 99, 100, 101), 1358, 1408),              # 水汽强吸收核心
    ("1.4 μm 邻近低信噪", (102, 103, 104, 105), 1425, 1476),           # 吸收带邻近低信噪
    ("1.812 μm 水汽吸收", (125,), 1812, 1812),                        # 单一吸收波段
    ("1.9 μm 水汽吸收", (126, 127, 128, 129, 130, 131), 1829, 1913),  # 1.9 μm 吸收带
    ("1.9 μm 邻近低信噪", (132, 133, 134), 1930, 1964),                # 吸收带邻近低信噪
    ("SWIR 末段低信噪", (164, 165, 166), 2468, 2501),                  # SWIR 末尾信噪比低
]
WATER_VAPOR_BANDS = [(1.35, 1.43), (1.80, 1.96)]  # 备选：按水汽区间剔除

CAL_HAS_BIAS = True  # 定标文件是否含 Bias 列
# 掩膜设置：如果确认影像中 DN=0 的像元为无效黑边，请将 MASK_NODATA 改为 True
MASK_NODATA = False   # 是否把 DN=0 当无效像元
NODATA_IN = 0         # 输入的无效像元值
NODATA_OUT = -9999.0  # 输出的无效值

# 正射参数
TARGET_EPSG = 32650           # 目标投影：UTM 50N
TARGET_RES = 30.0             # 目标分辨率（米）
RESAMPLING = "bilinear"       # 重采样算法
RPC_FIT_MAX_RESIDUAL = 2.0    # RPC 拟合残差阈值（像素）
REQUIRE_GEOREF = True         # 是否要求输入带地理参考
WRITE_GEOREF = True           # 是否输出 georef.json
KEEP_NATIVE_RADIANCE = False  # 是否保留未正射的辐亮度

# 批处理参数
OVERWRITE = False                     # 默认跳过已有成果
GDAL_CACHE_MB = 1024                  # GDAL 块缓存（MB）
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]  # GeoTIFF 基础创建选项
GEOTIFF_COMPRESS = "DEFLATE"          # 压缩算法
GEOTIFF_COMPRESS_LEVEL = 6            # 压缩级别 1~9

# 质检参数
MAKE_QC = True      # 是否生成质检图
QC_MAX_SIZE = 1400  # 质检图降采样最长边
QC_DPI = 150        # 质检图输出分辨率

# 输出文件名：以缩写标注定标、坏波段剔除和正射处理
OUT_ORTHO_SUFFIX = "_cal_badband_ortho.tif"     # 正射成果后缀
LEGACY_ORTHO_NAME = "radiance_ortho.tif"        # 旧版正射文件名，用于复用
OUT_NATIVE_NAME = "radiance_native.tif"         # 未正射辐亮度文件名
OUT_WAVELENGTHS_NAME = "wavelengths_clean.npy"  # 保留波段波长文件
OUT_GEOREF_NAME = "georef.json"                 # 几何信息文件名
OUT_STATS_NAME = "band_stats.csv"               # 逐波段统计文件名
OUT_QC_NAME = "preview_quality.png"             # 质检图文件名
LOG_NAME = "zy1e.step1"                         # 日志器名称

def _flatten_remove_bands(groups):
    bands = []                                         # 收集 1-based 波段号
    for _, indexes, _, _ in groups:                    # 只取每组的波段号元组
        bands.extend(int(index) for index in indexes)  # 逐个追加
    return bands                                       # 返回展平结果

REMOVE_BANDS_1BASED = _flatten_remove_bands(REMOVE_BAND_GROUPS)  # 展平全部剔除波段


# =================================================================
# 通用工具
# =================================================================
def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(description="ZY1E AHSI 步骤1 修复版")           # 创建参数解析器
    parser.add_argument("--src", default=SRC_DIR)                               # 源目录
    parser.add_argument("--out", default=OUT_ROOT)                              # 输出根目录
    parser.add_argument("--pattern", action="append", default=None)             # 自定义扫描通配符
    parser.add_argument("--scene", action="append", default=None)               # 只处理匹配景名
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE)  # 覆盖已有成果
    parser.add_argument("--limit", type=int, default=0)                         # 只处理前 N 景
    parser.add_argument("--keep-native", action="store_true", default=KEEP_NATIVE_RADIANCE)  # 保留未正射影像
    parser.add_argument("--no-qc", action="store_true")                         # 关闭质检图
    parser.add_argument("--dry-run", action="store_true")                       # 只打印计划不落盘
    parser.add_argument("--log", default=None)                                  # 自定义日志文件路径
    parser.add_argument("--verbose", action="store_true")                       # 控制台输出 DEBUG 日志
    return parser.parse_args(argv)                                              # 解析并返回参数命名空间

def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)                                        # 取具名日志器
    logger.setLevel(logging.DEBUG)                                              # 日志器总级别设为 DEBUG
    logger.handlers.clear()                                                     # 清空旧处理器，避免重复输出
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")  # 统一日志格式
    
    # 增加 flush=True 确保终端实时显示
    console = logging.StreamHandler(stream=sys.stdout)            # 控制台处理器
    console.setLevel(logging.DEBUG if verbose else logging.INFO)  # 级别随 --verbose
    console.setFormatter(fmt)                                     # 套用统一格式
    logger.addHandler(console)                                    # 挂载控制台处理器
    
    if log_path:  # 指定了日志文件才写文件
        try:
            os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)  # 确保日志目录存在
            file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")  # 追加写日志文件
            file_handler.setLevel(logging.DEBUG)                                # 文件记录最详细级别
            file_handler.setFormatter(fmt)                                      # 套用统一格式
            logger.addHandler(file_handler)                                     # 挂载文件处理器
        except OSError:                                                         # 日志目录不可写时忽略
            pass
    logger.propagate = False  # 不向根日志器传播
    return logger             # 返回配置好的日志器

def human_size(num_bytes):
    """把字节数格式化为 B/KB/MB/GB/TB 文本。"""
    size = float(num_bytes)                     # 转浮点便于逐级相除
    for unit in ("B", "KB", "MB", "GB", "TB"):  # 逐级换算单位
        if size < 1024.0 or unit == "TB":       # 不足一级或是最大单位
            return "%.1f %s" % (size, unit)     # 格式化后返回
        size /= 1024.0                          # 换算到下一级单位

def discover_scenes(src_dir, patterns):
    """扫描源目录下匹配后缀的压缩包，返回排序后的文件路径列表。"""
    if not os.path.isdir(src_dir):
        raise NotADirectoryError("源目录不存在: %s" % src_dir)                        # 源目录缺失直接报错
    found = set()                                                               # 存放去重后的路径
    for pattern in patterns:                                                    # 逐个通配符扫描
        found.update(glob.glob(os.path.join(src_dir, pattern), recursive=True))  # 递归匹配并合并
    return sorted(f for f in found if os.path.isfile(f))                        # 排序后只保留文件

def scene_name(tar_path):
    """由压缩包文件名推导景名（去掉 .tar/.tar.gz/.tgz 后缀）。"""
    name = os.path.basename(tar_path)                     # 只取文件名
    lower = name.lower()                                  # 小写化用于后缀判断
    for ext in (".tar.gz", ".tgz", ".tar"):               # 从长到短尝试后缀
        if lower.endswith(ext): return name[: -len(ext)]  # 命中则去掉后缀
    return os.path.splitext(name)[0]                      # 兜底：去掉最后扩展名

def find_tar_member(names, required, suffixes, label, logger):
    """在包内成员名中按“关键词全部包含 + 后缀命中”定位文件，返回首个匹配。"""
    hits = []                                                               # 命中的成员名
    for name in names:                                                      # 遍历包内成员
        base = os.path.basename(name).lower()                               # 只用文件名小写参与匹配
        if not base.endswith(suffixes): continue                            # 后缀不符则跳过
        if all(key.lower() in base for key in required): hits.append(name)  # 关键词全中含则记录
    if not hits:
        raise FileNotFoundError(f"压缩包中未找到 {label}（关键词 {required}，后缀 {suffixes}）")  # 无命中报错
    hits.sort()                                                                 # 排序保证选择稳定
    if len(hits) > 1:
        logger.warning("%s 匹配到 %d 个文件，使用第一个: %s", label, len(hits), os.path.basename(hits[0]))  # 多命中时告警
    return hits[0]                                                              # 返回首个匹配

def decode_text(raw):
    """按 UTF-8 宽松解码并统一换行、去掉 BOM。"""
    text = raw.decode("utf-8", errors="ignore")                             # 宽松解码，忽略非法字节
    return text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")  # 去 BOM 并统一换行

def read_member_text(tar, member, label):
    """读取 tar 内文本成员并解码；成员不可读时抛 FileNotFoundError。"""
    handle = tar.extractfile(member)                                       # 取得包内成员文件对象
    if handle is None: raise FileNotFoundError(f"{label} 无法读取: {member}")  # 取出失败即报错
    with handle: return decode_text(handle.read())                         # 读取并解码后返回


# =================================================================
# 产品元数据解析
# =================================================================
def parse_footprint(xml_text, logger):
    """从产品 XML 读取四个角点经纬度，返回 {"TL","TR","BR","BL": (lon, lat)}。"""
    """从产品 XML 解析四角经纬度，返回 {TL/TR/BR/BL: (lon, lat)}。"""
    footprint = {}                                                              # 角点缩写到经纬度的映射
    corner_map = {"TL": "TopLeft", "TR": "TopRight", "BR": "BottomRight", "BL": "BottomLeft"}  # 缩写对应 XML 标签
    for abbr, corner in corner_map.items():                                     # 逐个角点解析
        lat_match = re.search(r"<%sLatitude>\s*([-+0-9.eE]+)\s*</%sLatitude>" % (corner, corner), xml_text, re.IGNORECASE)  # 匹配纬度
        lon_match = re.search(r"<%sLongitude>\s*([-+0-9.eE]+)\s*</%sLongitude>" % (corner, corner), xml_text, re.IGNORECASE)  # 匹配经度
        if not lat_match or not lon_match: raise ValueError(f"XML 缺少 {corner} 经纬度")  # 缺字段即报错
        footprint[abbr] = (float(lon_match.group(1)), float(lat_match.group(1)))  # 按 (lon, lat) 存放
    return footprint                                                            # 返回四角坐标

def parse_rpb(rpb_text):
    """从 .rpb 文件中提取 RPC 参数，兼容 XML 和 key=value 格式。"""
    rpc_dict = {}                                                               # 归一化后的 RPC 参数字典
    for match in re.finditer(r"([a-zA-Z_]+)\s*=\s*(.*?);", rpb_text, re.DOTALL):  # 逐项匹配 key=value;
        key = match.group(1).lower()                                            # 键名统一小写
        val = match.group(2).strip()                                            # 值去掉首尾空白
        
        key_map = {  # rpb 键名到标准 RPC 键的映射
            "errbias": "ERR_BIAS", "errrand": "ERR_RAND",
            "lineoffset": "LINE_OFF", "sampleoffset": "SAMP_OFF", "sampoffset": "SAMP_OFF",
            "latoffset": "LAT_OFF", "longoffset": "LONG_OFF", "heightoffset": "HEIGHT_OFF",
            "linescale": "LINE_SCALE", "samplescale": "SAMP_SCALE", "sampscale": "SAMP_SCALE",
            "latscale": "LAT_SCALE", "longscale": "LONG_SCALE", "heightscale": "HEIGHT_SCALE",
            "linenumcoef": "LINE_NUM_COEFF", "linedencoef": "LINE_DEN_COEFF",
            "sampnumcoef": "SAMP_NUM_COEFF", "sampdencoef": "SAMP_DEN_COEFF"
        }
        if key in key_map:                                                      # 只保留已知参数
            if "coef" in key:                                                   # 系数项需清理分隔符
                val = val.replace("(", "").replace(")", "").replace("\n", " ").replace("\r", " ").replace(",", " ")  # 去括号与换行
                val = re.sub(r"\s+", " ", val).strip()                          # 压缩多余空白
            rpc_dict[key_map[key]] = val                                        # 写入标准键名
            
    if not rpc_dict:                                                            # 未解析到则尝试 XML 形式
        for key in ["ERR_BIAS", "ERR_RAND", "LINE_OFF", "SAMP_OFF", "LAT_OFF", "LONG_OFF",  # 逐个标准键查找
                    "HEIGHT_OFF", "LINE_SCALE", "SAMP_SCALE", "LAT_SCALE", "LONG_SCALE", "HEIGHT_SCALE",
                    "LINE_NUM_COEFF", "LINE_DEN_COEFF", "SAMP_NUM_COEFF", "SAMP_DEN_COEFF"]:
            match = re.search(r"<%s>\s*([^<]+)\s*</%s>" % (key, key), rpb_text, re.IGNORECASE)  # 匹配 <KEY>值</KEY>
            if match: rpc_dict[key] = match.group(1).strip()                    # 命中则记录
    return rpc_dict                                                             # 返回 RPC 字典

def parse_wavelengths(text, n_bands, label, logger):
    """解析光谱响应文件的中心波长（nm→μm），返回前 n_bands 个值。"""
    """解析光谱响应文件，取每行第一列（或第二列）中心波长并换算为 μm。"""
    values = []                                        # 收集中心波长（nm）
    for line in text.splitlines():                     # 逐行解析
        line = line.strip()                            # 去掉首尾空白
        if not line or line.startswith("#"): continue  # 跳过空行与注释行
        parts = re.split(r"[,\s]+", line)              # 按逗号或空白切分
        if not parts: continue                         # 无字段则跳过
        try:
            if len(parts) >= 2 and float(parts[0]).is_integer() and float(parts[0]) < 200:  # 首列疑似波段号
                values.append(float(parts[1]))                                  # 取第二列作为波长
            else:
                values.append(float(parts[0]))                               # 否则首列即波长
        except ValueError: continue                                          # 非数字行直接跳过
    if len(values) < n_bands: raise ValueError(f"{label} 数据不足 {n_bands} 行")  # 行数不足报错
    wavelengths = np.asarray(values[:n_bands], dtype=np.float64) / 1000.0    # 截取并 nm 转 μm
    if not np.all(np.isfinite(wavelengths)) or np.any(wavelengths <= 0):
        raise ValueError(f"{label} 波长异常（存在非正数或 NaN）")  # 非法波长报错
    return wavelengths                                 # 返回波长数组

def parse_calibration(text, n_bands, has_bias, label, logger):
    """解析辐射定标文件，返回 (gain, bias) 两个数组，长度均为 n_bands。"""
    """解析辐射定标系数，返回 (gain 数组, bias 数组)，各取前 n_bands 行。

    has_bias=True 时每行需为 [波段号?, Gain, Bias]；否则 Bias 置 0。
    """
    rows = []                                                                   # 每行解析出的数字
    for line in text.splitlines():                                              # 逐行解析定标文件
        line = line.strip()                                                     # 去掉首尾空白
        if not line or line.startswith("#"): continue                           # 跳过空行与注释
        tokens = re.split(r"[,\s]+", line)                                      # 按逗号或空白切分
        numbers = []                                                            # 本行数字
        for token in tokens:                                                    # 逐个字段转数字
            try: numbers.append(float(token))                                   # 尝试转浮点
            except ValueError: break                                            # 遇非数字即停止
        if numbers:                                                             # 有数字才处理
            if len(numbers) >= 3 and numbers[0].is_integer() and numbers[0] < 200:  # 首列是波段号
                rows.append(numbers[1:3])                                       # 取 Gain 与 Bias
            elif len(numbers) == 2 and numbers[0].is_integer() and numbers[0] < 200:  # 波段号 + Gain
                rows.append([numbers[1]])                                       # 先只记 Gain
            else:
                rows.append(numbers[:2])  # 无波段号，取前两列
                
    if not rows: raise ValueError(f"{label} 定标系数为空")                            # 空文件报错
    if has_bias is None:                                                        # 未显式指定则按行长度推断
        has_bias = all(len(row) >= 2 for row in rows)                           # 全都有两列才含 Bias
    if has_bias:                                                                # 需要 Bias 列
        bad_rows = [i for i, row in enumerate(rows) if len(row) < 2]            # 找出缺 Bias 的行
        if bad_rows: raise ValueError(f"{label} 定标文件第 {bad_rows[0]+1} 行缺 Bias")  # 缺列报错
        data = np.asarray([row[:2] for row in rows], dtype=np.float64)          # 组成 (N,2) 数组
    else:
        data = np.asarray([[row[0], 0.0] for row in rows], dtype=np.float64)  # Bias 统一补 0
        
    if len(data) < n_bands: raise ValueError(f"{label} 定标系数不足 {n_bands} 行")  # 行数不足报错
    return data[:n_bands, 0], data[:n_bands, 1]                              # 返回前 n 波段的 gain/bias


# =================================================================
# 坏波段剔除与几何处理
# =================================================================
def resolve_band_selection(wavelengths, logger):
    """按 REMOVE_BAND_GROUPS（1-based）生成保留波段索引、剔除列表与分组信息。"""
    """按剔除表（0-based 转换）返回 (保留波段索引, 被剔除索引, 分组说明)。"""
    n_total = int(wavelengths.size)                                            # 总波段数
    if REMOVE_BANDS_1BASED:                                                    # 优先使用 1-based 剔除表
        invalid = [b for b in REMOVE_BANDS_1BASED if not (1 <= b <= n_total)]  # 检查越界波段号
        if invalid: raise ValueError(f"波段号越界: {invalid[:10]}")                 # 越界即报错
        removed = sorted({int(b) - 1 for b in REMOVE_BANDS_1BASED})            # 转 0-based 并去重排序
        groups = []                                                            # 分组说明列表
        for label, indexes, wl_lo, wl_hi in REMOVE_BAND_GROUPS:                # 逐组转换
            idx0 = sorted({int(b) - 1 for b in indexes})                       # 组内波段号转 0-based
            if not idx0: continue                                              # 空组跳过
            groups.append({"label": label, "bands_1based": list(indexes),
                           "wavelength_nm": [float(wavelengths[idx0[0]]) * 1000.0, float(wavelengths[idx0[-1]]) * 1000.0]})  # 记录组首末波长
    else:
        removed = [idx for idx, wl in enumerate(wavelengths) if any(lo <= wl <= hi for lo, hi in WATER_VAPOR_BANDS)]  # 按水汽区间剔除
        groups = []                                                             # 无分组信息
    keep = np.asarray([i for i in range(n_total) if i not in set(removed)], dtype=np.int64)  # 计算保留索引
    if keep.size <= 0: raise ValueError("剔除后没有剩余波段")                            # 全部剔除则报错
    return keep, removed, groups                                                # 返回保留/剔除/分组

def open_vsitar_dataset(tar_path, internal_path, label):
    """通过 /vsitar/ 直接打开压缩包内影像，避免解压到磁盘。"""
    """用 /vsitar/ 直接打开压缩包内影像，避免整包解压。"""
    vsi_path = "/vsitar/" + os.path.abspath(tar_path).replace("\\", "/") + "/" + internal_path  # 拼 vsitar 虚拟路径
    dataset = gdal.Open(vsi_path)                                               # 打开包内影像
    if dataset is None: raise RuntimeError(f"无法通过 /vsitar/ 打开 {label}")         # 打开失败报错
    return dataset                                                              # 返回数据集对象

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
    rpc = ds_vn.GetMetadata("RPC")                # 读取数据集 RPC 元数据
    if not rpc: raise ValueError("影像缺少 RPC 元数据")  # 无 RPC 无法校正
    
    gt = ds_vn.GetGeoTransform()                                         # 读取仿射变换
    if not (abs(abs(gt[1]) - 1) < 1e-6 and abs(abs(gt[5]) - 1) < 1e-6):  # 检查是否单位阵
        raise ValueError(f"交付影像 GeoTransform 非单位阵，RPC 校正失效: {gt}")       # 非单位阵时校正无意义

    transformer = gdal.Transformer(ds_vn, None, ["METHOD=RPC"])                 # 建立 RPC 正反变换器
    width, height = ds_vn.RasterXSize, ds_vn.RasterYSize                        # 影像宽高
    pixel_corners = np.array([[0.0, 0.0], [width, 0.0], [width, height], [0.0, height]])  # 四角像素坐标
    raw = []                                                                    # 反变换得到的真实行列
    for corner in ("TL", "TR", "BR", "BL"):                                     # 按四角顺序处理
        lon, lat = footprint[corner]                                            # 取角点经纬度
        ok, point = transformer.TransformPoint(1, lon, lat)                     # 经纬度转像素行列
        if not ok: raise ValueError(f"RPC 反变换失败: {corner}")                     # 反变换失败报错
        raw.append([point[0], point[1]])                                        # 记录行列
    raw = np.asarray(raw, dtype=np.float64)                                     # 转成数组

    (sx, tx) = np.linalg.lstsq(np.column_stack([pixel_corners[:, 0], np.ones(4)]), raw[:, 0], rcond=None)[0]  # 列向最小二乘拟合
    (sy, ty) = np.linalg.lstsq(np.column_stack([pixel_corners[:, 1], np.ones(4)]), raw[:, 1], rcond=None)[0]  # 行向最小二乘拟合
    predicted = np.column_stack([pixel_corners[:, 0] * sx + tx, pixel_corners[:, 1] * sy + ty])  # 拟合预测行列
    residual = float(np.max(np.hypot(predicted[:, 0] - raw[:, 0], predicted[:, 1] - raw[:, 1])))  # 取最大残差
    
    if residual > RPC_FIT_MAX_RESIDUAL: raise ValueError(f"RPC 拟合残差过大: {residual}")  # 残差超限报错
    corrected = dict(rpc)                                                       # 复制一份 RPC 参数
    corrected["SAMP_OFF"] = str((float(rpc["SAMP_OFF"]) - tx) / sx)             # 修正列偏移
    corrected["SAMP_SCALE"] = str(float(rpc["SAMP_SCALE"]) / sx)                # 修正列比例
    corrected["LINE_OFF"] = str((float(rpc["LINE_OFF"]) - ty) / sy)             # 修正行偏移
    corrected["LINE_SCALE"] = str(float(rpc["LINE_SCALE"]) / sy)                # 修正行比例
    return corrected, {"sx": sx, "tx": tx, "sy": sy, "ty": ty, "residual": residual}  # 返回修正 RPC 与拟合参数

def write_cutline_vsimem(footprint, logger):
    """把四角点足迹写成 /vsimem/ 下的 GeoJSON，作为 gdal.Warp 的裁剪线。"""
    """把四角足迹写成 /vsimem/ 下的 GeoJSON 裁切线，供 gdal.Warp 使用。"""
    ring = [footprint[key] for key in ("TL", "TR", "BR", "BL")] + [footprint["TL"]]  # 闭合外环
    geojson = {"type": "FeatureCollection", "features": [{"type": "Feature", "properties": {},  # 组装单要素面
               "geometry": {"type": "Polygon", "coordinates": [[[float(lon), float(lat)] for lon, lat in ring]]}}]}  # 顶点转浮点
    path = "/vsimem/zy1e_step1_footprint.geojson"                               # 内存文件路径
    gdal.FileFromMemBuffer(path, json.dumps(geojson).encode("utf-8"))           # 写入内存文件
    return path                                                                 # 返回虚拟路径

def build_gtiff_creation_options(band_interleave=True):
    """组装 GeoTIFF 创建选项（压缩、预测器、分块等）。"""
    """组装 GeoTIFF 创建选项：BIGTIFF、按波段存储、DEFLATE+PREDICTOR=3。"""
    creation = list(GTIFF_OPTIONS)                                              # 复制基础选项
    if band_interleave: creation.append("INTERLEAVE=BAND")                      # 按波段独立存储
    if GEOTIFF_COMPRESS: creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS, "PREDICTOR=3", "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL])  # 压缩及预测器
    return creation                                                             # 返回创建选项

def orthorectify_radiance(native_path, footprint, ortho_tmp, logger):
    """用 RPC + 足迹裁剪线把辐亮度正射到 TARGET_EPSG/TARGET_RES，返回网格信息。"""
    cutline = write_cutline_vsimem(footprint, logger)  # 生成足迹裁剪线
    try:
        options = gdal.WarpOptions(                                             # 组装 Warp 参数
            format="GTiff", dstSRS="EPSG:%d" % TARGET_EPSG, xRes=TARGET_RES, yRes=TARGET_RES,  # 输出格式与网格
            resampleAlg=RESAMPLING, dstNodata=NODATA_OUT, targetAlignedPixels=True,  # 重采样与无效值
            multithread=True, creationOptions=build_gtiff_creation_options(),   # 多线程与创建选项
            rpc=True, cutlineDSName=cutline, cropToCutline=True, outputType=gdal.GDT_Float32)  # 启用 RPC 并裁边
        result = gdal.Warp(ortho_tmp, native_path, options=options)             # 执行正射写临时文件
        if result is None: raise RuntimeError("RPC 正射失败")                       # Warp 失败报错
        width, height = result.RasterXSize, result.RasterYSize                  # 成果宽高
        geo_transform = result.GetGeoTransform()                                # 成果仿射变换
        band_count = result.RasterCount                                         # 成果波段数
        result = None                                                           # 关闭数据集释放句柄
    finally:
        try: gdal.Unlink(cutline)  # 删除内存裁剪线
        except: pass               # 删除失败也不影响主流程

    bounds = [geo_transform[0], geo_transform[3] + geo_transform[5] * height, geo_transform[0] + geo_transform[1] * width, geo_transform[3]]  # 计算外接范围
    return {"width": int(width), "height": int(height), "band_count": int(band_count),
            "geotransform": [float(v) for v in geo_transform], "bounds": [float(v) for v in bounds]}  # 返回网格信息


# =================================================================
# 统计输出与质检成图
# =================================================================
def write_georef_json(output_path, scene, tar_path, vn_tif, sw_tif, native_size, ortho_info, band_count, bands_source, footprint, rpc_correction, wavelengths, keep, removed, remove_groups, logger):
    """写 georef.json（网格、足迹、RPC 修正、剔除波段），供步骤2/3 追溯。"""
    payload = {                                                                 # 汇总本景几何与剔除信息
        "scene": scene, "source_tar": os.path.basename(tar_path), "source_tiffs": [os.path.basename(vn_tif), os.path.basename(sw_tif)],  # 场景与来源
        "geometry": "RPC orthorectified", "target_crs_epsg": TARGET_EPSG, "target_res": TARGET_RES, "resampling": RESAMPLING,  # 目标投影与重采样
        "native_raster_size": [int(native_size[0]), int(native_size[1])], "raster_size": [int(ortho_info["width"]), int(ortho_info["height"])],  # 原图与成果尺寸
        "band_count": int(band_count), "bands_source": int(bands_source), "geotransform": ortho_info["geotransform"], "bounds": ortho_info["bounds"],  # 波段数与网格
        "footprint_wgs84": footprint, "rpc_grid_correction": rpc_correction, "remove_band_groups": remove_groups,  # 足迹、RPC 修正与剔除分组
        "removed_bands_1based": [int(i) + 1 for i in removed], "removed_wavelengths_um": [float(wavelengths[i]) for i in removed],  # 被剔除波段与波长
        "radiance_formula": "L = Gain * DN + Bias", "nodata": NODATA_OUT, "mask_nodata": bool(MASK_NODATA)  # 定标公式与无效值
    }
    with open(output_path + ".tmp", "w", encoding="utf-8") as f: json.dump(payload, f, ensure_ascii=False, indent=2)  # 先写临时文件
    os.replace(output_path + ".tmp", output_path)                               # 原子替换为正式文件

def write_band_stats_csv(stats, output_path):
    """把逐波段统计写成 CSV（band_stats.csv）。"""
    tmp_path = output_path + ".tmp"                                             # 临时文件路径
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:           # 以 UTF-8 新建
        writer = csv.writer(handle)                                             # 创建 CSV 写入器
        writer.writerow(["band", "wavelength_um", "min", "max", "mean", "std", "p1", "p99", "valid_pixels"])  # 写表头
        for item in stats:                                                      # 逐波段写一行
            writer.writerow([item["band"], "%.6f" % item["wavelength"], "%.6f" % item["min"], "%.6f" % item["max"],
                             "%.6f" % item["mean"], "%.6f" % item["std"], "%.6f" % item["p1"], "%.6f" % item["p99"], item["valid"]])  # 指标格式化输出
    os.replace(tmp_path, output_path)                                           # 原子替换为正式文件

def _nearest_band(wavelengths, target):
    """返回最接近 target（μm）的波段下标。"""
    return int(np.argmin(np.abs(wavelengths - target)))  # 取距离最近的波段
def select_qc_bands(wavelengths):
    """为质检图挑选真彩色、假彩色与 SWIR 所用波段下标。"""
    return {"R": _nearest_band(wavelengths, 0.66), "G": _nearest_band(wavelengths, 0.55), "B": _nearest_band(wavelengths, 0.47),
            "NIR": _nearest_band(wavelengths, 0.86), "SWIR": _nearest_band(wavelengths, 2.20)}, None  # 各通道对应波段

def _stretch(array, valid):
    """2%–98% 线性拉伸到 0~1，仅用 valid 像元估计分位数。"""
    values = array[valid]                                   # 只取有效像元
    if values.size == 0: return np.zeros_like(array)        # 无有效值返回全 0
    low, high = np.percentile(values, [2, 98])              # 估 2%/98% 分位数
    if high <= low: high = low + 1e-6                       # 防止除零
    return np.clip((array - low) / (high - low), 0.0, 1.0)  # 线性拉伸并裁剪到 0~1

def load_qc_arrays(tif_path, band_map, unique):
    """降采样读取质检所需波段，返回 (数组字典, 绘图范围信息)。"""
    dataset = gdal.Open(tif_path)                                          # 打开正射成果
    if dataset is None: return None                                        # 打不开直接返回
    width, height = dataset.RasterXSize, dataset.RasterYSize               # 原图宽高
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))  # 计算降采样倍数
    buf_w, buf_h = max(1, width // scale), max(1, height // scale)         # 缩略图尺寸
    
    arrays = {}                                                                 # 通道名到数组的映射
    for name, band_index in band_map.items():                                   # 逐个所需波段读取
        position = int(band_index) + 1                                          # GDAL 波段号从 1 开始
        arrays[name] = dataset.GetRasterBand(position).ReadAsArray(buf_xsize=buf_w, buf_ysize=buf_h).astype(np.float64)  # 降采样读入
        
    geo_transform = dataset.GetGeoTransform()  # 读取仿射变换
    extent = [geo_transform[0] / 1000.0, (geo_transform[0] + geo_transform[1] * width) / 1000.0,
              (geo_transform[3] + geo_transform[5] * height) / 1000.0, geo_transform[3] / 1000.0]  # 范围换算为 km
    info = {"extent": extent, "xlabel": "Easting (km)", "ylabel": "Northing (km)", "size_text": f"{width}x{height}"}  # 绘图信息
    return arrays, info                                                         # 返回数组与范围

def write_quality_figure(tif_path, wavelengths, removed_ranges, output_png, scene, logger, band_stats=None):
    """生成 preview_quality.png：真彩色 / 假彩色 / SWIR / 均值光谱（灰带为剔除波段）。"""
    try:
        import matplotlib
        matplotlib.use("Agg")  # 无界面后端，便于批量出图
        import matplotlib.pyplot as plt
    except ImportError: 
        logger.error("Matplotlib 未安装或导入失败，无法生成质检图。")  # 缺少依赖时告警
        return False
    
    band_map, _ = select_qc_bands(wavelengths)           # 选定质检波段
    prepared = load_qc_arrays(tif_path, band_map, None)  # 降采样读取波段
    if prepared is None: 
        logger.error("无法读取正射成果用于生成质检图。")  # 读不到成果
        return False
    arrays, info = prepared  # 拆包数组与范围
    
    def combine(names):
        channels = [arrays[name] for name in names]                          # 取各通道数组
        valid = np.ones_like(channels[0], dtype=bool)                        # 有效像元掩膜
        for ch in channels: valid &= np.isfinite(ch) & (ch != NODATA_OUT)    # 逐通道累积有效
        rgb = np.ones(channels[0].shape + (3,), dtype=np.float64)            # 白色底图
        for i, ch in enumerate(channels): rgb[..., i] = _stretch(ch, valid)  # 逐通道拉伸
        rgb[~valid] = 1.0                                                    # 无效像元填白
        return rgb, valid                                                    # 返回图像与掩膜

    true_rgb, valid = combine(("R", "G", "B"))             # 合成真彩色
    false_rgb, _ = combine(("NIR", "R", "G"))              # 合成假彩色
    swir = arrays["SWIR"].copy()                           # 取 SWIR 通道
    swir_valid = np.isfinite(swir) & (swir != NODATA_OUT)  # SWIR 有效掩膜
    swir_stretched = _stretch(swir, swir_valid)            # SWIR 拉伸
    swir_stretched[~swir_valid] = np.nan                   # 无效像元留空

    figure, axes = plt.subplots(2, 2, figsize=(13, 10)) if band_stats else plt.subplots(1, 3, figsize=(15, 5.5))  # 有统计时 2x2
    axes = axes.ravel()                                                         # 拉平轴数组
    axes[0].imshow(true_rgb, extent=info["extent"], origin="upper"); axes[0].set_title("True color")  # 真彩色子图
    axes[1].imshow(false_rgb, extent=info["extent"], origin="upper"); axes[1].set_title("False color")  # 假彩色子图
    axes[2].imshow(swir_stretched, extent=info["extent"], origin="upper", cmap="gray"); axes[2].set_title("SWIR")  # SWIR 子图
    
    if band_stats:                                                              # 有统计时追加光谱图
        wl_arr = np.array([i["wavelength"] for i in band_stats])                # 波长数组
        mean_arr = np.array([i["mean"] for i in band_stats])                    # 均值数组
        axes[3].plot(wl_arr, mean_arr, color="black", linewidth=1.0)            # 绘制均值光谱
        for low, high in removed_ranges: axes[3].axvspan(low, high, color="0.88", zorder=0)  # 灰带标注剔除区间
        axes[3].set_title("Mean spectrum (shaded = removed)")                   # 子图标题
        
    figure.savefig(output_png, dpi=QC_DPI, bbox_inches="tight")  # 保存 PNG
    plt.close(figure)                                            # 释放画布
    return True                                                  # 生成成功


# =================================================================
# 单景核心处理逻辑
# =================================================================
def resolve_scene_paths(scene_dir):
    """按场景目录拼出本步各成果文件路径。"""
    scene = os.path.basename(os.path.normpath(scene_dir))  # 取场景名
    return {k: os.path.join(scene_dir, v) for k, v in {    # 拼接各成果路径
        "ortho": scene + OUT_ORTHO_SUFFIX, "native": OUT_NATIVE_NAME, "wavelengths": OUT_WAVELENGTHS_NAME,
        "georef": OUT_GEOREF_NAME, "stats": OUT_STATS_NAME, "qc": OUT_QC_NAME}.items()}  # 键值对展开

def process_scene(tar_path, out_root, args, logger):
    """处理单景：①元数据 ②打开影像/RPC ③波段规划 ④定标写盘

    ⑤RPC 正射 ⑥原子替换成果 ⑦元数据与统计 ⑧质检图。
    单景异常向上抛出，由 main 记录堆栈并计入失败，不中断批处理。
    """
    scene = scene_name(tar_path)               # 由包名推导景名
    scene_dir = os.path.join(out_root, scene)  # 本景输出目录
    paths = resolve_scene_paths(scene_dir)     # 本景各成果路径
    make_qc = MAKE_QC and not args.no_qc       # 汇总质检开关

    required_files = ["ortho", "wavelengths", "georef", "stats"]  # 判定完成所需成果
    if make_qc:
        required_files.append("qc")                            # 需质检图时一并检查
    legacy_ortho = os.path.join(scene_dir, LEGACY_ORTHO_NAME)  # 旧版正射文件名
    if (not args.overwrite and not os.path.isfile(paths["ortho"])
            and os.path.isfile(legacy_ortho)
            and all(os.path.isfile(paths[key])
                    for key in required_files if key != "ortho")):     # 旧成果齐全则复用
        os.replace(legacy_ortho, paths["ortho"])                       # 旧名改为新名
        logger.info("复用已有完整成果并更新正射文件名：%s", paths["ortho"])             # 记录复用
    have_main = all(os.path.isfile(paths[p]) for p in required_files)  # 主成果是否齐全

    if not args.overwrite and have_main:
        logger.info("跳过（成果已存在）: %s", scene_dir)                                 # 成果齐全则跳过
        return {"scene": scene, "status": "跳过", "bands": None, "elapsed": 0.0, "message": "已存在"}  # 返回跳过结果

    logger.info("开始处理场景: %s", scene)                                           # 记录开始
    started = time.time()                                                      # 计时起点
    created_dir = not os.path.isdir(scene_dir)                                 # 记录目录是否本次新建
    os.makedirs(scene_dir, exist_ok=True)                                      # 创建输出目录
    tmp_native = os.path.join(scene_dir, "radiance_native.tmp.tif")            # 定标临时文件
    tmp_ortho = os.path.join(scene_dir, scene + "_cal_badband_ortho.tmp.tif")  # 正射临时文件
    tmp_wl = os.path.join(scene_dir, "wavelengths_clean.tmp.npy")              # 波长临时文件
    ds_vn = ds_swir = None                                                     # 占位，便于 finally 释放
    
    try:
        # ① 读取并解析元数据：光谱响应、定标系数、足迹与坏波段表
        logger.info("  - 正在读取并解析元数据...")                                        # 进度提示
        with tarfile.open(tar_path, "r:*") as tar:                              # 自动识别压缩格式
            names = tar.getnames()                                              # 包内成员名列表
            vn_tif = find_tar_member(names, ("_vn",), (".tif", ".tiff"), "VNIR 影像", logger)  # 定位 VNIR 影像
            sw_tif = find_tar_member(names, ("_sw",), (".tif", ".tiff"), "SWIR 影像", logger)  # 定位 SWIR 影像
            vn_rpb = find_tar_member(names, ("_vn",), (".rpb",), "VNIR RPC 文件", logger)  # 定位 VNIR RPC
            vn_cal = find_tar_member(names, ("vnir", "radcal"), (".raw",), "VNIR 定标", logger)  # 定位 VNIR 定标
            sw_cal = find_tar_member(names, ("swir", "radcal"), (".raw",), "SWIR 定标", logger)  # 定位 SWIR 定标
            vn_sr = find_tar_member(names, ("vnir", "spectralresponse"), (".raw",), "VNIR 光谱", logger)  # 定位 VNIR 光谱
            sw_sr = find_tar_member(names, ("swir", "spectralresponse"), (".raw",), "SWIR 光谱", logger)  # 定位 SWIR 光谱
            xml_member = find_tar_member(names, (), (".xml",), "元数据 XML", logger)  # 定位元数据 XML

            wavelengths = np.concatenate([                                      # 合并两段波长
                parse_wavelengths(read_member_text(tar, vn_sr, "VNIR 光谱"), N_VNIR, "VNIR", logger),  # VNIR 波长
                parse_wavelengths(read_member_text(tar, sw_sr, "SWIR 光谱"), N_SWIR, "SWIR", logger)])  # SWIR 波长
            vn_gain, vn_bias = parse_calibration(read_member_text(tar, vn_cal, "VNIR 定标"), N_VNIR, CAL_HAS_BIAS, "VNIR", logger)  # 解析 VNIR 系数
            sw_gain, sw_bias = parse_calibration(read_member_text(tar, sw_cal, "SWIR 定标"), N_SWIR, CAL_HAS_BIAS, "SWIR", logger)  # 解析 SWIR 系数
            keep, removed, remove_groups = resolve_band_selection(wavelengths, logger)  # 计算保留波段
            footprint = parse_footprint(read_member_text(tar, xml_member, "XML"), logger)  # 解析四角足迹

            vn_rpc_dict = {}  # 外部 .rpb 解析结果
            try:
                vn_rpb_text = read_member_text(tar, vn_rpb, "VNIR RPC")  # 读取 .rpb 文本
                vn_rpc_dict = parse_rpb(vn_rpb_text)                     # 解析 RPC 参数
            except Exception as e:
                logger.warning("手动解析外部 .rpb 失败，将依赖 TIFF 内自带的 RPC: %s", e)  # 解析失败仅告警

        # ② 打开 /vsitar/ 影像；外部 .rpb 解析失败时依赖 TIFF 自带 RPC
        logger.info("  - 正在打开影像并计算 RPC 校正参数...")                 # 进度提示
        ds_vn = open_vsitar_dataset(tar_path, vn_tif, "VNIR")    # 打开 VNIR 数据集
        ds_swir = open_vsitar_dataset(tar_path, sw_tif, "SWIR")  # 打开 SWIR 数据集
        
        if not ds_vn.GetMetadata("RPC") and vn_rpc_dict:  # 数据内无 RPC 时注入
            ds_vn.SetMetadata(vn_rpc_dict, "RPC")         # 写入 RPC 元数据
            logger.info("    已注入外部解析的 RPC 参数至数据集")        # 记录注入结果

        rpc_corrected, rpc_correction = compute_rpc_grid_correction(ds_vn, footprint, logger)  # 计算 RPC 网格修正
        
        n_rows, n_cols = ds_vn.RasterYSize, ds_vn.RasterXSize   # 原始行列数
        n_original, n_output = N_VNIR + N_SWIR, int(keep.size)  # 原波段数与输出波段数

        # ③ 规划输出波段：剔除波段直接跳过，按全局波段号取对应 VNIR/SWIR 定标系数
        band_plan = []                                                          # 输出波段处理计划
        for global_idx in keep:                                                 # 遍历保留的全局波段号
            if global_idx < N_VNIR:                                             # 属于 VNIR
                band_plan.append((ds_vn, int(global_idx) + 1, vn_gain[global_idx], vn_bias[global_idx], wavelengths[global_idx]))  # 记录来源与系数
            else:
                local = int(global_idx) - N_VNIR                                # 换算为 SWIR 内部序号
                band_plan.append((ds_swir, local + 1, sw_gain[local], sw_bias[local], wavelengths[global_idx]))  # 记录 SWIR 计划

        driver = gdal.GetDriverByName("GTiff")                                  # 获取 GeoTIFF 驱动
        if os.path.exists(tmp_native): os.remove(tmp_native)                    # 清理旧临时文件
        out_ds = driver.Create(tmp_native, n_cols, n_rows, n_output, gdal.GDT_Float32, options=build_gtiff_creation_options())  # 创建输出数据集
        out_ds.SetGeoTransform(ds_vn.GetGeoTransform())                         # 继承仿射变换
        out_ds.SetProjection(ds_vn.GetProjection())                             # 继承投影信息
        out_ds.SetMetadata(rpc_corrected, "RPC")                                # 写入修正后的 RPC

        # ④ 逐波段辐射定标：L = Gain * DN + Bias（可选 DN=0 掩膜），并统计
        logger.info("  - 开始逐波段执行辐射定标与坏波段剔除 (共 %d 个波段)...", n_output)            # 进度提示
        calibration_started = time.time()                                       # 定标计时起点
        stats = []                                                              # 逐波段统计结果
        for out_idx, (src_ds, src_band, gain, bias, wavelength) in enumerate(band_plan, start=1):  # 逐输出波段处理
            dn = src_ds.GetRasterBand(src_band).ReadAsArray()                   # 读取原始 DN
            radiance = dn.astype(np.float64) * float(gain) + float(bias)        # 辐射定标 L=G*DN+B
            if MASK_NODATA:                                                     # 需要掩膜时
                radiance[dn == NODATA_IN] = NODATA_OUT                          # DN=0 置为无效
            
            out_band = out_ds.GetRasterBand(out_idx)          # 取输出波段
            out_band.WriteArray(radiance.astype(np.float32))  # 按 Float32 写盘
            out_band.SetNoDataValue(NODATA_OUT)               # 设置无效值
            
            valid_values = radiance[radiance != NODATA_OUT] if MASK_NODATA else radiance[np.isfinite(radiance)]  # 提取有效像元
            if valid_values.size:                                               # 存在有效像元
                p1, p99 = np.percentile(valid_values[::4], [1, 99])             # 抽样估 1%/99% 分位
                stats.append({"band": out_idx, "wavelength": float(wavelength), "min": float(valid_values.min()),
                              "max": float(valid_values.max()), "mean": float(valid_values.mean()),
                              "std": float(valid_values.std()), "p1": float(p1), "p99": float(p99), "valid": int(valid_values.size)})  # 记录统计
            else:
                stats.append({"band": out_idx, "wavelength": float(wavelength), "min": 0.0, "max": 0.0, "mean": 0.0, "std": 0.0, "p1": 0.0, "p99": 0.0, "valid": 0})  # 全无效填 0
            
            # 每 10 个波段打印一次进度
            if out_idx == 1 or out_idx % 10 == 0 or out_idx == n_output:        # 首尾及每 10 个波段
                logger.info("    定标与写盘进度: %d/%d (%.1f%%)", out_idx, n_output, 100.0 * out_idx / n_output)  # 打印进度
                
        out_ds.FlushCache()  # 刷新缓存到磁盘
        out_ds = None        # 关闭数据集
        
        # ⑤ RPC 正射：用足迹裁剪线裁掉黑边，输出 EPSG:32650 / 30 m
        logger.info("  - 定标完成（%.1f s），开始执行 RPC 正射校正 (gdal.Warp 处理中，请耐心等待)...",
                    time.time() - calibration_started)                          # 记录定标耗时
        ortho_info = orthorectify_radiance(tmp_native, footprint, tmp_ortho, logger)  # 执行 RPC 正射
        logger.info("  - 正射完成，成果尺寸: %d x %d", ortho_info["width"], ortho_info["height"])  # 记录成果尺寸
        
        # ⑥ 原子替换成果文件：波长 → 正射 tif（先写 .tmp 再 os.replace）
        np.save(tmp_wl, wavelengths[keep])                            # 保存保留波段波长
        if args.keep_native: os.replace(tmp_native, paths["native"])  # 需要时保留未正射成果
        else: os.remove(tmp_native)                                   # 否则删除临时定标文件
        os.replace(tmp_ortho, paths["ortho"])                         # 正射成果转正
        os.replace(tmp_wl, paths["wavelengths"])                      # 波长文件转正
        
        # ⑦ 写 georef.json（格网/足迹/剔除信息）与 band_stats.csv
        logger.info("  - 正在写入几何与统计元数据文件...")                                    # 进度提示
        write_georef_json(paths["georef"], scene, tar_path, vn_tif, sw_tif, (n_cols, n_rows), ortho_info, n_output, n_original, footprint, rpc_correction, wavelengths, keep, removed, remove_groups, logger)  # 写几何 JSON
        write_band_stats_csv(stats, paths["stats"])                             # 写统计 CSV

        # ⑧ 质检图：真彩色 / 假彩色 / SWIR / 均值光谱
        if make_qc:                                                             # 需要质检图时
            logger.info("  - 正在生成质检图...")                                       # 进度提示
            removed_ranges = [(g["wavelength_nm"][0]/1000, g["wavelength_nm"][1]/1000) for g in remove_groups if g["bands_1based"]]  # 剔除波长区间（μm）
            qc_success = write_quality_figure(paths["ortho"], wavelengths[keep], removed_ranges, paths["qc"], scene, logger, band_stats=stats)  # 生成质检图
            if not qc_success:
                raise RuntimeError("质检图生成失败")  # 失败即报错

        logger.info("场景处理完成: %s（%d 波段，用时 %.1f s）",
                    scene, n_output, time.time() - started)                     # 记录完成与耗时
        return {"scene": scene, "status": "成功", "bands": n_output, "elapsed": time.time() - started, "message": "OK"}  # 返回成功结果
        
    except BaseException as e:
        for p in (tmp_native, tmp_ortho, tmp_wl):  # 清理所有临时文件
            try:
                if os.path.exists(p): os.remove(p)                              # 存在才删除
            except OSError: pass                                                # 删除失败忽略
        if created_dir and os.path.exists(scene_dir) and not os.listdir(scene_dir):  # 目录为本次新建且已空
            try: os.rmdir(scene_dir)                                            # 删除空目录
            except OSError: pass                                                # 删除失败忽略
        raise e                                                                 # 原样抛出，交给 main 记录
    finally:
        ds_vn = ds_swir = None  # 释放数据集


# =================================================================
# 批处理主流程
# =================================================================
def main(argv=None):
    """批处理入口：扫描压缩包 → （dry-run 仅打印计划）→ 逐景处理 → 汇总。"""
    args = parse_args(argv)                                           # 解析命令行参数
    out_root = os.path.abspath(args.out)                              # 输出根目录绝对路径
    os.makedirs(out_root, exist_ok=True)                              # 确保输出根目录存在
    log_path = args.log or os.path.join(out_root, "step1_batch.log")  # 日志文件路径
    logger = setup_logging(log_path, args.verbose)                    # 初始化日志
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)                     # 设置 GDAL 块缓存

    logger.info("=" * 60)                                # 打印分隔线
    logger.info("ZY1E AHSI 步骤1：辐射定标 + 坏波段剔除 + RPC 正射")   # 打印标题
    logger.info("源目录  : %s", os.path.abspath(args.src))  # 打印源目录
    logger.info("输出目录: %s", out_root)                    # 打印输出目录
    logger.info("日志文件: %s", log_path)                    # 打印日志文件
    logger.info("=" * 60)                                # 打印分隔线

    patterns = tuple(args.pattern) if args.pattern else TAR_PATTERNS  # 确定扫描模式
    try:
        scenes = discover_scenes(os.path.abspath(args.src), patterns)  # 扫描压缩包
    except Exception as e:
        logger.error("扫描失败: %s", e); return 2  # 扫描失败返回码 2
        
    if args.scene:  # 指定了景名过滤
        import fnmatch
        scenes = [p for p in scenes if any(fnmatch.fnmatch(scene_name(p), pat) for pat in args.scene)]  # 通配匹配景名
    if args.limit > 0: scenes = scenes[:args.limit]                             # 限制处理数量
    if not scenes: logger.error("未找到匹配压缩包"); return 2                           # 无匹配返回码 2

    if args.dry_run:                         # 仅打印计划
        for idx, p in enumerate(scenes, 1):  # 逐景打印
            scene = scene_name(p)            # 推导景名
            logger.info("[%d/%d] %s -> %s", idx, len(scenes),
                        os.path.basename(p),
                        os.path.join(out_root, scene, scene + OUT_ORTHO_SUFFIX))  # 打印预期输出路径
        logger.info("dry-run 完成：共 %d 景，未处理文件", len(scenes))                     # 结束提示
        return 0

    results = []                                          # 各景处理结果
    seen = {}                                             # 景名去重表
    started = time.time()                                 # 总计时起点
    for idx, tar_path in enumerate(scenes, start=1):      # 逐景处理
        scene = scene_name(tar_path)                      # 推导景名
        logger.info("="*50)                               # 打印分隔线
        logger.info("正在处理第 %d/%d 景数据", idx, len(scenes))  # 打印进度
        
        if scene in seen:                                # 景名重复
            logger.error("场景重名跳过: %s", scene); continue  # 跳过重复场景
        seen[scene] = tar_path                           # 登记已处理景名
        try:
            results.append(process_scene(tar_path, out_root, args, logger))  # 处理单景
        except Exception:
            logger.exception("处理失败: %s", os.path.basename(tar_path))            # 记录异常堆栈
            results.append({"scene": scene, "status": "失败", "bands": None, "elapsed": 0.0, "message": "见日志"})  # 计入失败

    success = sum(1 for r in results if r["status"] == "成功")                    # 统计成功景数
    skipped = sum(1 for r in results if r["status"].startswith("跳过"))           # 统计跳过景数
    failed = sum(1 for r in results if r["status"] == "失败")                     # 统计失败景数
    logger.info("="*50)                                                         # 打印分隔线
    logger.info("批处理结束: 成功 %d, 跳过 %d, 失败 %d，用时 %.1f s", success, skipped, failed, time.time() - started)  # 汇总结果
    for result in results:
        logger.info("  [%s] %s（%.1f s）", result["status"], result["scene"], result["elapsed"])  # 逐景结果
    if failed:
        logger.error("有 %d 景处理失败，详见日志: %s", failed, log_path)  # 有失败则提示
        return 1
    logger.info("全部完成，成果位于: %s", out_root)  # 全部完成
    return 0

if __name__ == "__main__":  # 脚本入口
    sys.exit(main())        # 以返回码退出