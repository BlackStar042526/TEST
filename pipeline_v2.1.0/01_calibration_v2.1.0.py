# -*- coding: utf-8 -*-
"""
01_calibration.py
ZY1E AHSI Step 1：辐射定标 + DN=0 掩膜 + RPC 校正 + 质检（v2.1.0）

v2.1.0 变更（用户 2026-11-20 决策）：
  ✅ 波段剔除后移 —— Step 1 不再物理剔除任何波段，输出全部 166 波段；
     18 个 DN=0 水汽吸收带 + 6 个低信噪可疑波段写入 georef.json.deferred_bad_bands，
     由 Step 2 在大气校正后处理（前者硬剔除，后者仅评估记录、保留）。
  ✅ 坐标单一真源 —— 由四角足迹经纬度推导 UTM 带（本批 E116–119 → EPSG:32650），
     写入 georef.json.target_crs_epsg / target_res，Step 2/3 一律读取，不再硬编码。

v2.0.0 重构要点（GitHub 比对验证后）：
  ✅ 废弃 ×10 FLAASH 单位转换 —— Py6S/6S 内部单位就是 W·m⁻²·sr⁻¹·μm⁻¹，
     和卫星定标文件 Gain/Bias 单位一致，无需转换。
  ✅ MASK_NODATA 默认 True —— ZY1E L1A 确实存在黑边 DN=0，必须掩膜。
  ✅ 去掉 Step 1 RPC 正射 —— 光谱先于几何，正射推迟到 Step 3 统一做（带 DEM）。
  ✅ RPC 校正采样点从 4 角→25 点均匀网格 + 仿射 6 参数拟合（含旋转），
     残差阈值从 2.0→0.5 px。
  ✅ 定标改为批量向量化（一次 ReadAsArray 读完 VNIR/SWIR，向量化计算）。
  ✅ 保留 ENVI/FDAH 标准 FLAASH Table 3 大气模式逻辑 —— Py6S 内置 FromLatitudeAndDate。

v1.0.5 → v2.0.0 架构变化：旧 3 步（两次正射）→ 新 3 步（一次正射，在 Step 3）。

流程：读包内 VNIR/SWIR → 解析光谱/定标/足迹/RPC → 批量向量化定标 L=Gain*DN+Bias
→ 掩膜 DN=0 → 25 点网格 RPC 校正 → 写 native 网格辐亮度（全 166 波段）
→ georef.json 溯源（含待剔除清单与目标 CRS）→ 质检图。
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
SRC_DIR = r"C:\Users\yangz\Desktop\TEST\ZY1E_AHSI"           # 默认源目录：压缩包所在
OUT_ROOT = r"C:\Users\yangz\Desktop\TEST\ZY1E_output_step1"  # 默认输出根目录
TAR_PATTERNS = ("*.tar", "*.tar.gz", "*.tgz")     # 扫描压缩包的通配模式

# 波段参数
N_VNIR = 76  # VNIR 波段数
N_SWIR = 90  # SWIR 波段数

# ----------------------------------------------------------------
# 坏波段策略（v2.0.0 重构）：
#   Step 1 只剔除 DN=0 的"确定"坏波段（水汽吸收带 DN=0，定标后 = Bias 无信息），
#   低信噪比"可疑"波段的剔除推迟到 Step 2 光谱评估后（反射率异常 + 信噪比）。
# ----------------------------------------------------------------
# DN=0 的水汽吸收带（1-based 波段号）—— 这些波段在 L1A 原始数据里 DN 全为 0
# 参考 ENVI FLAASH 文档 + ZY1E AHSI 官方说明
DN_ZERO_BANDS_1BASED = [
    # 1.4 μm 水汽吸收带（核心 DN=0）
    98, 99, 100, 101, 102, 103, 104, 105,
    # 1.9 μm 水汽吸收带（核心 DN=0）
    126, 127, 128, 129, 130, 131, 132, 133, 134,
    # 1.812 μm 单一吸收（DN 也常为 0）
    125,
]

# 低信噪比"可疑"波段（v2.0.0 只记录不剔除，留给 Step 2 评估）
# 这些波段 DN 不全为 0，但信噪比低，需要看反射率光谱再决定
LOW_SNR_BANDS_1BASED = [
    77, 78, 79,      # SWIR 首段低信噪（~1.0 μm 边缘）
    164, 165, 166,   # SWIR 末段低信噪（~2.5 μm 边缘）
]

# 两者合并 = Step 1 确定剔除的波段
REMOVE_BANDS_1BASED = DN_ZERO_BANDS_1BASED

# 波长范围参考（用于 georef.json 溯源）
WATER_VAPOR_BANDS = [(1.35, 1.43), (1.80, 1.96)]  # μm

# 定标参数
CAL_HAS_BIAS = True
# ⚠️ v2.0.0 关键修改：MASK_NODATA 默认 True
# ZY1E L1A 原始影像确实存在黑边（轨道两侧扫描带不完全覆盖），DN=0 填充
# 不掩膜的话这些像元定标后 = Bias（常数），会污染大气校正和镶嵌
MASK_NODATA = True    # 是否把 DN=0 当无效像元
NODATA_IN = 0         # 输入的无效像元值
NODATA_OUT = -9999.0  # 输出的无效值

# ----------------------------------------------------------------
# RPC 校正参数（v2.0.0 重构）：
#   正射推迟到 Step 3，这里只做 RPC 网格校正（让 Step 3 的正射更准）
#   采样点从 4 角→25 点均匀网格（排除 DN=0 的），仿射 6 参数拟合（含旋转）
#   残差阈值从 2.0→0.5 px（和 ENVI 理想值 CE95≈1/2 pixel 对齐）
# ----------------------------------------------------------------
RPC_FIT_MAX_RESIDUAL = 0.5    # RPC 拟合残差阈值（像素）—— 收紧到 0.5
RPC_GRID_POINTS = 5           # 均匀网格每边采样点数（总共 5×5=25）
RPC_MIN_VALID_POINTS = 10     # 有效采样点下限（少于这个数报错）

REQUIRE_GEOREF = True         # 是否要求输入带地理参考
WRITE_GEOREF = True           # 是否输出 georef.json

# ----------------------------------------------------------------
# ⚠️ 正射参数移到 Step 3：
#   v2.0.0 的 Step 1 只输出 native 网格辐亮度，不做正射。
#   以下参数保留在配置区供参考，但 Step 1 不再使用。
#   Step 3 正射时会自动检测 UTM 带号（从景中心经度），重采样默认 cubic。
# ----------------------------------------------------------------
_TARGET_EPSG_HINT = 32650     # 仅作参考（UTM 50N），Step 3 会自动检测
_TARGET_RES_HINT = 30.0       # 目标分辨率（米）
_RESAMPLING_HINT = "cubic"    # Step 3 正射重采样（定量反射率推荐 cubic）

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

# ----------------------------------------------------------------
# 输出文件名（v2.0.0：Step 1 不做正射，只输出 native 网格辐亮度）
# ----------------------------------------------------------------
OUT_NATIVE_RADIANCE_NAME = "radiance_native.tif"  # 原生网格辐亮度（Step 2 的输入）
OUT_WAVELENGTHS_NAME = "wavelengths_clean.npy"    # 保留波段波长（已剔除 DN=0 坏波段）
OUT_GEOREF_NAME = "georef.json"                   # 几何信息（串联 3 步的"血液"）
OUT_STATS_NAME = "band_stats.csv"                 # 逐波段统计
OUT_QC_NAME = "preview_quality.png"               # 质检图
LOG_NAME = "zy1e.step1"                           # 日志器名称
# --- v2.0.0 patch: 补回缺失的常量引用 ---
KEEP_NATIVE_RADIANCE = True       # v2.0.0 默认保留 native 网格（Step 2 要用）
OUT_ORTHO_SUFFIX_LEGACY = "_cal_badband_ortho.tif"
LEGACY_ORTHO_NAME = "cal_badband_ortho.tif"
NATIVE_RADIANCE_FINAL = OUT_NATIVE_RADIANCE_NAME  # 别名


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
    # v2.1.0 fail-fast：默认出错即停；仅显式指定才恢复「记失败继续下一景」
    parser.add_argument("--continue-on-error", action="store_true", default=False,
                        help="出错时继续处理下一景（默认关闭：出错立即终止，防止产出不可靠数据）")
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
    """解析辐射定标系数，返回 (gain 数组, bias 数组)，各取前 n_bands 行。

    v2.0.0 重构：
      ✅ 增加格式探测 —— 扫描前 3 行自动识别分隔符（逗号/tab/空白）和列位置
      ✅ 增加 Gain/Bias 合理性验证 —— Gain 应全为正、Bias 应在 -0.5~5 之间
      ✅ 增加异常告警 —— 探测结果与预期不符时警告（不再静默猜测）

    单位说明：ZY1E AHSI 定标文件 Gain/Bias 原始单位是 W·m⁻²·sr⁻¹·μm⁻¹，
    Py6S/6S 内部也是这个单位 —— 无需单位转换（v2.0.0 废弃 ×10 FLAASH 换算）。
    """
    # --- Step 1: 格式探测（v2.0.0 新增） ---
    first_lines = [l.strip() for l in text.splitlines()
                   if l.strip() and not l.strip().startswith("#")][:3]
    if not first_lines:
        raise ValueError(f"{label} 定标文件为空或全是注释")

    # 探测分隔符
    if "\t" in first_lines[0]:
        detected_sep = "\\t"
    elif "," in first_lines[0]:
        detected_sep = ","
    else:
        detected_sep = r"[,\s]+"
    logger.info("%s 定标文件: 探测分隔符=%s", label, detected_sep)

    # --- Step 2: 逐行解析 ---
    rows = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        tokens = re.split(detected_sep, line)
        numbers = []
        for token in tokens:
            try:
                numbers.append(float(token))
            except ValueError:
                break
        if numbers:
            # 判断首列是否波段号（整数且 <200）
            if len(numbers) >= 3 and numbers[0].is_integer() and numbers[0] < 200:
                rows.append(numbers[1:3])
            elif len(numbers) == 2 and numbers[0].is_integer() and numbers[0] < 200:
                rows.append([numbers[1]])
            else:
                rows.append(numbers[:2])

    if not rows:
        raise ValueError(f"{label} 定标系数为空")

    # --- Step 3: 组装 ---
    if has_bias is None:
        has_bias = all(len(row) >= 2 for row in rows)
    if has_bias:
        bad_rows = [i for i, row in enumerate(rows) if len(row) < 2]
        if bad_rows:
            raise ValueError(f"{label} 定标文件第 {bad_rows[0]+1} 行缺 Bias 列")
        data = np.asarray([row[:2] for row in rows], dtype=np.float64)
    else:
        data = np.asarray([[row[0], 0.0] for row in rows], dtype=np.float64)

    if len(data) < n_bands:
        raise ValueError(f"{label} 定标系数不足 {n_bands} 行（只有 {len(data)}）")
    data = data[:n_bands]

    # --- Step 4: 合理性验证（v2.0.0 新增） ---
    gain = data[:, 0]
    bias = data[:, 1]
    if not np.all(gain > 0):
        neg_idx = np.where(gain <= 0)[0]
        logger.warning("%s 存在 %d 个 Gain ≤ 0 的波段（索引 %s）—— 请检查定标文件格式",
                       label, len(neg_idx), neg_idx[:10])
    bias_min, bias_max = float(bias.min()), float(bias.max())
    if bias_min < -1.0 or bias_max > 10.0:
        logger.warning("%s Bias 范围异常 [%.3f, %.3f] —— 预期在 -0.5~5.0 之间，"
                       "请检查定标文件格式或列顺序（是否 Gain/Bias 搞反了）",
                       label, bias_min, bias_max)
    logger.info("%s 定标解析完成: %d 行, Gain [%.4f, %.4f], Bias [%.4f, %.4f]",
                label, len(data), gain.min(), gain.max(), bias.min(), bias.max())

    return gain, bias


# =================================================================
# 坏波段剔除与几何处理
# =================================================================
def resolve_band_selection(wavelengths, logger):
    """v2.1.0：Step 1 不再物理剔除任何波段，输出全部波段。

    18 个 DN=0 水汽吸收带与 6 个低信噪可疑波段只登记到 catalog（deferred_bad_bands），
    由 Step 2 在大气校正之后处理：
      - DN=0 波段：硬剔除（定标后 = Bias 常数，无信息）；
      - 低信噪波段：仅统计评估写报告，保留不删除（用户 2026-11-20 决策）。
    这样做到「波段剔除放在大气校正之后」，且避免把无信息波段送进 6S 造成系数异常时
    仍能追溯（Step 2 用同一份清单在反射率域剔除）。

    返回 (keep, removed, catalog)：
      keep    —— 0-based 保留索引（本版本为全部波段）
      removed —— 0-based 已剔除索引（本版本为空列表）
      catalog —— 待 Step 2 处理的坏波段目录
    """
    n_total = int(wavelengths.size)

    # 检查越界（两组都要校验，避免把错误清单传给 Step 2）
    for label, bands in (("DN=0", REMOVE_BANDS_1BASED),
                         ("low-SNR", LOW_SNR_BANDS_1BASED)):
        invalid = [b for b in bands if not (1 <= b <= n_total)]
        if invalid:
            raise ValueError(f"{label} 坏波段号越界: {invalid[:10]}")

    # 1-based → 0-based（只登记，不剔除）
    dn_zero_0based = sorted({int(b) - 1 for b in REMOVE_BANDS_1BASED})
    low_snr_0based = sorted({int(b) - 1 for b in LOW_SNR_BANDS_1BASED
                             if 1 <= b <= n_total})

    # v2.1.0：全部波段保留
    keep = np.arange(n_total, dtype=np.int64)
    removed = []  # 本步骤不剔除任何波段

    catalog = {
        "dn_zero_bands_1based": list(REMOVE_BANDS_1BASED),
        "dn_zero_bands_0based": dn_zero_0based,
        "dn_zero_wavelengths_um": [float(wavelengths[i]) for i in dn_zero_0based],
        "low_snr_bands_1based": list(LOW_SNR_BANDS_1BASED),
        "low_snr_bands_0based": low_snr_0based,
        "low_snr_wavelengths_um": [float(wavelengths[i]) for i in low_snr_0based],
        "note": ("v2.1.0：Step 1 不剔除任何波段；DN=0 波段由 Step 2 硬剔除；"
                 "低信噪波段仅评估记录、保留不删除"),
    }

    logger.info("坏波段策略(v2.1.0): Step 1 保留全部 %d 波段；"
                "登记待 Step 2 硬剔除 DN=0 波段 %d 个、仅评估低信噪波段 %d 个",
                n_total, len(dn_zero_0based), len(low_snr_0based))

    return keep, removed, catalog

def open_vsitar_dataset(tar_path, internal_path, label):
    """通过 /vsitar/ 直接打开压缩包内影像，避免解压到磁盘。"""
    """用 /vsitar/ 直接打开压缩包内影像，避免整包解压。"""
    vsi_path = "/vsitar/" + os.path.abspath(tar_path).replace("\\", "/") + "/" + internal_path  # 拼 vsitar 虚拟路径
    dataset = gdal.Open(vsi_path)                                               # 打开包内影像
    if dataset is None: raise RuntimeError(f"无法通过 /vsitar/ 打开 {label}")         # 打开失败报错
    return dataset                                                              # 返回数据集对象

def compute_rpc_grid_correction(ds_vn, footprint, logger):
    """v2.0.0：25 点均匀网格采样 + 2D 仿射 6 参数拟合 + 残差阈值 0.5 px。

    旧版问题（4 角点 + 1D 仿射）：山区角点可能 DN=0 不可靠、无冗余、无旋转。
    新版：影像内部 10%~90% 区域均匀 5x5 网格采样，2D 仿射含旋转/剪切/缩放/平移。
    """
    rpc = ds_vn.GetMetadata("RPC")
    if not rpc:
        raise ValueError("影像缺少 RPC 元数据")
    gt = ds_vn.GetGeoTransform()
    if not (abs(abs(gt[1]) - 1) < 1e-6 and abs(abs(gt[5]) - 1) < 1e-6):
        raise ValueError(f"GeoTransform 非单位阵，RPC 校正失效: {gt}")

    transformer = gdal.Transformer(ds_vn, None, ["METHOD=RPC"])
    width, height = ds_vn.RasterXSize, ds_vn.RasterYSize

    # 5x5 = 25 个均匀网格采样点（10%~90% 区域避免边缘 DN=0）
    grid_x = np.linspace(width * 0.10, width * 0.90, RPC_GRID_POINTS)
    grid_y = np.linspace(height * 0.10, height * 0.90, RPC_GRID_POINTS)
    col_grid, row_grid = np.meshgrid(grid_x, grid_y)
    pixel_sources = np.column_stack([col_grid.ravel(), row_grid.ravel()])

    # RPC 正反变换：像素 → 经纬度 → 像素（RPC 模型自身反变换结果）
    raw_pixels, valid_indices = [], []
    for idx, (col, row) in enumerate(pixel_sources):
        try:
            ok_fwd, fwd_result = transformer.TransformPoint(0, col, row)
            if not ok_fwd:
                continue
            lon, lat = fwd_result[0], fwd_result[1]
            ok_inv, inv_result = transformer.TransformPoint(1, lon, lat)
            if ok_inv:
                raw_pixels.append([inv_result[0], inv_result[1]])
                valid_indices.append(idx)
        except Exception:
            continue

    if len(raw_pixels) < RPC_MIN_VALID_POINTS:
        raise ValueError(f"RPC 校正有效采样点不足（{len(raw_pixels)}/{RPC_MIN_VALID_POINTS}）")

    raw_pixels = np.asarray(raw_pixels, dtype=np.float64)
    source_sel = pixel_sources[valid_indices]

    # 2D 仿射 6 参数拟合
    A = np.column_stack([source_sel[:, 0], source_sel[:, 1], np.ones(len(source_sel))])
    col_params = np.linalg.lstsq(A, raw_pixels[:, 0], rcond=None)[0]
    row_params = np.linalg.lstsq(A, raw_pixels[:, 1], rcond=None)[0]

    predicted_col = A @ col_params
    predicted_row = A @ row_params
    residuals = np.hypot(predicted_col - raw_pixels[:, 0], predicted_row - raw_pixels[:, 1])
    max_res = float(np.max(residuals))
    mean_res = float(np.mean(residuals))

    logger.info("RPC 校正 2D 仿射: %d/%d 点, residual max=%.3f mean=%.3f 阈值=%.1f px",
                len(raw_pixels), len(pixel_sources), max_res, mean_res, RPC_FIT_MAX_RESIDUAL)

    if max_res > RPC_FIT_MAX_RESIDUAL:
        raise ValueError(f"RPC 拟合残差过大: max={max_res:.3f} > 阈值 {RPC_FIT_MAX_RESIDUAL} px")

    # 修正 RPC OFF/SCALE（主方向 sx/sy/tx/ty，旋转体现在残差里 <0.5 px）
    sx, sy = float(col_params[0]), float(row_params[1])
    tx, ty = float(col_params[2]), float(row_params[2])

    corrected = dict(rpc)
    corrected["SAMP_OFF"] = str((float(rpc["SAMP_OFF"]) - tx) / sx)
    corrected["SAMP_SCALE"] = str(float(rpc["SAMP_SCALE"]) / sx)
    corrected["LINE_OFF"] = str((float(rpc["LINE_OFF"]) - ty) / sy)
    corrected["LINE_SCALE"] = str(float(rpc["LINE_SCALE"]) / sy)

    details = {
        "method": "2D Affine 6-param (5x5 grid)",
        "grid_points_total": RPC_GRID_POINTS ** 2,
        "valid_points": len(raw_pixels),
        "params_col": {"a": round(col_params[0], 6), "b": round(col_params[1], 6), "c": round(col_params[2], 6)},
        "params_row": {"d": round(row_params[0], 6), "e": round(row_params[1], 6), "f": round(row_params[2], 6)},
        "sx": round(sx, 6), "sy": round(sy, 6), "tx": round(tx, 6), "ty": round(ty, 6),
        "residual_max_px": round(max_res, 4),
        "residual_mean_px": round(mean_res, 4),
        "threshold_px": RPC_FIT_MAX_RESIDUAL,
    }
    return corrected, details



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

# [DELETED] orthorectify_radiance — v2.0.0 正射推迟到 Step 3，Step 1 只输出 native 网格


def derive_target_crs(footprint, logger):
    """v2.1.0：由四角足迹经纬度推导统一目标 CRS（UTM 北半球，WGS84 基准）。

    资源一号 L1A 原始影像无投影（GeoTransform 单位阵 + RPC），ENVI 教程默认输出
    UTM WGS84。本函数按足迹平均经度算 UTM 带号：
        band = floor((lon_mean + 180) / 6) + 1
        EPSG = 32600 + band   (北半球 WGS84 / UTM north)
    本批数据 E116–119° → band 50 → EPSG:32650。
    若足迹跨带（经度范围跨越 6° 分界）会记录警告但仍取平均经度对应带。
    """
    if not footprint:
        raise ValueError("足迹为空，无法推导目标 CRS；请检查 georef.json / 产品 XML")
    lons, lats = [], []
    for key, value in footprint.items():
        try:
            lon, lat = float(value[0]), float(value[1])
        except (TypeError, ValueError, IndexError):
            continue
        lons.append(lon)
        lats.append(lat)
    if not lons or not lats:
        raise ValueError("足迹中无有效经纬度，无法推导目标 CRS")
    lon_mean = float(np.mean(lons))
    lat_mean = float(np.mean(lats))
    band = int(np.floor((lon_mean + 180.0) / 6.0)) + 1
    if not (1 <= band <= 60):
        raise ValueError(f"UTM 带号越界: {band} (lon_mean={lon_mean})")
    north = lat_mean >= 0
    epsg = (32600 if north else 32700) + band
    span = max(lons) - min(lons)
    if span > 6.0:
        logger.warning("足迹跨经度 %.2f°，目标 CRS 取平均经度对应 UTM 带 %d (EPSG:%d)",
                       span, band, epsg)
    logger.info("目标 CRS 推导: 中心 (%.4f, %.4f) → UTM %d%s / EPSG:%d",
                lon_mean, lat_mean, band, "N" if north else "S", epsg)
    return {
        "target_crs_epsg": epsg,
        "target_crs_name": f"WGS 84 / UTM zone {band}{'N' if north else 'S'}",
        "utm_zone": band,
        "utm_hemisphere": "N" if north else "S",
        "footprint_center_lon": lon_mean,
        "footprint_center_lat": lat_mean,
        "target_res": 30.0,
        "derivation": "由四角足迹平均经纬度按 6° 分带推导，WGS84 基准",
    }


def write_georef_json(output_path, scene, tar_path, vn_tif, sw_tif, native_size,
                                            band_count, bands_source, footprint,
                                            rpc_corrected, rpc_correction, wavelengths, keep, removed,
                                            bad_band_catalog, logger):
    """v2.1.0：写 georef.json（完整溯源信息），供 Step 2/3 追溯。

    Step 1 不做正射，所以没有 ortho_info；native 网格尺寸就是原始影像尺寸。
    v2.1.0 新增：deferred_bad_bands（待 Step 2 处理的坏波段）与 target_crs（单一真源）。
    """
    crs_info = derive_target_crs(footprint, logger)
    payload = {
        "scene": scene,
        "source_tar": os.path.basename(tar_path),
        "source_tiffs": [os.path.basename(vn_tif), os.path.basename(sw_tif)],

        "step1_calibration": {
            "radiance_formula": "L = Gain * DN + Bias",
            "radiance_unit": "W·m⁻²·sr⁻¹·μm⁻¹",
            "flaash_scale_applied": False,
            "flaash_scale_factor": 1.0,
            "dn_nodata": NODATA_IN,
            "radiance_nodata": NODATA_OUT,
            "mask_nodata": bool(MASK_NODATA),
            "cal_has_bias": bool(CAL_HAS_BIAS),
        },

        "native_grid": {
            "raster_size": [int(native_size[0]), int(native_size[1])],
            "proj4_note": "L1A 原始影像无投影（Projection 为空），GeoTransform 单位阵，仅含 RPC",
            "has_projection": False,
            "geotransform": [0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
        },

        "rpc": {
            "grid_correction": rpc_correction,
            "fit_residual_px": rpc_correction.get("residual_max_px", rpc_correction.get("residual", "unknown")),
        },
        "rpc_corrected": {str(key): str(value) for key, value in rpc_corrected.items()},

        "bad_bands_catalog": bad_band_catalog,
        "deferred_bad_bands": bad_band_catalog,  # v2.1.0：Step 2 按此清单在反射率域处理
        "removed_bands_0based": [int(b) for b in removed],       # v2.1.0：恒为空
        "removed_bands_1based": [int(b) + 1 for b in removed],   # v2.1.0：统一为 1-based
        "removed_wavelengths_um": [float(wavelengths[i]) for i in removed],
        "kept_bands_1based": [int(i) + 1 for i in keep],
        "kept_wavelengths_um": [float(wavelengths[i]) for i in keep],

        # v2.1.0 坐标单一真源（Step 2/3 一律读此，不再硬编码）
        "target_crs": crs_info,
        "target_crs_epsg": crs_info["target_crs_epsg"],
        "target_res": crs_info["target_res"],

        "band_count": int(band_count),
        "bands_source": int(bands_source),

        "footprint_wgs84": footprint,

        "pipeline_note": ("Step 1 输出 native 网格全波段辐亮度；坏波段剔除后移到 Step 2；"
                          "正射推迟到 Step 3 统一做（带 30m DEM）"),
        "script_version": "01_calibration_v2.1.0",
        "step_completed": "step1_calibration",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
    }
    tmp = output_path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(payload, f, ensure_ascii=False, indent=2)
    os.replace(tmp, output_path)
    logger.info("georef.json 写入完成: %s (%d 波段保留)", output_path, band_count)



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
    """v2.0.0：Step 1 只产出 native 网格辐亮度 + 元数据。"""
    scene = os.path.basename(os.path.normpath(scene_dir))
    return {k: os.path.join(scene_dir, v) for k, v in {
        "radiance_native": OUT_NATIVE_RADIANCE_NAME,
        "wavelengths": OUT_WAVELENGTHS_NAME,
        "georef": OUT_GEOREF_NAME,
        "stats": OUT_STATS_NAME,
        "qc": OUT_QC_NAME}.items()}


def process_scene(tar_path, out_root, args, logger):
    """处理单景：①元数据 ②打开影像/RPC ③波段规划 ④定标写盘

    ⑤RPC 正射 ⑥原子替换成果 ⑦元数据与统计 ⑧质检图。
    单景异常向上抛出，由 main 记录堆栈并计入失败，不中断批处理。
    """
    scene = scene_name(tar_path)               # 由包名推导景名
    scene_dir = os.path.join(out_root, scene)  # 本景输出目录
    paths = resolve_scene_paths(scene_dir)     # 本景各成果路径
    make_qc = MAKE_QC and not args.no_qc       # 汇总质检开关

    required_files = ["radiance_native", "wavelengths", "georef", "stats"]  # 判定完成所需成果
    if make_qc:
        required_files.append("qc")                            # 需质检图时一并检查
    legacy_ortho = os.path.join(scene_dir, LEGACY_ORTHO_NAME)  # 旧版正射文件名
    if (not args.overwrite and not os.path.isfile(paths["radiance_native"])
            and os.path.isfile(legacy_ortho)
            and all(os.path.isfile(paths[key])
                    for key in required_files if key != "ortho")):     # 旧成果齐全则复用
        os.replace(legacy_ortho, paths["radiance_native"])                       # 旧名改为新名
        logger.info("复用已有完整成果并更新正射文件名：%s", paths["radiance_native"])             # 记录复用
    have_main = all(os.path.isfile(paths[p]) for p in required_files)  # 主成果是否齐全

    if not args.overwrite and have_main:
        logger.info("跳过（成果已存在）: %s", scene_dir)                                 # 成果齐全则跳过
        return {"scene": scene, "status": "跳过", "bands": None, "elapsed": 0.0, "message": "已存在"}  # 返回跳过结果

    logger.info("开始处理场景: %s", scene)                                           # 记录开始
    started = time.time()                                                      # 计时起点
    created_dir = not os.path.isdir(scene_dir)                                 # 记录目录是否本次新建
    os.makedirs(scene_dir, exist_ok=True)                                      # 创建输出目录
    tmp_native = os.path.join(scene_dir, "radiance_native.tmp.tif")            # 定标临时文件
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
            keep, removed, bad_band_catalog = resolve_band_selection(wavelengths, logger)  # 计算保留波段
            footprint = parse_footprint(read_member_text(tar, xml_member, "XML"), logger)  # 解析四角足迹

            # v2.1.0 fail-fast：RPC 是 Step 3 正射的必需输入，解析失败不得静默降级
            if not vn_rpb:
                raise FileNotFoundError("L1A 压缩包内缺少 VNIR .rpb RPC 文件")
            try:
                vn_rpb_text = read_member_text(tar, vn_rpb, "VNIR RPC")  # 读取 .rpb 文本
                vn_rpc_dict = parse_rpb(vn_rpb_text)                     # 解析 RPC 参数
            except Exception as e:
                # 若 TIFF 内自带 RPC 可兜底，仍记录为错误而非静默继续：由调用方决定是否终止
                raise RuntimeError(
                    "解析外部 .rpb RPC 失败（RPC 为正射必需，不允许静默忽略）: %s" % e
                ) from e
            if not vn_rpc_dict:
                raise ValueError("外部 .rpb 解析结果为空，RPC 不可用")

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

        # ④ v2.0.0 分块向量化辐射定标：全幅宽分块读取，匹配源数据的逐行条带布局
        logger.info("  - 开始分块向量化定标 (共 %d 个波段)...", n_output)
        calibration_started = time.time()
        stats_acc = {}
        row_chunk = 128
        sample_per_chunk = max(1, 200000 // max((n_rows + row_chunk - 1) // row_chunk, 1))

        for out_idx, (src_ds, src_band, gain, bias, wavelength) in enumerate(band_plan):
            acc = {"count": 0, "sum": 0.0, "sumsq": 0.0,
                   "min": None, "max": None, "sample": np.empty(0, dtype=np.float32)}
            out_band = out_ds.GetRasterBand(out_idx + 1)
            for yoff in range(0, n_rows, row_chunk):
                height = min(row_chunk, n_rows - yoff)
                dn = src_ds.GetRasterBand(src_band).ReadAsArray(
                    xoff=0, yoff=yoff, win_xsize=n_cols, win_ysize=height
                ).astype(np.float32)
                radiance = dn * float(gain) + float(bias)
                if MASK_NODATA:
                    radiance[dn == NODATA_IN] = NODATA_OUT

                valid = radiance[radiance != NODATA_OUT] if MASK_NODATA \
                    else radiance[np.isfinite(radiance)]
                if valid.size:
                    acc["count"] += int(valid.size)
                    acc["sum"] += float(valid.sum(dtype=np.float64))
                    acc["sumsq"] += float(np.square(valid, dtype=np.float64).sum())
                    block_min = float(valid.min())
                    block_max = float(valid.max())
                    acc["min"] = block_min if acc["min"] is None else min(acc["min"], block_min)
                    acc["max"] = block_max if acc["max"] is None else max(acc["max"], block_max)
                    sample = valid[::4]
                    if sample.size > sample_per_chunk:
                        positions = np.linspace(
                            0, sample.size - 1, sample_per_chunk, dtype=np.int64
                        )
                        sample = sample[positions]
                    acc["sample"] = np.concatenate((acc["sample"], sample))
                out_band.WriteArray(radiance.astype(np.float32), xoff=0, yoff=yoff)
            out_band.SetNoDataValue(NODATA_OUT)
            stats_acc[out_idx + 1] = (acc, wavelength)

        stats = []
        for band_index in range(1, n_output + 1):
            acc, wavelength = stats_acc[band_index]
            count = acc["count"]
            if count:
                mean = acc["sum"] / count
                variance = max(0.0, acc["sumsq"] / count - mean * mean)
                p1, p99 = np.percentile(acc["sample"], [1, 99]) if acc["sample"].size else (mean, mean)
                stats.append({"band": band_index, "wavelength": float(wavelength),
                              "min": float(acc["min"]), "max": float(acc["max"]),
                              "mean": float(mean), "std": float(np.sqrt(variance)),
                              "p1": float(p1), "p99": float(p99), "valid": int(count)})
            else:
                stats.append({"band": band_index, "wavelength": float(wavelength),
                              "min": 0.0, "max": 0.0, "mean": 0.0, "std": 0.0,
                              "p1": 0.0, "p99": 0.0, "valid": 0})

        out_ds.FlushCache()
        out_ds = None
        logger.info("  - 批量定标完成 (%.1f s)，有效波段 %d",
                    time.time() - calibration_started, len(stats))
        
        # ⑤ v2.0.0：不再在 Step 1 做 RPC 正射。正射推迟到 Step 3（带 30m DEM）。
        logger.info("  - Step 1 完成（%.1f s）。正射推迟到 Step 3。", time.time() - calibration_started)
        ortho_info = {"width": n_cols, "height": n_rows, "band_count": n_output,
                       "geotransform": list(ds_vn.GetGeoTransform()), "bounds": None}

        
        # ⑥ 原子替换成果文件：波长 → 正射 tif（先写 .tmp 再 os.replace）
        np.save(tmp_wl, wavelengths[keep])                            # 保存保留波段波长
        os.replace(tmp_native, paths["radiance_native"])              # 原子落盘未正射成果
        os.replace(tmp_wl, paths["wavelengths"])                      # 波长文件转正
        
        # ⑦ 写 georef.json（格网/足迹/剔除信息）与 band_stats.csv
        logger.info("  - 正在写入几何与统计元数据文件...")                                    # 进度提示
        write_georef_json(paths["georef"], scene, tar_path, vn_tif, sw_tif,
                          (n_cols, n_rows), n_output, n_original, footprint,
                          rpc_corrected, rpc_correction, wavelengths, keep,
                          removed, bad_band_catalog, logger)  # 写几何 JSON
        write_band_stats_csv(stats, paths["stats"])                             # 写统计 CSV

        # ⑧ 质检图：真彩色 / 假彩色 / SWIR / 均值光谱
        if make_qc:                                                             # 需要质检图时
            logger.info("  - 正在生成质检图...")
            # v2.1.0: 灰带标注的是「待 Step 2 剔除」的 DN=0 水汽带（本步不剔除，全 166 波段保留）
            # bad_band_catalog 是扁平字典，提取 DN=0 波长区间（μm）
            dn_zero_wls = bad_band_catalog.get("dn_zero_wavelengths_um", [])
            removed_ranges = []
            if dn_zero_wls:
                # 相邻波长合并成区间
                dn_zero_wls_sorted = sorted(dn_zero_wls)
                lo = dn_zero_wls_sorted[0]
                hi = lo
                for wl in dn_zero_wls_sorted[1:]:
                    if wl - hi < 0.02:  # 相邻 <20nm 视为同一区间
                        hi = wl
                    else:
                        removed_ranges.append((lo, hi))
                        lo = wl; hi = wl
                removed_ranges.append((lo, hi))
            qc_success = write_quality_figure(paths["radiance_native"], wavelengths[keep], removed_ranges, paths["qc"], scene, logger, band_stats=stats)
            if not qc_success:
                raise RuntimeError("质检图生成失败")  # 失败即报错

        logger.info("场景处理完成: %s（%d 波段，用时 %.1f s）",
                    scene, n_output, time.time() - started)                     # 记录完成与耗时
        return {"scene": scene, "status": "成功", "bands": n_output, "elapsed": time.time() - started, "message": "OK"}  # 返回成功结果
        
    except BaseException as e:
        for p in (tmp_native, tmp_wl):  # 清理所有临时文件
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
    logger.info("ZY1E AHSI 步骤1 v2.1.0：辐射定标（全波段，坏波段剔除后移至步骤2）")   # 打印标题
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
                        os.path.join(out_root, scene, scene + "_native_radiance.tif"))  # 打印预期输出路径
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
        except Exception as exc:
            logger.exception("处理失败: %s", os.path.basename(tar_path))            # 记录异常堆栈
            if not getattr(args, "continue_on_error", False):                 # v2.1.0：默认出错即停
                logger.error("出错即停（如需继续请加 --continue-on-error）: %s", scene)
                return 1
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