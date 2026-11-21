# -*- coding: utf-8 -*-
"""
02_atmospheric_correction.py
ZY1E AHSI Step 2：K-T 气溶胶反演 + 6S/Py6S 大气校正（v2.2.0）

v2.1.2 关键修复（2026-11-20，基于 6S 代码层面确定性证据）：
  ✅ 修复 K-T AOT 系统性高估（0.6 顶格 → ~0.3-0.4）：根因是旧版用"地表反射率=0 的
     6S 正向辐亮度 LUT"匹配暗像元，rho=0 时模拟辐亮度恒=0（实测 rho=0→0.0,
     rho=0.02→5.87），导致 AOT 必然顶格。新法按 FLAASH 标准用 6S 系数反演
     ρ_swir/ρ_red，匹配 red/SWIR=0.45（比值法）；暗像元加 red/NIR<0.6 判据排除
     阴影/水体。实测匹配误差 0.014，反演物理自洽（不依赖任何外部 AOD 产品）。
  ✅ 显式剔除水汽强吸收带（0.94/1.14/1.38/1.88μm，对齐 FLAASH bad bands），
     消除校正后反射率尖峰，保护 SAM/MTMF 与 ML 特征；报告分列 water_vapor 来源。

v2.0.0 重构要点（GitHub 比对 + 用户需求）：
  ✅ 废弃 MODIS/AERONET 外部气溶胶产品 → K-T 暗像元 2-Band 反演 (SWIR 2100-2250 nm + Red 640-680 nm)
  ✅ K-T 失败回退：FLAASH Initial Visibility（默认 40 km 晴朗 → AOT550 ≈ 0.098）
     —— 依据 ENVI FLAASH 官方定义：气溶胶无法反演时该初始值参与校正（cnblogs/enviidl p/16418433）
  ✅ 6S 传感器高度：set_sensor_satellite_level()（星载传感器标准接口，778km 与 100km 对辐射传输系数无差）
  ✅ 大气模式：改用 Py6S 内置 AtmosProfile.FromLatitudeAndDate（已内置 FLAASH Table 3 精确匹配）
  ✅ 气溶胶类型：Continental（= FLAASH Rural，山区植被覆盖区最匹配）
  ✅ DEM 分辨率 90m → 30m；海拔 > 3km 警告；海拔跨度 > 500m 时分块校正
  ✅ 反射率反演后 Savitzky-Golay 光谱平滑（FLAASH Spectral Polishing 简化版）
  ✅ 坏波段最终确定推迟到 Step 2 光谱评估后（反射率异常 + 信噪比）
  ✅ 保持 native 网格，不做正射（正射推迟到 Step 3 统一做，带 30m DEM）
  ✅ 输入从 Step 1 的 georef.json 读取 radiance_native.tif（旧版读 radiance_ortho.tif）
  ✅ 多源数据坐标统一：DEM/矢量重投影到 ZY1E 影像的 UTM 带

架构变化：旧 3 步（Step1 含正射 + Step2 + Step3 含正射）→ 新 3 步（Step1 native → Step2 native → Step3 唯一正射）
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

try:                                                 # 检查 numpy 是否可用
    import numpy as np                               # 导入数值计算库
except ImportError as exc:                           # 解释器缺少 numpy
    raise SystemExit(                                # 直接退出并提示
        "当前 Python 解释器缺少 numpy：%s\n"          # 错误说明
        "请用装好 GDAL + Py6S 的解释器运行本脚本，例如：\n"  # 使用提示
        "    E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe 02_atmospheric_correction_v1.1.2.py"  # 示例命令
        % exc)                                       # 附异常信息

try:                                                 # 检查 GDAL 是否可用
    from osgeo import gdal, osr                      # 导入栅格与坐标系模块
except ImportError as exc:                           # 解释器缺少 GDAL
    raise SystemExit(                                # 直接退出并提示
        "当前 Python 解释器缺少 GDAL（osgeo）：%s\n"  # 错误说明
        "请用装好 GDAL + Py6S 的解释器运行本脚本，例如：\n"  # 使用提示
        "    E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe 02_atmospheric_correction_v1.1.2.py"  # 示例命令
        % exc)                                       # 附异常信息

gdal.UseExceptions()                              # GDAL 出错时抛异常而不是返回空
gdal.PushErrorHandler("CPLQuietErrorHandler")     # 屏蔽 GDAL 非致命错误的输出

# v2.2.0：路径/环境统一收口到仓库根 config.py
import os as _os, sys as _sys
_REPO_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
if _REPO_ROOT not in _sys.path:
    _sys.path.insert(0, _REPO_ROOT)
import config as _C

# =================================================================
# 配置区（按需修改）
# =================================================================

# 输入步骤1 正射成果目录；输出大气校正成果根目录
IN_ROOT_STEP1 = _C.OUT_ROOT_STEP1                               # 步骤1成果输入根目录
OUT_ROOT_STEP2 = _C.OUT_ROOT_STEP2                              # 本步大气校正成果输出根目录

# 原始 TAR 目录；读取观测几何和日期，优先使用 georef.json 的 source_tar
TAR_DIR = _C.SRC_DIR                                              # 原始TAR目录，用于读观测几何

# MODIS MOD04/MYD04 AOD 产品目录（按日期及场景足迹选有效像元）

# ASTGTM DEM 目录；可放置 GeoTIFF/IMG 或包含栅格的 ZIP
DEM_DIR = _C.DEM_DIR                                               # DEM目录，用于算海拔

# 6S 可执行程序路径；留空时依次查环境变量、脚本 _bin 目录和 PATH
SIXS_EXE = _C.SIXS_EXE                                  # 6S程序路径，空则自动查找

# v2.1.0 新增常量：此前被引用却未定义，运行期必抛 NameError，现补齐
SENSOR_ALTITUDE_KM = 778.0          # ZY1E 真实轨道高度，仅作元数据；不传给 custom_altitude
# 初始能见度（FLAASH Initial Visibility）：气溶胶无法反演时作为初始值参与 6S 校正。
# ENVI FLAASH 官方默认 40 km（见 https://www.cnblogs.com/enviidl/p/16418433.html：
# "根据影像获取时大气情况设置……此处保持默认即可"，默认值即 FLAASH 的 40 km）。
# 40 km 对应 Koschmieder AOT550≈0.0978（3.912/40），属晴朗水平；可用 --visibility 覆盖。
# 回退过程会明确记入 used_parameters.json 且打 initial_visibility_fallback 标记，非静默。
VISIBILITY_MANUAL_KM = 40.0

# v2.1.0 6S 运行库 PATH：sixs.exe 依赖 libgfortran-3.dll，位于 conda 的 mingw-w64/bin，
# 默认不在 PATH，缺失会导致 exit=127、输出 0 行。启动时自动注入。
SIXS_DLL_DIR = _C.SIXS_DLL_DIR

# v2.1.0 偏移一致性校验阈值：微调成功后比对邻带系数，超阈值说明邻域不平坦，抛错而非静默接受
COEF_CONSISTENCY_MAX_REL_XA = 0.05   # xa 相对偏差上限 5%
COEF_CONSISTENCY_MAX_REL_XBXC = 0.10 # xb/xc 相对偏差上限 10%

# v2.1.0 反射率域波段检测判据（大气校正后判定无信息波段，任一命中即剔）
BAND_DETECT_MIN_VALID_RATIO = 0.05   # 有效像元比例下限，低于 5% 判无信息
BAND_DETECT_MIN_STD = 0.002          # 标准差下限，低于 0.002 判近似常数
BAND_DETECT_MIN_DYNAMIC_RANGE = 0.005  # 动态范围(p99-p1)下限，低于 0.005 判无区分能力

# v2.1.2：水汽强吸收带剔除（对齐 FLAASH bad bands 惯例，见教程博客园 p/16418433）
# 这些波段透过率极低，6S 逐波段校正会放大噪声产生反射率尖峰，污染 SAM/MTMF 与 ML 特征。
# 区间来自高光谱大气校正共识（ScienceDirect: 0.94/1.14/1.38/1.88μm 四个水汽带）。
# 单位 μm，闭区间；None 表示不剔除水汽带。
WATER_VAPOR_BANDS_UM = [
    (0.930, 0.965),   # 0.94 μm 水汽带
    (1.120, 1.160),   # 1.14 μm 水汽带
    (1.350, 1.420),   # 1.38 μm 水汽带（最强）
    (1.820, 1.960),   # 1.88 μm 水汽带（含 Step1 登记的 DN=0 段 1.812-1.964）
]
DROP_WATER_VAPOR_BANDS = True        # 是否剔除上述水汽带（用户 20:59 确认采用）
OUT_BAND_DETECT_REPORT = "band_detection_report.csv"  # 波段检测报告
OUT_COEF_REPAIR_LOG = "coef_repair_log.json"          # NaN 修复台账

# 自动按日期和纬度选择大气模式；关闭时使用手工模式
AUTO_DETECT_ATMOS = True                     # 按日期和纬度自动选大气模式
ATMOS_PROFILE_MANUAL = "MidlatitudeWinter"  # Tropical/MidlatitudeSummer/MidlatitudeWinter/SubarcticSummer/SubarcticWinter/USStandard

# 气溶胶类型（无土地覆盖数据时默认 Continental）
AERO_PROFILE_MANUAL = "Continental"  # Continental/Urban/Maritime/Desert
AUTO_DETECT_AERO = False             # 启用时按经纬度粗略判断

# 暗像元反射比 AOT 模式：None / KT_DARK_PIXEL / 2-BAND_OVER_WATER
AOT_METHOD = "KT_DARK_PIXEL"                 # 直接指定暗像元反射比算法
AOT_METHOD_CHOICES = ("None", "2-Band (K-T)", "2-Band Over Water")

# v2.1.2 FLAASH 标准 K-T 参数（修正 v2.1.0 正向LUT地表反射率=0 致 AOT 顶格 0.6 的 bug）
KT_DARK_PERCENTILE = 5.0          # 暗像元候选：SWIR 最低百分位（FLAASH 用低反射率 cutoff）
KT_RED_NIR_RATIO_MAX = 0.6        # 阴影/水体排除：暗植被 red/NIR 低（<0.6），阴影 red≈NIR（≈1）被剔
KT_RED_SWIR_EXPECTED = 0.45       # FLAASH 暗植被 red/SWIR 地表反射率经验比值（干湿气候均值，论文 L195）
KT_RED_SWIR_TOL = 0.12            # 比值容差：|red/swir − 0.45| 最小者定 AOT
KT_SWIR_REF_MIN = 0.005           # 反演出的 SWIR 地表反射率下限，低于此视为无效暗像元
KT_NIR_TARGET_UM = 0.860          # 阴影排除用的 NIR 波段中心波长（μm）
AOT_SOURCE = "manual"                        # compatibility：manual/auto 已废弃为算法来源标识
AOT550_MANUAL = 0.2                  # 自动来源无有效值时使用的 AOT550
AOT_MIN_VALID = 0.0                  # MODIS 有效 AOD 开区间下限
AOT_MAX_VALID = 5.0                  # MODIS 有效 AOD 闭区间上限
AERONET_ENABLED = True               # 是否允许联网查询 AERONET
AERONET_MAX_DISTANCE_KM = 300        # AERONET 站点最远距离（km）
AERONET_SEARCH_DAYS = 3              # AERONET 日期搜索范围（天）

# 自动从足迹内 DEM 计算平均海拔；失败时使用固定兜底值
AUTO_DETECT_DEM = True               # 启用以从足迹内DEM自动算平均海拔
DEM_RES = 30.0                       # v2.0.0: 30m DEM（旧版 90m）
ELEVATION_KM_FALLBACK = 0.2          # DEM 不可用时的海拔（km）
# v2.2.0：高程分块 6S——把影像按海拔分档，每档分别算 6S 系数（山区精度提升）。
# 线性覆盖 [min, max] 高程范围；ELEV_BAND_COUNT=1 等价于单海拔（关闭分块）。
ELEV_BAND_COUNT = 3
ELEV_BAND_COARSE_STEP = 120          # 原生网格降采样步长（约 120x120 粗网格反算海拔）

# 水汽为空时使用 6S 默认廓线；指定水汽后可同时指定臭氧
WATER_VAPOR_CM = None               # 水汽含量（g/cm²），None 表示使用大气模式廓线
OZONE_CM = None                     # 臭氧含量（cm-atm），None 时使用下方默认值
OZONE_CM_DEFAULT = 0.3              # 手工水汽廓线下的默认臭氧量（cm-atm）

# 反射率裁剪范围；None 表示不裁剪
REFLECTANCE_CLIP = (0.0, 1.0)        # 反射率裁剪上下限
MASK_INVALID_RADIANCE = True         # 是否掩膜低于阈值的输入辐亮度
RADIANCE_MIN_VALID = 0.0             # 输入辐亮度有效下限
NODATA_OUT = -9999.0                 # 输出无效像元值

# XML 缺少观测几何时使用的手工兜底值（角度）
SOLAR_ZENITH_MANUAL = None           # 手工太阳天顶角（度），None表示不用
SOLAR_AZIMUTH_MANUAL = None          # 手工太阳方位角（度）
VIEW_ZENITH_MANUAL = None            # 手工观测天顶角（度）
VIEW_AZIMUTH_MANUAL = None           # 手工观测方位角（度）

# 批处理和 GDAL 参数
OVERWRITE = False                    # 是否重做已有成果
GDAL_CACHE_MB = 1024                 # GDAL 缓存大小（MB）

# GeoTIFF 写出设置
GTIFF_OPTIONS = ["BIGTIFF=IF_SAFER"]  # 大文件时自动启用BigTIFF
GEOTIFF_COMPRESS = "DEFLATE"         # DEFLATE/LZW/None
GEOTIFF_COMPRESS_LEVEL = 6           # 压缩等级

# 质检图设置
MAKE_QC = True                       # 是否生成 preview_quality.png（需要 matplotlib）
QC_GEOREF_PREVIEW = True             # 是否生成带地理参考的降采样预览
QC_PREVIEW_SCALE = 2.0               # 预览像元大小相对输入的倍数
QC_MAX_SIZE = 1400                   # 普通预览图长边最大像素数
QC_DPI = 150                         # 质检图分辨率

# 步骤1 输入文件名：新命名优先，同时兼容旧版固定名称
IN_RADIANCE_SUFFIX = "_native_radiance.tif"      # v2.0.0: Step 1 输出 native 网格辐亮度
LEGACY_RADIANCE_NAME = "radiance_native.tif"    # v2.0.0: Step 1 原生网格辐亮度
IN_WAVELENGTHS_NAME = "wavelengths_clean.npy"   # 步骤1输出的波长数组文件
IN_GEOREF_NAME = "georef.json"                  # 步骤1输出的几何元数据文件

# 本步输出文件名：保留场景名，并以缩写标注正射与大气校正处理
OUT_REFLECTANCE_SUFFIX = "_cal_badband_ortho_atmcorr.tif"  # 输出反射率文件名后缀
LEGACY_REFLECTANCE_NAMES = (                               # 旧版反射率文件名列表
    "_ortho_atmcorr.tif",                                  # 旧命名一
    "reflectance_atm.tif",                                 # 旧命名二
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
AERO_PROFILE_CHOICES = ("Continental", "Urban", "Maritime", "Desert")  # 可选气溶胶类型
ATMOS_PROFILE_CHOICES = ("Tropical", "MidlatitudeSummer", "MidlatitudeWinter",  # 可选大气模式
                         "SubarcticSummer", "SubarcticWinter", "USStandard")     # 其余大气模式


# =================================================================
# 通用工具
# =================================================================

def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(                              # 创建命令行参数解析器
        description="ZY1E AHSI 步骤2：6S 大气校正（正射格网，直接承接步骤1 正射辐亮度）")  # 程序说明
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP1,  # 输入根目录参数
                        help="步骤1 成果根目录（默认: %(default)s）")   # 参数帮助文本
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP2,  # 输出根目录参数
                        help="本步成果根目录，每景生成同名子文件夹（默认: %(default)s）")  # 参数帮助文本
    parser.add_argument("--scene", action="append", default=None,  # 场景名过滤，可多次指定
                        help="只处理匹配的场景名，可多次指定，支持通配符（如 \"*E117.20*\"）")  # 参数帮助文本
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE,  # 覆盖已有成果
                        help="已有成果时覆盖重做（默认跳过）")       # 参数帮助文本
    parser.add_argument("--limit", type=int, default=0,            # 限制处理的景数
                        help="只处理前 N 景，便于调试（默认 0 = 全部）")  # 参数帮助文本
    parser.add_argument("--dry-run", action="store_true",          # 只演练不写盘
                        help="只列出待处理场景与目标目录，不执行")   # 参数帮助文本
    parser.add_argument("--no-qc", action="store_true",            # 关闭质检图生成
                        help="跳过质检图生成（保留统计 CSV 与预览 tif 开关不受影响）")  # 参数帮助文本
    parser.add_argument("--log", default=None,                     # 自定义日志文件路径
                        help="日志文件路径（默认 <out>/step2_batch.log）")  # 参数帮助文本
    parser.add_argument("--verbose", action="store_true", help="输出调试日志")  # 打开DEBUG日志
    # v2.1.3：与 Step1/Step3 统一——默认出错即停，显式 --continue-on-error 才继续
    parser.add_argument("--continue-on-error", action="store_true", default=False,
                        help="单景失败时继续下一景（默认关闭=出错即停；"
                             "无视错误会导致数据不可用，仅明确知情时开启）")
    # v2.2.0：原始 tar / DEM 目录可由命令行覆盖，便于换机运行
    parser.add_argument("--tar-dir", default=TAR_DIR,
                        help="原始 L1A 压缩包目录（读观测几何，默认: %(default)s）")
    parser.add_argument("--dem-dir", default=DEM_DIR,
                        help="ASTGTM DEM 目录（默认: %(default)s）")

    # 手工覆盖 6S 关键参数（便于快速试验，不写死文件）
    parser.add_argument("--aot", type=float, default=None,         # 手工指定AOT550
                        help="直接指定 AOT550（覆盖自动 AOD 检索）")  # 参数帮助文本
    parser.add_argument("--aot-method", choices=AOT_METHOD_CHOICES,  # 暗像元反射比算法
                        default=None,
                        help="AOT 算法：None / 2-Band (K-T) / 2-Band Over Water")
    parser.add_argument("--visibility", type=float, default=None,  # v2.1.0 初始能见度
                        help="初始能见度 km（FLAASH Initial Visibility）：气溶胶无法反演时"
                             "作为初始值参与校正，默认 %g km" % VISIBILITY_MANUAL_KM)
    parser.add_argument("--elev", type=float, default=None,        # 手工指定海拔
                        help="直接指定目标海拔（km，覆盖 DEM 自动计算）")  # 参数帮助文本
    parser.add_argument("--atmos", default=None, choices=ATMOS_PROFILE_CHOICES,  # 指定大气模式
                        help="直接指定 6S 大气模式")                # 参数帮助文本
    parser.add_argument("--aero", default=None, choices=AERO_PROFILE_CHOICES,  # 指定气溶胶类型
                        help="直接指定 6S 气溶胶类型")              # 参数帮助文本
    return parser.parse_args(argv)                                 # 解析并返回参数对象


def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)                 # 取得具名日志器，避免重名
    logger.setLevel(logging.DEBUG)                       # 记录器总级别设为DEBUG
    logger.handlers.clear()                              # 清空旧处理器，防止重复输出
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s",  # 统一日志格式
                            datefmt="%Y-%m-%d %H:%M:%S")  # 时间格式
    console = logging.StreamHandler(stream=sys.stdout)    # 控制台处理器输出到标准输出
    console.setLevel(logging.DEBUG if verbose else logging.INFO)  # 按--verbose定级别
    console.setFormatter(fmt)                            # 套用统一格式
    logger.addHandler(console)                           # 挂到日志器
    if log_path:                                         # 指定了日志文件才写文件
        try:                                             # 写文件失败不影响主流程
            parent = os.path.dirname(os.path.abspath(log_path))  # 取日志所在目录
            if parent:                                   # 目录非空才创建
                os.makedirs(parent, exist_ok=True)       # 确保日志目录存在
            file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")  # 追加写文件
            file_handler.setLevel(logging.DEBUG)         # 文件始终记录最详细级别
            file_handler.setFormatter(fmt)               # 套用统一格式
            logger.addHandler(file_handler)              # 挂到日志器
        except OSError as exc:                           # 日志文件无法创建
            console.handle(logging.LogRecord(            # 手动构造记录输出告警
                LOG_NAME, logging.WARNING, __file__, 0,  # 记录器名、级别、来源行
                "无法写入日志文件 %s: %s", (log_path, exc), None))  # 告警内容与参数
    logger.propagate = False                             # 不再向根日志器传播
    return logger                                        # 返回配置好的日志器


def human_size(num_bytes):
    """把字节数格式化为 B/KB/MB/GB/TB 文本。"""
    size = float(num_bytes)                              # 转成浮点便于逐级换算
    for unit in ("B", "KB", "MB", "GB", "TB"):           # 依次按单位降级
        if size < 1024.0 or unit == "TB":                # 小于1024或已到最大单位
            return "%.1f %s" % (size, unit)              # 返回带单位的文本
        size /= 1024.0                                   # 继续除以1024换下一单位
    return "%.1f TB" % size                              # 兜底返回TB格式


def discover_scene_dirs(in_root):
    """返回步骤1 成果子目录列表；若传入的就是单景目录则直接返回它。"""
    if not os.path.isdir(in_root):                       # 输入目录必须存在
        raise NotADirectoryError("步骤1 成果目录不存在: %s" % in_root)  # 不存在直接报错
    if resolve_radiance_path(in_root):                   # 本身就是单景目录
        return [in_root]                                 # 直接作为唯一候选返回
    scene_dirs = sorted(glob.glob(os.path.join(in_root, "*")))  # 列出根目录下全部条目
    return [path for path in scene_dirs                  # 过滤出真正的景目录
            if os.path.isdir(path) and resolve_radiance_path(path)]  # 需含辐亮度文件


def resolve_radiance_path(scene_dir):
    """优先返回新命名的步骤1辐亮度，其次兼容旧名 radiance_ortho.tif。"""
    renamed = sorted(glob.glob(os.path.join(             # 优先查找新命名辐亮度
        scene_dir, "*" + IN_RADIANCE_SUFFIX)))           # 按后缀通配匹配
    if renamed:                                          # 找到新命名文件
        return renamed[0]                                # 取第一个返回
    legacy = os.path.join(scene_dir, LEGACY_RADIANCE_NAME)  # 回退到旧固定文件名
    return legacy if os.path.isfile(legacy) else None    # 存在则返回，否则返回None


def fmt_elapsed(seconds):
    """把秒数格式化为“分:秒.十分位”，便于日志阅读。"""
    seconds = max(0.0, float(seconds))                   # 负值归零，统一为浮点
    return "%d:%04.1f" % (int(seconds // 60), seconds % 60.0)  # 拼成“分:秒.十分位”


def log_stage_start(logger, number, total, description):
    """打印阶段开始日志并返回计时起点。"""
    logger.info("[步骤 %d/%d] 开始：%s", number, total, description)  # 打印阶段开始信息
    return time.time()                                   # 返回计时起点


def log_stage_done(logger, number, total, description, started):
    """打印阶段完成日志与耗时。"""
    logger.info("[步骤 %d/%d] 完成：%s，用时 %s",          # 打印阶段完成信息
                number, total, description, fmt_elapsed(time.time() - started))  # 附带耗时


def atomic_write_json(payload, path):
    """先写 .tmp 再 os.replace，避免中断留下半个 JSON。"""
    tmp_path = path + ".tmp"                             # 先写同目录临时文件
    with open(tmp_path, "w", encoding="utf-8") as handle:  # 以UTF-8打开临时文件
        json.dump(payload, handle, ensure_ascii=False, indent=2)  # 写入缩进JSON
    os.replace(tmp_path, path)                           # 原子替换为目标文件


def atomic_write_text(text, path):
    """先写 .tmp 再 os.replace，原子写出文本文件。"""
    tmp_path = path + ".tmp"                             # 先写同目录临时文件
    with open(tmp_path, "w", encoding="utf-8") as handle:  # 以UTF-8打开临时文件
        handle.write(text)                               # 写入文本内容
    os.replace(tmp_path, path)                           # 原子替换为目标文件


def copy_file_atomic(src, dst):
    """复制到 .tmp 后原子替换目标文件。"""
    tmp_path = dst + ".tmp"                              # 先复制到临时文件
    shutil.copyfile(src, tmp_path)                       # 复制源文件到临时路径
    os.replace(tmp_path, dst)                            # 原子替换目标文件


def load_georef_json(path):
    """读取步骤1 的 georef.json。"""
    with open(path, "r", encoding="utf-8") as handle:    # 以UTF-8读取元数据文件
        return json.load(handle)                         # 解析JSON并返回字典


# =================================================================
# 产品元数据解析
# =================================================================

def locate_tar(scene, georef, tar_dir, logger):
    """定位原始 L1A tar：优先 georef.json 记录的 source_tar，再按景名匹配。"""
    """定位原始 tar 压缩包：先按 georef.json 记录，再按场景名匹配。"""
    candidates = []                                      # 候选 tar 路径列表
    if georef and georef.get("source_tar"):              # georef 记录了来源 tar
        candidates.append(os.path.join(tar_dir, georef["source_tar"]))  # 拼到TAR目录下
        candidates.append(georef["source_tar"])          # 也考虑记录中的绝对路径
    for ext in (".tar", ".tar.gz", ".tgz"):              # 尝试三种压缩包后缀
        candidates.append(os.path.join(tar_dir, scene + ext))  # 按景名拼路径
    for path in candidates:                              # 逐个验证候选路径
        if path and os.path.isfile(path):                # 路径非空且文件存在
            return os.path.abspath(path)                 # 返回绝对路径

    # 最后兜底：按场景名在 TAR_DIR 里模糊匹配
    hits = []                                            # 模糊匹配结果列表
    for ext in (".tar", ".tar.gz", ".tgz"):              # 同样三种后缀
        hits.extend(glob.glob(os.path.join(tar_dir, "*" + scene + "*" + ext)))  # 通配匹配
    hits = sorted(set(os.path.abspath(path) for path in hits))  # 去重并排序
    if hits:                                             # 至少匹配到一个
        if len(hits) > 1:                                # 匹配到多个时提示
            logger.warning("场景 %s 匹配到 %d 个压缩包，使用第一个: %s",  # 多包告警
                           scene, len(hits), os.path.basename(hits[0]))  # 告警参数
        return hits[0]                                   # 取第一个匹配包
    return None                                          # 完全找不到返回None


def decode_text(raw):
    """按 UTF-8 宽松解码并统一换行、去掉 BOM。"""
    text = raw.decode("utf-8", errors="ignore")          # 宽松解码，坏字节忽略
    return text.lstrip("\ufeff").replace("\r\n", "\n").replace("\r", "\n")  # 去BOM并统一换行


def read_product_xml(tar_path, logger):
    """从 L1A 包内读取产品元数据 XML 文本。"""
    """从 tar 中读取产品元数据 XML 文本。"""
    with tarfile.open(tar_path, "r:*") as tar:           # 自动识别压缩方式打开tar
        hits = [name for name in tar.getnames()          # 收集包内XML成员名
                if name.lower().endswith(".xml")         # 只看xml后缀
                and not name.endswith("/")]              # 排除目录条目
        if not hits:                                     # 一个XML都没有
            raise FileNotFoundError("压缩包中没有 XML 元数据: %s" % tar_path)  # 直接报错
        # 优先包含 ProductMetaData 的产品 XML
        preferred = None                                 # 选中的(名,文本)对
        for name in sorted(hits):                        # 按名称排序逐个尝试
            handle = tar.extractfile(name)               # 取出成员文件对象
            if handle is None:                           # 目录等无效成员
                continue                                 # 跳过
            with handle:                                 # 用完自动关闭
                text = decode_text(handle.read())        # 读内容并解码
            if "ProductMetaData" in text or "SolarZenith" in text:  # 命中产品元数据特征
                preferred = (name, text)                 # 记为选中的XML
                break                                    # 不再继续找
        if preferred is None:                            # 没有特征命中的XML
            name = sorted(hits)[0]                       # 退回第一个XML
            handle = tar.extractfile(name)               # 取出该成员
            if handle is None:                           # 无法读取
                raise FileNotFoundError("无法读取 XML: %s" % name)  # 明确报错
            with handle:                                 # 用完自动关闭
                preferred = (name, decode_text(handle.read()))  # 读取并解码
        if len(hits) > 1:                                # 包内有多份XML
            logger.debug("压缩包内 XML %d 个，使用: %s", len(hits), preferred[0])  # 记录选中项
        return preferred[1]                              # 返回XML文本


def _find_float(text, tag):
    match = re.search(r"<%s>\s*([-+0-9.eE]+)\s*</%s>" % (tag, tag),  # 正则匹配数值标签
                      text, re.IGNORECASE)               # 忽略大小写
    if not match:                                        # 没有匹配
        return None                                      # 返回None
    try:                                                 # 数值转换可能失败
        return float(match.group(1))                     # 返回浮点值
    except ValueError:                                   # 内容不是合法数字
        return None                                      # 视为缺失


def parse_product_metadata(xml_text, scene, logger):
    """解析成像日期、四角经纬度与卫星/传感器等产品元数据。"""
    """解析产品 XML：观测几何、日期、中心位置、四角足迹。"""
    meta = {                                             # 汇总元数据字典
        "satellite": (_find_text(xml_text, "SatelliteID") or "ZY1E"),  # 卫星，缺省ZY1E
        "sensor": (_find_text(xml_text, "SensorID") or "AHSI"),        # 传感器，缺省AHSI
        "solar_z": _find_float(xml_text, "SolarZenith"),   # 太阳天顶角
        "solar_a": _find_float(xml_text, "SolarAzimuth"),  # 太阳方位角
        "view_z": _find_float(xml_text, "SatelliteZenith"),  # 观测天顶角
        "view_a": _find_float(xml_text, "SatelliteAzimuth"),  # 观测方位角
        "center_lat": _find_float(xml_text, "CenterLatitude"),   # 中心纬度
        "center_lon": _find_float(xml_text, "CenterLongitude"),  # 中心经度
        "date": None,                                    # 成像日期，稍后解析
        "footprint": {},                                 # 四角足迹
    }

    match = re.search(r"<StartTime>\s*(\d{4})-(\d{2})-(\d{2})"  # 匹配成像开始时间
                      r"(?:[ T](\d{2}):(\d{2}):(\d{2}))?", xml_text)  # 时分秒可选
    if match:                                            # 找到开始时间
        parts = [int(value) for value in match.groups()[:6] if value is not None]  # 转整数
        if len(parts) == 3:                              # 只有日期没有时间
            parts += [0, 0, 0]                           # 补零为时分秒
        try:                                             # 日期可能非法
            meta["date"] = datetime.datetime(*parts)     # 构造成像时间
        except ValueError:                               # 越界等非法时间
            logger.warning("StartTime 解析异常: %s", match.group(0))  # 记录告警

    corner_map = {"TL": "TopLeft", "TR": "TopRight",     # 角点缩写到标签前缀
                  "BR": "BottomRight", "BL": "BottomLeft"}  # 其余两角
    for abbr, corner in corner_map.items():              # 逐个角点解析经纬度
        lat = _find_float(xml_text, corner + "Latitude")  # 取该角纬度
        lon = _find_float(xml_text, corner + "Longitude")  # 取该角经度
        if lat is not None and lon is not None:          # 经纬度都有效
            meta["footprint"][abbr] = (lon, lat)         # 记录(经,纬)坐标

    missing = [key for key in ("solar_z", "solar_a", "view_z", "view_a")  # 检查几何字段
               if meta[key] is None]                     # 缺失字段清单
    if missing:                                          # 有缺失则告警
        logger.warning("产品 XML 缺少观测几何: %s", ", ".join(missing))  # 拼接缺失项
    return meta                                          # 返回解析结果


def _find_text(text, tag):
    match = re.search(r"<%s>\s*([^<]+?)\s*</%s>" % (tag, tag), text, re.IGNORECASE)  # 匹配文本标签
    return match.group(1).strip() if match else None     # 去空白后返回，未命中返回None


def parse_scene_name(scene):
    """从景名中解析成像日期与经纬度（解析失败返回 None 字段）。"""
    """从场景名兜底解析日期/中心经纬度。"""
    info = {"date": None, "center_lat": None, "center_lon": None}  # 待填的日期与经纬度
    match = re.search(r"_(\d{8})_", scene)               # 匹配景名中的8位日期
    if match:                                            # 找到日期片段
        try:                                             # 非法日期可能抛错
            info["date"] = datetime.datetime.strptime(match.group(1), "%Y%m%d")  # 解析日期
        except ValueError:                               # 日期格式非法
            pass                                         # 保持None
    match = re.search(r"_([EW])(\d{1,3}\.\d+)_([NS])(\d{1,3}\.\d+)", scene)  # 匹配经纬度
    if match:                                            # 找到经纬度片段
        lon = float(match.group(2)) * (1.0 if match.group(1) == "E" else -1.0)  # 东经为正
        lat = float(match.group(4)) * (1.0 if match.group(3) == "N" else -1.0)  # 北纬为正
        info["center_lon"] = lon                         # 记录中心经度
        info["center_lat"] = lat                         # 记录中心纬度
    return info                                          # 返回解析结果


def resolve_metadata(scene, scene_dir, georef, tar_dir, logger):
    """汇总观测几何：XML → georef.json → 景名/手工兜底，返回统一字典。"""
    """汇总本景 6S 计算所需的元数据：观测几何、日期、中心位置、足迹、tar 路径。"""
    name_info = parse_scene_name(scene)                  # 先从景名解析兜底信息
    meta = {                                             # 本景元数据字典
        "scene": scene,                                  # 场景名
        "tar": None,                                     # 原始压缩包路径
        "solar_z": None, "solar_a": None, "view_z": None, "view_a": None,  # 观测几何
        "date": name_info["date"],                       # 成像日期（来自景名）
        "center_lat": name_info["center_lat"],           # 中心纬度（来自景名）
        "center_lon": name_info["center_lon"],           # 中心经度（来自景名）
        "footprint": {},                                 # 四角足迹
        "satellite": "ZY1E", "sensor": "AHSI",           # 默认卫星与传感器
    }

    tar_path = locate_tar(scene, georef, tar_dir, logger)  # 定位原始压缩包
    if tar_path:                                         # 找到压缩包
        meta["tar"] = tar_path                           # 记录压缩包路径
        try:                                             # XML 读取失败可降级
            xml_text = read_product_xml(tar_path, logger)  # 读产品XML文本
            xml_meta = parse_product_metadata(xml_text, scene, logger)  # 解析XML元数据
            for key in ("solar_z", "solar_a", "view_z", "view_a",  # 覆盖观测几何
                        "satellite", "sensor"):              # 以及卫星/传感器
                if xml_meta.get(key) is not None:        # XML 有值才覆盖
                    meta[key] = xml_meta[key]            # 写入元数据
            if xml_meta.get("date") is not None:         # XML 提供了日期
                meta["date"] = xml_meta["date"]          # 采用XML日期
            for key in ("center_lat", "center_lon"):     # 覆盖中心经纬度
                if xml_meta.get(key) is not None:        # XML 有值才覆盖
                    meta[key] = xml_meta[key]            # 写入元数据
            if xml_meta.get("footprint"):                # XML 提供了足迹
                meta["footprint"] = dict(xml_meta["footprint"])  # 复制足迹字典
        except Exception as exc:                         # 任何读取异常都降级
            logger.warning("读取产品 XML 失败，改用 georef.json/场景名兜底: %s", exc)  # 告警
    else:                                                # 没找到压缩包
        logger.warning("未找到原始压缩包，观测几何将使用手工兜底值")  # 提示用手工值

    # 足迹优先用步骤1 的 georef.json 的 footprint_wgs84（WGS84，与步骤1/步骤3 一致）
    if georef and georef.get("footprint_wgs84"):         # 步骤1提供了WGS84足迹
        footprint = {}                                   # 新的足迹字典
        for key, value in georef["footprint_wgs84"].items():  # 遍历四角
            footprint[key] = (float(value[0]), float(value[1]))  # 转成浮点(经,纬)
        meta["footprint"] = footprint                    # 覆盖为步骤1足迹

    # 手工几何兜底
    manual = {"solar_z": SOLAR_ZENITH_MANUAL, "solar_a": SOLAR_AZIMUTH_MANUAL,  # 太阳角兜底
              "view_z": VIEW_ZENITH_MANUAL, "view_a": VIEW_AZIMUTH_MANUAL}       # 观测角兜底
    for key, value in manual.items():                    # 逐个检查手工配置
        if meta[key] is None and value is not None:      # 缺失且有手工值
            meta[key] = float(value)                     # 填入手工值
            logger.warning("观测几何 %s 使用手工兜底值 %s", key, value)  # 记录告警

    missing = [key for key in ("solar_z", "solar_a", "view_z", "view_a")  # 检查几何字段
               if meta[key] is None]                     # 仍缺失的字段
    if missing:                                          # 有缺失则终止
        raise ValueError(                                # 缺少几何无法做6S
            "缺少观测几何参数 %s；请在脚本配置区填写 *_MANUAL 兜底值，"  # 提示补手工值
            "或确认 TAR_DIR 指向原始压缩包" % ", ".join(missing))       # 拼缺失项

    if not meta["footprint"]:                            # 足迹为空
        raise ValueError(                                # 无足迹无法筛选像元
            "缺少四角足迹（georef.json 与产品 XML 均未提供）；"        # 提示原因
            "请确认步骤1 正常完成（georef.json）")                    # 提示检查步骤1

    if meta["center_lat"] is None or meta["center_lon"] is None:  # 中心经纬度缺失
        lons = [value[0] for value in meta["footprint"].values()]  # 收集经度
        lats = [value[1] for value in meta["footprint"].values()]  # 收集纬度
        meta["center_lon"] = float(np.mean(lons))        # 用足迹经度均值兜底
        meta["center_lat"] = float(np.mean(lats))        # 用足迹纬度均值兜底
        logger.warning("中心经纬度缺失，用足迹均值兜底: (%.4f, %.4f)",  # 记录兜底值
                       meta["center_lat"], meta["center_lon"])  # 打印经纬度

    if meta["date"] is None:                             # 成像日期缺失
        meta["date"] = datetime.datetime(2000, 1, 1)     # 用固定日期兜底
        logger.warning("成像日期缺失，6S 月份/日期使用 2000-01-01 兜底")  # 记录告警

    return meta                                          # 返回汇总后的元数据


def footprint_ring(footprint):
    """把四角点足迹转成闭合环 [[lon,lat],...]，供点在多边形判断。"""
    """返回闭合的 WGS84 足迹环 [[lon, lat], ...]。"""
    order = [key for key in ("TL", "TR", "BR", "BL") if key in footprint]  # 按顺时针取角
    if len(order) < 3:                                   # 角点不足三个
        order = list(footprint.keys())                   # 退回全部已有键
    ring = [[float(footprint[key][0]), float(footprint[key][1])] for key in order]  # 转浮点坐标
    if ring:                                             # 环非空
        ring.append(list(ring[0]))                       # 首点追加到末尾闭合
    return ring                                          # 返回闭合环


def points_in_polygon(x, y, ring):
    """射线法判断像元中心 (x, y) 是否落在足迹环内，返回布尔数组。"""
    """射线法判断点是否在闭合多边形内（ring 首尾可重复）。"""
    x = np.asarray(x, dtype=np.float64)                  # 统一转为浮点数组
    y = np.asarray(y, dtype=np.float64)                  # 统一转为浮点数组
    inside = np.zeros(x.shape, dtype=bool)               # 初始都在外部
    points = np.asarray(ring, dtype=np.float64)          # 足迹环转数组
    if points.shape[0] > 1 and np.allclose(points[0], points[-1]):  # 首尾重复
        points = points[:-1]                             # 去掉重复的末点
    count = points.shape[0]                              # 环上的点数
    for idx in range(count):                             # 逐条边做射线判断
        x1, y1 = points[idx]                             # 当前边起点
        x2, y2 = points[(idx + 1) % count]               # 当前边终点（环状取模）
        cond = (y1 > y) != (y2 > y)                      # 射线是否穿过该边
        if not np.any(cond):                             # 该边不影响任何点
            continue                                     # 跳过此边
        with np.errstate(divide="ignore", invalid="ignore"):  # 忽略除零告警
            xin = (x2 - x1) * (y - y1) / (y2 - y1) + x1  # 求交点横坐标
        inside ^= cond & (x < xin)                       # 交点右侧则翻转内外
    return inside                                        # 返回布尔掩膜


# =================================================================
# AOT550：MODIS MOD04_L2 / AERONET / 手工
# =================================================================


# ============================================================
# v2.0.0: K-T 暗像元 2-Band 气溶胶反演（替代旧 MODIS/AERONET）
# Kaufman & Tanré (1992): 用 SWIR 和红光波段在暗像元上的辐亮度差反演 AOT550
# FLAASH 标准波段：SWIR 2100-2250 nm (中心 2175 μm) + Red 640-680 nm (中心 660 μm)
# ============================================================

def _nearest_band(wavelengths_um, target_um):
    """找最接近目标波长的波段索引。"""
    return int(np.argmin(np.abs(np.asarray(wavelengths_um) - target_um)))


def visibility_to_aot550(visibility_km, wavelength=0.550):
    """Koschmieder 公式：能见度 → AOT550。
    AOT = 3.912 / V_km × (λ/550)^(-1.3)
    晴朗 40-100km → AOT 0.04-0.1；中等雾霾 20-30km → AOT 0.13-0.2
    """
    return 3.912 / visibility_km * (wavelength / 0.550) ** (-1.3)


def resolve_aot550_kt(wavelengths_um, radiance_cube, meta, sixs_builder_func, logger):
    """
    K-T 暗像元反演 AOT550（整景均匀值）—— v2.1.2 按 FLAASH 标准重构。

    ★ 修正 v2.1.0 的致命 bug：旧版用"地表反射率=0 的 6S 正向辐亮度 LUT"匹配暗像元红光
      辐亮度。地表 rho=0 时 6S 模拟辐亮度≈0（实测 rho=0→0.0, rho=0.02→5.87），
      任何真实辐亮度都远超 LUT 上限 → AOT 恒顶到 0.6。这是 AOT 系统性高估的根因。

    新流程（FLAASH 2-Band K-T，论文 Adler-Golden 2008 §2.2 正统实现）：
      1. 定位 SWIR(~2.175μm)、Red(~0.66μm)、NIR(~0.86μm) 波段
      2. 暗像元 = SWIR 辐亮度最低 KT_DARK_PERCENTILE% 且 red/NIR < KT_RED_NIR_RATIO_MAX
         （red/NIR 判据排除阴影与水体：阴影 red≈NIR≈1 被剔，暗植被 NIR>>red 保留）
      3. 对候选 AOT ∈ [0, 0.6]，用 6S 系数(xa,xb,xc)把暗像元 SWIR/Red 辐亮度反演为
         地表反射率 ρ_swir、ρ_red（ρ=(xa·L−xb)/(1+xc·(xa·L−xb))）
      4. 取使 ρ_red/ρ_swir 最接近 FLAASH 经验比值 0.45 的 AOT（比值法消除绝对定标影响）
      5. 反演失败/暗像元不足 → 回退初始能见度（由 resolve_aot550 统一处理）

    参数：
      wavelengths_um: 波段中心波长数组（μm），须含 SWIR/Red/NIR
      radiance_cube: (n_bands, H, W) 辐亮度（含 SWIR/Red/NIR 三波段）
      sixs_builder_func: 构造 Py6S 实例的函数（传 aot/wavelength）

    返回 (aot550, info_dict)
    """
    # --- Step 1: 定位 K-T 波段（SWIR + Red + NIR）---
    swir_idx = _nearest_band(wavelengths_um, 2.175)       # SWIR ~2.175 μm
    red_idx = _nearest_band(wavelengths_um, 0.660)        # Red ~0.66 μm
    nir_idx = _nearest_band(wavelengths_um, KT_NIR_TARGET_UM)  # NIR ~0.86 μm（阴影判据）
    swir_wl = float(wavelengths_um[swir_idx])
    red_wl = float(wavelengths_um[red_idx])
    nir_wl = float(wavelengths_um[nir_idx])

    kt_ok = (2.10 <= swir_wl <= 2.25) and (0.64 <= red_wl <= 0.68)
    logger.info("K-T 波段(FLAASH 2-Band): SWIR[%d]=%.4fμm, Red[%d]=%.4fμm, "
                "NIR[%d]=%.4fμm(阴影判据)", swir_idx + 1, swir_wl,
                red_idx + 1, red_wl, nir_idx + 1, nir_wl)
    if not kt_ok:
        logger.warning("K-T 主波段偏离 FLAASH 标准范围，反演精度可能下降")

    n_bands = int(radiance_cube.shape[0])
    if max(swir_idx, red_idx, nir_idx) >= n_bands:
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "reason": "kt_bands_missing",
        }

    # --- Step 2: 暗像元选择（SWIR 低百分位 + red/NIR 排除阴影水体）---
    swir_l = radiance_cube[swir_idx].astype(np.float64)
    red_l = radiance_cube[red_idx].astype(np.float64)
    nir_l = radiance_cube[nir_idx].astype(np.float64)
    # 有效：有限、正值、排除 nodata(-9999) 与离群超大值
    valid = (np.isfinite(swir_l) & np.isfinite(red_l) & np.isfinite(nir_l)
             & (swir_l > 0.1) & (red_l > 0) & (nir_l > 0)
             & (swir_l < 100) & (red_l < 200) & (nir_l < 200))
    n_valid = int(valid.sum())
    if n_valid < 1000:
        logger.warning("有效像元不足（%d），K-T 反演失效", n_valid)
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "visibility_km": VISIBILITY_MANUAL_KM,
            "reason": "insufficient_valid_pixels",
        }

    swir_thr = float(np.percentile(swir_l[valid], KT_DARK_PERCENTILE))
    candidate = valid & (swir_l <= swir_thr)
    red_nir = red_l / np.maximum(nir_l, 1e-9)
    dark = candidate & (red_nir < KT_RED_NIR_RATIO_MAX)  # 排除阴影(red≈NIR)
    n_shadow_rejected = int((candidate & ~dark).sum())
    if int(dark.sum()) < 10:                              # 判据过严时放宽（保留候选）
        logger.warning("阴影判据后暗像元不足（%d），放宽为纯 SWIR 低值候选", int(dark.sum()))
        dark = candidate
    n_dark = int(dark.sum())
    if n_dark < 10:
        logger.warning("暗像元不足（%d 个），K-T 反演失效", n_dark)
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "visibility_km": VISIBILITY_MANUAL_KM,
            "reason": "no_dark_pixels",
        }

    L_swir = float(np.median(swir_l[dark]))
    L_red = float(np.median(red_l[dark]))
    logger.info("K-T 暗像元: %d 个（阴影/水体剔除 %d 个），SWIR p%.0f=%.3f, "
                "L_swir=%.3f, L_red=%.3f, red/NIR中位=%.3f",
                n_dark, n_shadow_rejected, KT_DARK_PERCENTILE, swir_thr,
                L_swir, L_red, float(np.median(red_nir[dark])))

    # --- Step 3: 对候选 AOT 用 6S 系数反演暗像元地表反射率，匹配 red/SWIR=0.45 ---
    def rho_from_radiance(aot, wavelength, radiance):
        """给定 AOT/波长，用 6S 系数从辐亮度反演地表反射率。"""
        s = sixs_builder_func(aot=aot, wavelength=float(wavelength))
        s.run()
        xa = float(s.outputs.coef_xa)
        xb = float(s.outputs.coef_xb)
        xc = float(s.outputs.coef_xc)
        if not (np.isfinite(xa) and np.isfinite(xb) and np.isfinite(xc)) or xa <= 0:
            return None
        y = xa * radiance - xb
        return y / (1.0 + xc * y)

    aot_candidates = np.arange(0.0, 0.65, 0.02)          # 0 ~ 0.64
    best_aot, best_err, best_ratio = None, None, None
    lut_info = []
    for aot_hypo in aot_candidates:
        try:
            rho_s = rho_from_radiance(aot_hypo, swir_wl, L_swir)
            rho_r = rho_from_radiance(aot_hypo, red_wl, L_red)
        except Exception as exc:
            logger.warning("K-T 6S 系数失败 AOT=%.2f: %s", aot_hypo, exc)
            continue
        if rho_s is None or rho_r is None or rho_s < KT_SWIR_REF_MIN:
            continue
        ratio = rho_r / rho_s
        err = abs(ratio - KT_RED_SWIR_EXPECTED)
        lut_info.append((round(float(aot_hypo), 3), round(float(rho_s), 4),
                         round(float(rho_r), 4), round(float(ratio), 3)))
        if best_err is None or err < best_err:
            best_aot, best_err, best_ratio = float(aot_hypo), err, float(ratio)

    if best_aot is None or best_err > KT_RED_SWIR_TOL:
        logger.warning("K-T 反射率比值无法匹配 FLAASH 期望值 %.2f（最好 %.3f，误差 %.3f），"
                       "反演不可靠，回退初始能见度",
                       KT_RED_SWIR_EXPECTED, best_ratio if best_ratio is not None else float("nan"),
                       best_err if best_err is not None else float("nan"))
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "reason": "ratio_no_match",
            "best_ratio": best_ratio,
        }

    aot_estimated = float(np.clip(best_aot, 0.01, 1.0))
    logger.info("K-T 反演成功: AOT550=%.4f（能见度 %.2f km），red/swir地表反射率比=%.3f"
                "（FLAASH期望 %.2f）", aot_estimated, 3.912 / aot_estimated,
                best_ratio, KT_RED_SWIR_EXPECTED)

    return aot_estimated, {
        "source": "Kaufman-Tanré 2-Band (FLAASH coefficient ratio)",
        "method": "6S coefficient rho retrieval + red/swir=0.45 matching; shadow/water excluded",
        "aot550": aot_estimated,
        "visibility_km": float(3.912 / aot_estimated),
        "swir_band_1based": int(swir_idx + 1),
        "swir_wavelength_um": swir_wl,
        "red_band_1based": int(red_idx + 1),
        "red_wavelength_um": red_wl,
        "nir_band_1based": int(nir_idx + 1),
        "nir_wavelength_um": nir_wl,
        "dark_pixel_count": n_dark,
        "shadow_water_rejected": n_shadow_rejected,
        "dark_swir_percentile": KT_DARK_PERCENTILE,
        "dark_swir_radiance": L_swir,
        "dark_red_radiance": L_red,
        "expected_red_swir_ratio": KT_RED_SWIR_EXPECTED,
        "matched_red_swir_ratio": best_ratio,
        "match_error": best_err,
        "lut_sample": lut_info[::5],  # 抽样记录 AOT→rho_swir/rho_red/ratio
    }


def resolve_aot550_over_water(wavelengths_um, radiance_cube,
                              meta, sixs_builder_func, logger):
    """2-Band Over Water：用近红外与红光暗水体像元估算 AOT550。"""
    nir_idx = _nearest_band(wavelengths_um, 0.865)
    red_idx = _nearest_band(wavelengths_um, 0.660)
    nir_wl = float(wavelengths_um[nir_idx])
    red_wl = float(wavelengths_um[red_idx])

    if not (0.84 <= nir_wl <= 0.89 and 0.64 <= red_wl <= 0.68):
        raise ValueError(
            "2-Band Over Water 缺少标准 NIR (~865 nm) 和 Red (~660 nm) 波段"
        )

    nir_l = radiance_cube[nir_idx].astype(np.float64)
    red_l = radiance_cube[red_idx].astype(np.float64)
    valid = (
        np.isfinite(nir_l) & np.isfinite(red_l)
        & (nir_l > 0.0) & (red_l > 0.0)
        & (nir_l < red_l)
    )
    n_valid = int(valid.sum())
    if n_valid < 1000:
        raise ValueError("有效水体暗像元不足: %d" % n_valid)

    nir_sorted = np.percentile(nir_l[valid], [0.1, 1.0])
    red_sorted = np.percentile(red_l[valid], [0.1, 1.0])
    measured_nir = float(np.median(nir_l[valid & (nir_l <= nir_sorted[1])]))
    measured_red = float(np.median(red_l[valid & (red_l <= red_sorted[1])]))

    candidates = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6]
    lut_nir = []
    lut_red = []
    for aot in candidates:
        try:
            nir_sixs = sixs_builder_func(aot=aot, wavelength=nir_wl)
            red_sixs = sixs_builder_func(aot=aot, wavelength=red_wl)
            nir_sixs.run()
            red_sixs.run()
            lut_nir.append(float(nir_sixs.outputs.pixel_radiance))
            lut_red.append(float(red_sixs.outputs.pixel_radiance))
        except Exception as exc:
            logger.warning("Over-Water 6S 模拟失败 AOT=%.2f: %s", aot, exc)
            lut_nir.append(np.nan)
            lut_red.append(np.nan)

    valid_lut = [
        (aot, nir, red)
        for aot, nir, red in zip(candidates, lut_nir, lut_red)
        if np.isfinite(nir) and np.isfinite(red)
    ]
    if len(valid_lut) < 3:
        raise ValueError("Over-Water 6S 查找表有效点不足")

    lut_aot = np.asarray([item[0] for item in valid_lut])
    lut_nir = np.asarray([item[1] for item in valid_lut])
    lut_red = np.asarray([item[2] for item in valid_lut])
    aot_nir = float(np.interp(measured_nir, lut_nir, lut_aot))
    aot_red = float(np.interp(measured_red, lut_red, lut_aot))
    aot_estimated = float(np.clip((aot_nir + aot_red) / 2.0, 0.01, 1.0))

    logger.info(
        "Over-Water K-T 反演: AOT550=%.4f, 有效像元=%d, "
        "NIR=%.3f, Red=%.3f",
        aot_estimated, n_valid, measured_nir, measured_red
    )
    return aot_estimated, {
        "source": "2-Band Over Water",
        "method": "dark water 865/660 nm forward LUT",
        "aot550": aot_estimated,
        "visibility_km": float(3.912 / aot_estimated),
        "nir_band_1based": int(nir_idx + 1),
        "nir_wavelength_um": nir_wl,
        "red_band_1based": int(red_idx + 1),
        "red_wavelength_um": red_wl,
        "dark_water_pixel_count": n_valid,
        "lut_aot": lut_aot.tolist(),
    }


def resolve_aot550(wavelengths_um, radiance_cube, meta, sixs_builder_func, args, logger):
    """统一入口：None / 2-Band (K-T) / 2-Band Over Water。"""
    method = getattr(args, "aot_method", None) or AOT_METHOD
    method_map = {
        "None": "None",
        "2-Band (K-T)": "KT_DARK_PIXEL",
        "2-Band Over Water": "KT_OVER_WATER",
    }
    method_key = method_map.get(method, method)

    if args.aot is not None:
        aot = float(np.clip(args.aot, AOT_MIN_VALID, AOT_MAX_VALID))
        logger.info("AOT550 手工覆盖: %.4f", aot)
        return aot, {"source": "manual", "method": method, "aot550": aot}

    if method_key == "None":
        aot = float(np.clip(AOT550_MANUAL, AOT_MIN_VALID, AOT_MAX_VALID))
        visibility_km = float(3.912 / aot) if aot > 0 else None
        logger.info("AOT 方法 None：使用固定 AOT550=%.4f，估算能见度=%s km",
                    aot, "%.2f" % visibility_km if visibility_km else "unknown")
        return aot, {
            "source": "manual", "method": "None", "aot550": aot,
            "visibility_km": visibility_km,
        }

    # v2.1.0：气溶胶无法反演时，用「初始能见度」（FLAASH Initial Visibility）作为初始值
    # 参与 6S 校正。该值根据成像时大气情况设置，可用 --visibility 覆盖，默认 40 km（FLAASH 默认）。
    # 与旧版"静默回退"的区别：来源写入 aot_detail.source=initial_visibility_fallback，
    # 并产生 warning，最终记入 used_parameters.json，用户可追溯，绝非无视错误。
    _vis = getattr(args, "visibility", None)
    initial_vis = float(VISIBILITY_MANUAL_KM if _vis is None else _vis)
    if initial_vis <= 0:
        raise ValueError("初始能见度必须为正数，当前 %g km" % initial_vis)

    if method_key in ("KT_DARK_PIXEL", "KT_OVER_WATER"):
        detail = None
        failure_reason = None
        try:
            if method_key == "KT_DARK_PIXEL":
                aot, detail = resolve_aot550_kt(
                    wavelengths_um, radiance_cube, meta, sixs_builder_func, logger
                )
            else:
                aot, detail = resolve_aot550_over_water(
                    wavelengths_um, radiance_cube, meta, sixs_builder_func, logger
                )
        except Exception as exc:
            failure_reason = str(exc)

        # 反演成功且未落到内部兜底 → 直接采用
        used_fallback = str((detail or {}).get("source", "")).startswith("visibility_fallback")
        if failure_reason is None and not used_fallback:
            return aot, detail

        # 反演失败或未满足暗像元条件 → 用初始能见度作为初始值（可追溯、非静默）
        reason = failure_reason or (detail or {}).get("reason", "反演条件不满足")
        aot = float(np.clip(visibility_to_aot550(initial_vis),
                            AOT_MIN_VALID, AOT_MAX_VALID))
        logger.warning(
            "气溶胶反演未能得到可靠结果（%s）；按 FLAASH 约定改用初始能见度 %g km "
            "作为初始值参与大气校正（AOT550=%.4f）。该回退已记入参数文件。",
            reason, initial_vis, aot)
        return aot, {
            "source": "initial_visibility_fallback", "method": method,
            "visibility_km": initial_vis, "aot550": aot,
            "reason": reason, "initial_visibility_user_set":
                getattr(args, "visibility", None) is not None,
        }

    raise RuntimeError("未知的 AOT 方法: %s" % method)



def find_dem_datasets(dem_dir, logger):
    """扫描 DEM 目录，返回可直接读取的栅格（含 ZIP 内栅格）列表。"""
    """列出 DEM 目录下的栅格数据集路径（zip 内 .img 通过 /vsizip/ 读取）。"""
    if not dem_dir or not os.path.isdir(dem_dir):        # 目录未配置或不存在
        return []                                        # 返回空列表
    paths = []                                           # 收集可用栅格路径
    for path in sorted(glob.glob(os.path.join(dem_dir, "*"))):  # 遍历目录条目
        lower = path.lower()                             # 小写便于判后缀
        if lower.endswith(".zip"):                       # 压缩包内可能含栅格
            try:                                         # 读取zip可能失败
                with zipfile.ZipFile(path) as archive:   # 打开压缩包
                    members = [name for name in archive.namelist()  # 挑出栅格成员
                               if name.lower().endswith((".img", ".tif", ".tiff"))  # 栅格后缀
                               and not name.endswith("/")]  # 排除目录项
                if not members:                          # 包内没有栅格
                    logger.debug("DEM 压缩包内没有栅格: %s", path)  # 调试记录
                    continue                             # 换下一个条目
                members.sort()                           # 成员名排序保证稳定
                paths.append("/vsizip/" + os.path.abspath(path).replace("\\", "/")  # 拼vsizip路径
                             + "/" + members[0])         # 取第一个栅格成员
            except Exception as exc:                     # 压缩包损坏等异常
                logger.warning("无法读取 DEM 压缩包: %s (%s)", path, exc)  # 记录告警
        elif lower.endswith((".tif", ".tiff", ".img")):  # 直接是栅格文件
            paths.append(path)                           # 直接加入列表
    return paths                                         # 返回栅格路径列表


def compute_mean_elevation(ring, dem_dir, logger):
    """在足迹范围内按 DEM_RES 网格采样 DEM，返回平均海拔（km）。"""
    """仅统计整景足迹多边形内的有效 DEM 像元并返回高程统计。"""
    started = time.time()                                # 记录起始时刻
    logger.info("开始从 DEM 计算整景足迹平均海拔：%s", dem_dir)  # 打印DEM目录
    paths = find_dem_datasets(dem_dir, logger)           # 扫描DEM栅格
    if not paths:                                        # 没有可用DEM
        logger.warning("未找到 DEM 数据: %s", dem_dir)    # 记录告警
        return None, None                                # 返回无结果
    logger.info("DEM 数据: %d 个瓦片", len(paths))        # 打印瓦片数量

    vrt_path = "/vsimem/zy1e_step2_dem.vrt"              # 内存虚拟VRT路径
    vrt = None                                           # VRT句柄初值
    try:                                                 # 构建与裁剪可能失败
        # v2.2.0：DEM 目录含多个 UTM 带（E113/E114=32649，E115~E118=32650），不能取"首个
        # 瓦片"的 CRS；改由足迹中心推导目标 UTM 带，并只保留与足迹外接矩形相交的瓦片，
        # 再允许 VRT 混合投影，避免异构瓦片导致 warp 全为 nodata。
        lon_mean = float(np.mean([float(lon) for lon, _ in ring]))
        lat_mean = float(np.mean([float(lat) for _, lat in ring]))
        utm_band = int(np.floor((lon_mean + 180.0) / 6.0)) + 1
        epsg = (32600 if lat_mean >= 0 else 32700) + utm_band

        srs_dem = osr.SpatialReference()                 # 目标 UTM 坐标系
        srs_dem.ImportFromEPSG(epsg)
        srs_dem.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        srs_wgs84 = osr.SpatialReference()               # WGS84坐标系对象
        srs_wgs84.ImportFromEPSG(4326)                   # 设为地理坐标系
        srs_wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
        transformer = osr.CoordinateTransformation(srs_wgs84, srs_dem)  # 建坐标转换
        points = [transformer.TransformPoint(lon, lat) for lon, lat in ring]  # 足迹转目标UTM
        min_x = min(point[0] for point in points)         # 西边界
        max_x = max(point[0] for point in points)         # 东边界
        min_y = min(point[1] for point in points)         # 南边界
        max_y = max(point[1] for point in points)         # 北边界

        # 只保留与足迹外接矩形相交的 DEM 瓦片
        overlapping = []
        for path in paths:
            _ds = gdal.Open(path)
            if _ds is None:
                continue
            _srs = osr.SpatialReference()
            _srs.ImportFromWkt(_ds.GetProjection())
            _srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
            _tr = osr.CoordinateTransformation(_srs, srs_dem)
            _gt = _ds.GetGeoTransform()
            _w, _h = _ds.RasterXSize, _ds.RasterYSize
            _cx = [_gt[0], _gt[0] + _gt[1] * _w, _gt[0] + _gt[1] * _w, _gt[0]]
            _cy = [_gt[3], _gt[3], _gt[3] + _gt[5] * _h, _gt[3] + _gt[5] * _h]
            _xs, _ys = [], []
            for _x0, _y0 in zip(_cx, _cy):
                _pt = _tr.TransformPoint(_x0, _y0)
                _xs.append(_pt[0]); _ys.append(_pt[1])
            _ds = None
            if max(_xs) < min_x or min(_xs) > max_x or max(_ys) < min_y or min(_ys) > max_y:
                continue
            overlapping.append(path)
        if not overlapping:
            logger.warning("没有与足迹相交的 DEM 瓦片")
            return None, None
        logger.info("DEM 与足迹相交瓦片: %d / %d", len(overlapping), len(paths))

        vrt = gdal.BuildVRT(
            vrt_path, overlapping,
            options=gdal.BuildVRTOptions(allowProjectionDifference=True))
        if vrt is None:                                  # 构建失败
            raise RuntimeError("DEM VRT 构建失败")        # 报错

        dataset = gdal.Warp(                             # 内存中裁剪重投影
            "", vrt, format="MEM",                       # 输出到内存
            dstSRS=srs_dem.ExportToWkt(),                # 目标坐标系（足迹 UTM 带）
            outputBounds=[min_x, min_y, max_x, max_y],   # 裁剪范围
            xRes=DEM_RES, yRes=DEM_RES,                  # 采样分辨率
            resampleAlg="bilinear",                      # 双线性重采样
            dstNodata=32767,                             # 目标无效值
            multithread=True,                            # 启用多线程
            warpOptions=["NUM_THREADS=ALL_CPUS"])        # 使用全部CPU
        if dataset is None:                              # 裁剪失败
            raise RuntimeError("DEM 范围裁剪失败")        # 报错
        elevation = dataset.GetRasterBand(1).ReadAsArray().astype(np.float64)  # 读高程数组
        output_gt = dataset.GetGeoTransform()            # 读取裁剪后仿射变换
        geo_transform = output_gt                        # v2.2.0：记录实际采样网格 GT
        logger.info("DEM 足迹矩形网格读取完成：%d×%d @ %g m（EPSG:%d）",
                    dataset.RasterXSize, dataset.RasterYSize, DEM_RES, epsg)
        dataset = None                                   # 释放内存数据集
    finally:                                             # 无论成败都清理临时VRT
        vrt = None                                       # 释放VRT句柄
        try:                                             # 删除可能失败
            gdal.Unlink(vrt_path)                        # 删除内存VRT
        except Exception:                                # 忽略删除错误
            pass                                         # 继续退出

    # 将 DEM 像元中心投影到同一坐标系，并排除足迹外接矩形中的足迹外像元。
    ring_dem = [[point[0], point[1]]                     # 足迹转DEM坐标后的环
                for point in (transformer.TransformPoint(lon, lat)  # 逐点转换
                              for lon, lat in ring)]     # 遍历原始足迹
    cols = np.arange(elevation.shape[1], dtype=np.float64) + 0.5  # 列中心坐标
    rows = np.arange(elevation.shape[0], dtype=np.float64) + 0.5  # 行中心坐标
    x_grid = (output_gt[0] + cols[None, :] * output_gt[1]  # 像元中心X坐标
              + rows[:, None] * output_gt[2])            # 叠加行列项
    y_grid = (output_gt[3] + cols[None, :] * output_gt[4]  # 像元中心Y坐标
              + rows[:, None] * output_gt[5])            # 叠加行列项
    footprint_mask = points_in_polygon(x_grid, y_grid, ring_dem)  # 足迹内掩膜

    # 排除常见 DEM 填充值及超出合理地形高程范围的值。
    valid = footprint_mask & np.isfinite(elevation)      # 足迹内且数值有限
    valid &= (elevation != 32767) & (elevation != -32768)  # 排除常见填充值
    valid &= (elevation > -1000) & (elevation < 9000)    # 排除不合理高程
    total = int(footprint_mask.sum())                    # 足迹内像元总数
    count = int(valid.sum())                             # 有效像元数
    if count == 0:                                       # 足迹内没有有效像元
        logger.warning("DEM 整景足迹多边形内没有有效像元")  # 记录原因
        return None, None                                # 返回无结果
    mean_m = float(elevation[valid].mean())              # 有效像元平均高程（m）
    coverage = 100.0 * count / max(total, 1)             # 有效覆盖率（%）
    detail = {                                           # 海拔统计明细
        "source": "ASTGTM DEM (%d tiles)" % len(overlapping),  # 数据来源与瓦片数
        "tiles": overlapping,                            # 实际参与计算的瓦片路径
        "elevation_km": mean_m / 1000.0,                 # 平均海拔（km）
        "mean_m": mean_m,                                # 平均高程（m）
        "std_m": float(elevation[valid].std()),          # 高程标准差
        "min_m": float(elevation[valid].min()),          # 最小高程
        "max_m": float(elevation[valid].max()),          # 最大高程
        "valid_pixels": count,                           # 有效像元数
        "coverage_percent": coverage,                    # 有效覆盖率
        "res_m": DEM_RES,                                # 采样分辨率
        "dem_gt": [float(value) for value in geo_transform[:6]],  # DEM地理变换
    }
    logger.info("DEM 整景足迹平均海拔 = %.1f m（范围 %.1f~%.1f m，有效覆盖 %.1f%%，"  # 打印统计
                "%d 像元 @ %g m）", mean_m, detail["min_m"], detail["max_m"],  # 统计参数
                coverage, count, DEM_RES)                # 覆盖率、像元数与分辨率
    logger.info("DEM 海拔统计完成，用时 %s",                # 记录统计完成
                fmt_elapsed(time.time() - started))      # 打印耗时
    # v2.2.0：附带高程图（裁剪到足迹的 DEM 数组 + GT + EPSG），供高程分块 6S 使用
    return detail, {"elevation": elevation, "gt": output_gt, "epsg": epsg}  # 返回明细与高程图


def resolve_elevation(meta, ring, args, logger):
    """确定目标海拔（km）：命令行 → DEM 平均 → 手工兜底，并记录来源。"""
    """确定目标海拔（km）：命令行 > DEM > 手工兜底。"""
    if args.elev is not None:                            # 命令行直接给了海拔
        logger.info("海拔 = %.3f km（命令行指定）", args.elev)  # 记录来源
        return float(args.elev), {"source": "命令行 --elev",  # 返回海拔与来源
                                  "elevation_km": float(args.elev)}, None  # 无高程图
    if AUTO_DETECT_DEM:                                  # 启用了DEM自动计算
        detail, elev_map = compute_mean_elevation(ring, getattr(args, "dem_dir", DEM_DIR), logger)  # 用DEM算平均海拔
        if detail is not None:                           # 计算成功
            return detail["elevation_km"], detail, elev_map  # 返回海拔、明细与高程图
    logger.warning("DEM 自动计算失败，使用兜底海拔 %.3f km", ELEVATION_KM_FALLBACK)  # 告警
    return float(ELEVATION_KM_FALLBACK), {"source": "ELEVATION_KM_FALLBACK",  # 返回兜底海拔
                                          "elevation_km": float(ELEVATION_KM_FALLBACK)}, None  # 无高程图


def build_native_elevation_bands(source_ds, elev_map, logger):
    """v2.2.0：用 RPC 把原生网格像元反算到地面，采样 DEM 得到逐像元海拔并分档。

    返回 (elev_band_map, band_elevs_km)：
      - elev_band_map: int16 数组 (rows, cols)，值为 0..N-1 的档号；
      - band_elevs_km: 各档代表海拔（km）。
    无法计算时返回 (None, None)，调用方回退单海拔。
    """
    if not elev_map or ELEV_BAND_COUNT <= 1:
        return None, None
    try:
        elevation = elev_map["elevation"]
        gt = elev_map["gt"]
        epsg = int(elev_map["epsg"])
    except Exception:
        return None, None
    n_rows, n_cols = source_ds.RasterYSize, source_ds.RasterXSize
    try:
        transformer = gdal.Transformer(source_ds, None, ["METHOD=RPC"])
    except Exception:
        logger.warning("RPC 变换不可用，跳过逐像元海拔分块")
        return None, None
    srs_utm = osr.SpatialReference(); srs_utm.ImportFromEPSG(epsg)
    srs_utm.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    srs_wgs = osr.SpatialReference(); srs_wgs.ImportFromEPSG(4326)
    srs_wgs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    trans = osr.CoordinateTransformation(srs_wgs, srs_utm)

    step = max(1, max(n_rows, n_cols) // ELEV_BAND_COARSE_STEP)
    ys = np.arange(0, n_rows, step, dtype=np.int64)
    xs = np.arange(0, n_cols, step, dtype=np.int64)
    ny, nx = len(ys), len(xs)
    elev_coarse = np.full((ny, nx), np.nan, dtype=np.float64)
    dh, dw = elevation.shape
    for i in range(ny):
        for j in range(nx):
            col, row = float(xs[j]), float(ys[i])
            ok, pt = transformer.TransformPoint(0, col, row)
            if not ok:
                continue
            x, y, _ = trans.TransformPoint(pt[0], pt[1])
            c = (x - gt[0]) / gt[1]
            r = (gt[3] - y) / (-gt[5]) if gt[5] < 0 else (y - gt[3]) / gt[5]
            c0 = int(np.floor(c)); r0 = int(np.floor(r))
            if 0 <= r0 < dh - 1 and 0 <= c0 < dw - 1:
                fc = c - c0; fr = r - r0
                v = ((elevation[r0, c0] * (1 - fc) + elevation[r0, c0 + 1] * fc) * (1 - fr)
                     + (elevation[r0 + 1, c0] * (1 - fc) + elevation[r0 + 1, c0 + 1] * fc) * fr)
                elev_coarse[i, j] = v
    valid = np.isfinite(elev_coarse)
    if int(valid.sum()) < 50:
        logger.warning("逐像元海拔有效点不足，退回单海拔")
        return None, None
    qs = np.linspace(float(elev_coarse[valid].min()), float(elev_coarse[valid].max()),
                     ELEV_BAND_COUNT)
    qs = np.unique(np.round(qs, 2))                       # 去重（极平坦地形档数可能减少）
    n_bands = len(qs)
    band_coarse = np.full((ny, nx), -1, dtype=np.int16)
    ev = elev_coarse[valid]
    band_coarse[valid] = np.argmin(np.abs(ev[:, None] - qs[None, :]), axis=1)

    band_map = np.zeros((n_rows, n_cols), dtype=np.int16)
    for i in range(ny):
        r0 = int(ys[i]); r1 = int(ys[i + 1]) if i + 1 < ny else n_rows
        for j in range(nx):
            c0 = int(xs[j]); c1 = int(xs[j + 1]) if j + 1 < nx else n_cols
            band_map[r0:r1, c0:c1] = int(band_coarse[i, j])
    band_map = np.clip(band_map, 0, n_bands - 1)          # 足迹外像元(-1)夹到有效档
    logger.info("高程分块 6S：%d 档，代表海拔 %s m", n_bands,
                ", ".join("%.0f" % q for q in qs))
    return band_map, [float(q) / 1000.0 for q in qs]


# =================================================================
# 6S 大气模式、气溶胶及系数计算
# =================================================================

def auto_select_atmos_profile(date, center_lat, logger):
    """按月份与纬度粗略选择 6S 大气模式。"""
    """按成像日期和中心纬度选择 6S 大气模式（南半球月份翻转）。"""
    if center_lat is None:                               # 缺少中心纬度
        logger.warning("无中心纬度，使用 %s", ATMOS_PROFILE_MANUAL)  # 记录改用配置值
        return ATMOS_PROFILE_MANUAL                      # 返回配置区手工模式
    lat_abs = abs(center_lat)                            # 取纬度绝对值
    is_summer = 4 <= date.month <= 9                     # 北半球4~9月为夏季
    if center_lat < 0:                                   # 南半球月份相反
        is_summer = not is_summer                        # 反转季节判断
    if lat_abs < 23.0:                                   # 低纬度为热带
        profile = "Tropical"                             # 热带模式
    elif lat_abs <= 55.0:                                # 中纬度
        profile = "MidlatitudeSummer" if is_summer else "MidlatitudeWinter"  # 夏/冬
    else:                                                # 高纬度
        profile = "SubarcticSummer" if is_summer else "SubarcticWinter"      # 夏/冬
    logger.info("大气模式: 纬度 %.2f°, 月份 %d → %s（自动）",  # 记录选择结果
                lat_abs, date.month, profile)            # 纬度、月份与模式
    return profile                                       # 返回大气模式名


def resolve_atmos_profile(meta, args, logger):
    """确定 6S 大气模式：命令行 → 自动判断 → 配置区手工值。"""
    if args.atmos is not None:                           # 命令行指定了大气模式
        return args.atmos, "命令行 --atmos"               # 直接采用并标注来源
    if WATER_VAPOR_CM is not None:                       # 配置了自定义水汽
        name = "UserWaterAndOzone(water=%.2f, ozone=%.2f)" % (  # 构造自定义模式名
            float(WATER_VAPOR_CM),                       # 水汽含量
            float(OZONE_CM if OZONE_CM is not None else OZONE_CM_DEFAULT))  # 臭氧量
        return name, "配置区 WATER_VAPOR_CM"              # 返回模式名与来源
    if AUTO_DETECT_ATMOS and meta.get("date") is not None:  # 可自动判断且有日期
        return auto_select_atmos_profile(meta["date"], meta.get("center_lat"),  # 自动选模式
                                         logger), "自动（纬度+月份）"  # 标注来源
    return ATMOS_PROFILE_MANUAL, "配置区 ATMOS_PROFILE_MANUAL"  # 退回配置区手工值


def guess_aero_profile(meta, logger):
    """按经纬度/下垫面粗略猜测气溶胶类型。"""
    """无离线土地覆盖数据时的经纬度粗判（默认关闭，仅预留）。"""
    lat = meta.get("center_lat")                         # 中心纬度
    lon = meta.get("center_lon")                         # 中心经度
    if lat is None or lon is None:                       # 缺少中心坐标
        return AERO_PROFILE_MANUAL                       # 退回配置区默认类型
    lat_abs = abs(lat)                                   # 取纬度绝对值
    desert = (75 < lon < 110) and (35 < lat_abs < 45)    # 粗略判定沙漠区
    big_city = ((115 < lon < 120 and 38 < lat_abs < 41)  # 华北城市群
                or (118 < lon < 123 and 29 < lat_abs < 33)  # 长三角城市群
                or (112 < lon < 115 and 21 < lat_abs < 24))  # 珠三角城市群
    if desert:                                           # 落在沙漠区
        profile = "Desert"                               # 沙漠型气溶胶
    elif big_city:                                       # 落在城市群
        profile = "Urban"                                # 城市型气溶胶
    else:                                                # 其他区域
        profile = AERO_PROFILE_MANUAL                    # 使用配置区默认类型
    logger.info("气溶胶类型: 经纬度粗判 → %s", profile)     # 记录判断结果
    return profile                                       # 返回气溶胶类型


def resolve_aero_profile(meta, args, logger):
    """确定 6S 气溶胶类型：命令行 → 自动判断 → 配置区手工值。"""
    if args.aero is not None:                            # 命令行指定了气溶胶类型
        return args.aero, "命令行 --aero"                 # 直接采用并标注来源
    if AUTO_DETECT_AERO:                                 # 启用了自动判断
        return guess_aero_profile(meta, logger), "自动（经纬度粗判）"  # 粗判结果
    return AERO_PROFILE_MANUAL, "配置区 AERO_PROFILE_MANUAL"  # 退回配置区手工值


def resolve_sixs_exe(logger):
    """定位 6S 可执行程序：配置区 → 环境变量 → 脚本 _bin → PATH。"""
    """查找 6S 可执行文件。"""
    # v2.1.0：注入 6S 的 gfortran 运行库目录到 PATH，否则 sixs.exe exit=127、输出 0 行
    if os.path.isdir(SIXS_DLL_DIR) and SIXS_DLL_DIR not in os.environ.get("PATH", ""):
        os.environ["PATH"] = SIXS_DLL_DIR + os.pathsep + os.environ.get("PATH", "")
        logger.info("已注入 6S 运行库目录到 PATH: %s", SIXS_DLL_DIR)
    candidates = []                                      # 候选可执行文件列表
    if SIXS_EXE:                                         # 配置区指定了路径
        candidates.append(SIXS_EXE)                      # 加入候选
    env_path = os.environ.get("SIXS_EXE")                # 读取同名环境变量
    if env_path:                                         # 环境变量有值
        candidates.append(env_path)                      # 加入候选
    # Conda on Windows usually installs executables in the active environment's Library/bin.
    conda_bin = os.path.join(sys.prefix, "Library", "bin")  # 当前环境的bin目录
    for name in ("sixs.exe", "sixs", "sixsV1.1.exe", "sixsV1.1"):  # 常见程序名
        candidates.append(os.path.join(conda_bin, name))  # 环境bin下逐个候选
    script_dir = os.path.dirname(os.path.abspath(__file__))  # 本脚本所在目录
    for name in ("sixs.exe", "sixs", "sixsV1.1.exe", "sixsV1.1"):  # 常见程序名
        candidates.append(os.path.join(script_dir, "_bin", name))  # 脚本_bin下候选
    for name in ("sixs.exe", "sixs", "sixsV1.1.exe", "sixsV1.1"):  # 常见程序名
        found = shutil.which(name)                       # 到PATH中查找
        if found:                                        # PATH中找到了
            candidates.append(found)                     # 加入候选

    for candidate in candidates:                         # 逐个验证候选
        if candidate and os.path.isfile(candidate):      # 文件存在
            logger.info("6S 可执行文件: %s", candidate)   # 记录找到的路径
            return os.path.abspath(candidate)            # 返回绝对路径
    raise FileNotFoundError(                             # 全都找不到则报错
        "未找到 6S 可执行文件 sixs.exe。\n"                # 错误说明
        "请把 sixs.exe（及 libgfortran/libgcc/libquadmath/libwinpthread DLL）放到 %s，\n"  # 放置提示
        "或设置环境变量 SIXS_EXE 指向可执行文件。"         # 环境变量提示
        % os.path.join(script_dir, "_bin"))              # 拼入_bin目录


def build_sixs(sixs_exe, meta, atmos_name, aero_name, aot550,
               elevation_km, logger):
    """构造 Py6S 对象（复用同一个对象逐波段 run，避免重复查找可执行文件）。"""
    try:                                                 # Py6S 为可选依赖
        from Py6S import (SixS, Geometry, AtmosProfile, AeroProfile,  # 导入6S核心类
                          Wavelength, AtmosCorr, GroundReflectance)  # 导入其余组件
    except ImportError as exc:                           # 未安装 Py6S
        raise ImportError(                               # 给出明确安装提示
            "未安装 Py6S，请用装有 GDAL+Py6S 的 Python 运行本脚本"  # 提示信息
            "（例如 E:\\ProgramData\\anaconda3\\envs\\hsi\\python.exe）: %s" % exc)  # 示例解释器

    atmos_map = {                                        # 大气模式名到Py6S枚举
        "Tropical": AtmosProfile.Tropical,               # 热带
        "MidlatitudeSummer": AtmosProfile.MidlatitudeSummer,  # 中纬夏季
        "MidlatitudeWinter": AtmosProfile.MidlatitudeWinter,  # 中纬冬季
        "SubarcticSummer": AtmosProfile.SubarcticSummer,  # 副极地夏季
        "SubarcticWinter": AtmosProfile.SubarcticWinter,  # 副极地冬季
        "USStandard": AtmosProfile.USStandard1962,       # 美国标准大气
    }
    aero_map = {                                         # 气溶胶类型到Py6S枚举
        "Continental": AeroProfile.Continental,          # 大陆型
        "Urban": AeroProfile.Urban,                      # 城市型
        "Maritime": AeroProfile.Maritime,                # 海洋型
        "Desert": AeroProfile.Desert,                    # 沙漠型
    }

    sixs = SixS(sixs_exe)                                # 用指定程序创建6S对象
    sixs.geometry = Geometry.User()                      # 使用自定义几何
    sixs.geometry.solar_z = float(meta["solar_z"])       # 太阳天顶角
    sixs.geometry.solar_a = float(meta["solar_a"])       # 太阳方位角
    sixs.geometry.view_z = float(meta["view_z"])         # 观测天顶角
    sixs.geometry.view_a = float(meta["view_a"])         # 观测方位角
    if meta.get("date") is not None:                     # 有成像日期
        sixs.geometry.day = int(meta["date"].day)        # 设置日
        sixs.geometry.month = int(meta["date"].month)    # 设置月

    if isinstance(atmos_name, str) and atmos_name.startswith("UserWaterAndOzone"):  # 自定义水汽
        sixs.atmos_profile = AtmosProfile.UserWaterAndOzone(  # 用指定水汽与臭氧
            float(WATER_VAPOR_CM),                       # 水汽含量
            float(OZONE_CM if OZONE_CM is not None else OZONE_CM_DEFAULT))  # 臭氧量
    else:                                                # 使用预定义大气模式
        if atmos_name not in atmos_map:                  # 模式名不在映射表
            raise ValueError("不支持的大气模式: %s（可选: %s）"  # 报错并列出可选值
                             % (atmos_name, ", ".join(sorted(atmos_map))))  # 拼参数
        sixs.atmos_profile = AtmosProfile.PredefinedType(atmos_map[atmos_name])  # 设模式

    if aero_name not in aero_map:                        # 类型名不在映射表
        raise ValueError("不支持的气溶胶类型: %s（可选: %s）"  # 报错并列出可选值
                         % (aero_name, ", ".join(sorted(aero_map))))  # 拼参数
    sixs.aero_profile = AeroProfile.PredefinedType(aero_map[aero_name])  # 设气溶胶
    sixs.aot550 = float(aot550)                          # 设置AOT550
    # v2.1.0：星载传感器必须用 satellite_level（写入 -1000，6S 卫星层标记）。
    # 实测 satellite_level 与 90/99 km 系数偏差 0.000%；降到 10 km（入大气层内）才偏 4~25%。
    # 778 km 与 100 km 间为近真空，卫星具体高度不进辐射传输方程，无精度损失。
    # 严禁改回 set_sensor_custom_altitude(778)：Py6S 硬性限制 <100 km，必抛 ValueError。
    sixs.altitudes.set_sensor_satellite_level()
    sixs.altitudes.set_target_custom_altitude(float(elevation_km))  # 目标海拔
    sixs.atmos_corr = AtmosCorr.AtmosCorrLambertianFromRadiance(0)  # 辐亮度反演模式
    # 地面反射率对 Lambertian 大气校正系数没有影响，显式设为 0 便于阅读 6S 输入文件
    sixs.ground_reflectance = GroundReflectance.HomogeneousLambertian(0.0)  # 均匀朗伯面
    return sixs                                          # 返回配置好的6S对象


# 6S（1.1.x 的某些构建）在个别波长会输出 NaN 系数（本数据实测 0.5506/0.6708/
# 1.2408 μm 附近约 1 nm 宽度内）。这些波段的系数与相邻波长几乎一致，因此按
# 微小平移重试，失败时再用相邻有效波段线性插值补齐。
SIXS_WL_OFFSETS = (0.001, -0.001, 0.002, -0.002, 0.003, -0.003,  # 重试波长偏移量
                   0.005, -0.005, -0.008, 0.008)                  # 其余偏移量


def _sixs_coefficients_once(sixs, wavelength):
    """在指定波长跑一次 6S，返回 (xa, xb, xc)。"""
    from Py6S import Wavelength                          # 局部导入波长类
    sixs.wavelength = Wavelength(float(wavelength))      # 设置本次计算波长
    sixs.run()                                           # 执行一次6S计算
    return np.asarray([float(sixs.outputs.coef_xa),      # 取xa系数
                       float(sixs.outputs.coef_xb),      # 取xb系数
                       float(sixs.outputs.coef_xc)], dtype=np.float64)  # 取xc系数


def _coef_from_neighbor_wavelength(sixs, wavelengths, band_index, logger):
    """v2.1.0：取目标波段最近的一个正常邻带系数，作为偏移一致性校验的参照。

    只调用 6S 一次（邻带波长），成本可忽略；若邻带本身算失败则返回 None 跳过校验。
    """
    n = int(len(wavelengths))
    for offset in (-1, 1):                              # 先左后右找邻带
        j = band_index + offset
        if not (0 <= j < n):                            # 越界跳过
            continue
        try:
            values = _sixs_coefficients_once(sixs, float(wavelengths[j]))
        except Exception:
            continue
        if np.all(np.isfinite(values)) and values[0] > 0:
            return values
    return None


def compute_sixs_coefficients(sixs, wavelengths, logger):
    """逐波段运行 6S，返回 (xa,xb,xc) 系数表、版本、耗时与波长扰动记录。

    反射率反演公式（6S 大气校正）：
        y  = xa * L - xb
        ρ  = y / (1 + xc * y)
    其中 L 为表观辐亮度，xa/xb/xc 由 6S 的气体透过率、散射项与球面反照率推出。
    """
    """逐波段运行 6S，返回系数矩阵 (n, 3)、6S 版本、耗时与重试波段信息。"""
    n_bands = int(len(wavelengths))                      # 波段总数
    logger.info("开始计算 6S 系数：共 %d 个波段", n_bands)  # 打印波段数
    coefficients = np.full((n_bands, 3), np.nan, dtype=np.float64)  # 系数表初始化为NaN
    version = None                                       # 6S版本号，稍后读取
    perturbed = {}                                       # 记录被扰动波段的明细
    started = time.time()                                # 记录起始时刻

    for index, wavelength in enumerate(wavelengths, start=1):  # 逐波段（序号从1开始）
        # v2.2.0：水汽强吸收带波段会在反射率域被剔除，无需精确 6S 系数；跳过 6S，
        # 留 NaN 由末段相邻波段插值补齐（这些波段的系数不影响最终成果）。否则 6S 在
        # 水汽带常输出 NaN 且 ±1nm 偏移后系数与邻带差异巨大，触发一致性校验误判。
        if DROP_WATER_VAPOR_BANDS and any(
                lo <= float(wavelength) <= hi for lo, hi in WATER_VAPOR_BANDS_UM):
            continue
        values = None                                    # 本次计算结果
        error = None                                     # 本次异常
        try:                                             # 单波段可能失败
            values = _sixs_coefficients_once(sixs, wavelength)  # 跑一次6S取系数
        except Exception as exc:                         # 记录异常稍后处理
            error = exc                                  # 暂存异常
        if values is not None and np.all(np.isfinite(values)) and values[0] > 0:  # 系数有效
            coefficients[index - 1] = values             # 直接写入系数表
        else:                                            # 系数无效
            # 6S 在个别波长输出 NaN/异常：按微小平移重试（系数随波长变化平缓）
            for offset in SIXS_WL_OFFSETS:               # 依次尝试各偏移量
                candidate = float(wavelength) + offset   # 平移后的波长
                if candidate <= 0.0:                     # 波长须为正
                    continue                             # 非法偏移跳过
                try:                                     # 重试仍可能失败
                    values = _sixs_coefficients_once(sixs, candidate)  # 用平移波长重算
                except Exception:                        # 本次偏移失败
                    continue                             # 换下一个偏移
                if np.all(np.isfinite(values)) and values[0] > 0:  # 系数有效
                    # v2.1.0 偏移一致性校验：微调成功不等于可信。比对偏移后系数与同波段
                    # 原始波长在重试失败前的邻近可用估计；若邻域不平坦（未来可能遇到
                    # 非孤立坏点），相对偏差超阈值即停，避免静默套用平滑假设。
                    # 实测坏点 ±1nm 邻域 xa 偏差 <1%，故 5%/10% 阈值留有余量。
                    neighbor_coef = _coef_from_neighbor_wavelength(
                        sixs, wavelengths, index - 1, logger)
                    if neighbor_coef is not None:        # 有可比对邻带
                        rel_xa = abs(values[0] - neighbor_coef[0]) / max(abs(neighbor_coef[0]), 1e-12)
                        rel_b = abs(values[1] - neighbor_coef[1]) / max(abs(neighbor_coef[1]), 1e-12)
                        rel_c = abs(values[2] - neighbor_coef[2]) / max(abs(neighbor_coef[2]), 1e-12)
                        if (rel_xa > COEF_CONSISTENCY_MAX_REL_XA
                                or rel_b > COEF_CONSISTENCY_MAX_REL_XBXC
                                or rel_c > COEF_CONSISTENCY_MAX_REL_XBXC):
                            raise ValueError(
                                "波段 %d (%.6f μm) 的 6S 系数偏移一致性校验失败："
                                "xa 相对偏差 %.1f%%，xb %.1f%%，xc %.1f%%（阈值 5%%/10%%）。"
                                "邻域不平坦，疑似非孤立坏点，拒绝静默修复。"
                                % (index, float(wavelength), rel_xa * 100,
                                   rel_b * 100, rel_c * 100))
                    coefficients[index - 1] = values     # 写入系数表
                    perturbed[index] = (float(candidate), offset)  # 记录扰动信息
                    logger.warning(                      # 提示发生波长平移
                        "6S 在 %.6f μm 无有效输出，改用 %.6f μm"  # 原波长与新波长
                        "（偏移 %+.1f nm）", float(wavelength), candidate,  # 偏移量
                        offset * 1000.0)                 # 换算为nm
                    break                                # 成功即停止重试
            if np.any(~np.isfinite(coefficients[index - 1])):  # 重试后仍无效
                logger.error("6S 在 %.6f μm 多次重试仍失败%s，稍后用相邻波段插值",  # 记录失败
                             float(wavelength),              # 失败波长
                             "（%s）" % error if error else "")  # 附带异常信息

        if version is None and getattr(sixs.outputs, "version", None):  # 首次取得版本号
            version = sixs.outputs.version               # 记录6S版本
        if index == 1 or index % 10 == 0 or index == n_bands:  # 按进度打印
            elapsed = time.time() - started              # 已用时间
            eta = elapsed / index * (n_bands - index)    # 按平均速度估算剩余
            logger.info("  6S 计算进度 %3d/%d 波段（%.1f%%），已用 %s，预计剩余 %s",  # 进度
                        index, n_bands, 100.0 * index / n_bands,  # 序号与百分比
                        fmt_elapsed(elapsed), fmt_elapsed(eta))  # 已用与剩余时间

    # 仍失败的波段：用相邻有效波段按波长线性插值补齐
    invalid = np.where(~np.isfinite(coefficients).all(axis=1))[0]  # 无效波段下标
    if invalid.size:                                     # 存在无效波段
        valid = np.where(np.isfinite(coefficients).all(axis=1))[0]  # 有效波段下标
        if valid.size < 2:                               # 有效波段太少无法插值
            raise ValueError("6S 系数几乎全部无效（有效波段 %d 个）" % valid.size)  # 报错
        order = np.argsort(np.asarray(wavelengths)[valid])  # 有效波段按波长排序
        sorted_wl = np.asarray(wavelengths)[valid][order]   # 排序后的波长
        for column in range(3):                          # 逐列(xa/xb/xc)插值
            sorted_values = coefficients[valid, column][order]  # 排序后的系数值
            coefficients[invalid, column] = np.interp(    # 线性插值补齐
                np.asarray(wavelengths)[invalid], sorted_wl, sorted_values)  # 插值输入
        logger.warning("6S 系数用相邻有效波段线性插值补齐 %d 个波段: %s",  # 记录补齐情况
                       invalid.size,                     # 补齐波段数
                       ", ".join(str(int(pos) + 1) for pos in invalid))  # 波段序号列表

    if not np.all(np.isfinite(coefficients)):            # 仍存在NaN/Inf
        raise ValueError("6S 系数仍包含 NaN/Inf，结果不可用")  # 结果不可用报错
    if np.any(coefficients[:, 0] <= 0):                  # xa 必须为正
        bad = np.where(coefficients[:, 0] <= 0)[0]       # 找出非法波段
        raise ValueError("6S 系数 xa <= 0（波段 %s），结果不可用"  # 报错
                         % ", ".join(str(int(pos) + 1) for pos in bad[:5]))  # 列出前5个

    if perturbed:                                        # 有波段做过平移重试
        logger.warning("共 %d 个波段使用了微小平移重试: %s", len(perturbed),  # 汇总提示
                       ", ".join("波段%d" % band for band in sorted(perturbed)))  # 波段列表
    elapsed = time.time() - started                      # 总耗时
    logger.info("6S 系数计算完成: %d 波段，用时 %s", n_bands, fmt_elapsed(elapsed))  # 记录完成
    return coefficients, version, elapsed, perturbed     # 返回结果与过程信息


# =================================================================
# 反射率写出与统计
# =================================================================

def build_gtiff_creation_options():
    creation = list(GTIFF_OPTIONS)                       # 复制基础创建选项
    creation.append("INTERLEAVE=BAND")   # 高光谱立方体用 BSQ：逐波段写入效率高
    if GEOTIFF_COMPRESS:                                 # 配置了压缩算法
        creation.extend(["COMPRESS=%s" % GEOTIFF_COMPRESS,  # 压缩方式
                         "PREDICTOR=3",                  # 浮点预测器
                         "ZLEVEL=%d" % GEOTIFF_COMPRESS_LEVEL])  # 压缩等级
    return creation                                      # 返回创建选项列表


def write_reflectance(src_ds, coefficients, wavelengths, tmp_path, logger,
                      elev_band_coeffs=None, elev_band_map=None):
    """逐波段把辐亮度按 6S 系数换算为反射率，写临时 GeoTIFF 并返回统计。

    ρ = (xa*L - xb) / (1 + xc*(xa*L - xb))；无效辐亮度与越界反射率置 NODATA。
    """
    """逐波段把正射辐亮度换算成反射率并写出同网格 GeoTIFF，返回逐波段统计。

    输出直接沿用输入正射辐亮度的仿射变换与投影（CRS/尺寸完全一致），
    不做重采样、不写 RPC；步骤1 已完成正射。
    """
    n_bands, n_rows, n_cols = (src_ds.RasterCount,       # 波段数
                               src_ds.RasterYSize, src_ds.RasterXSize)  # 行数与列数
    if n_bands != int(len(wavelengths)):                 # 波段数与波长数须一致
        raise ValueError("辐亮度波段数 %d != 波长数 %d"     # 报错
                         % (n_bands, len(wavelengths)))  # 拼参数
    if n_bands != coefficients.shape[0]:                 # 波段数与系数行数须一致
        raise ValueError("辐亮度波段数 %d != 6S 系数行数 %d"  # 报错
                         % (n_bands, coefficients.shape[0]))  # 拼参数

    logger.info("开始反射率逐波段处理：%d 波段，影像 %d×%d，临时输出 %s",  # 打印任务信息
                n_bands, n_cols, n_rows, tmp_path)       # 参数
    driver = gdal.GetDriverByName("GTiff")               # 取得GeoTIFF驱动
    if os.path.exists(tmp_path):                         # 临时文件已存在
        os.remove(tmp_path)                              # 先删除避免冲突
    out_ds = driver.Create(tmp_path, n_cols, n_rows, n_bands,  # 创建输出数据集
                           gdal.GDT_Float32,             # 以Float32存储
                           options=build_gtiff_creation_options())  # 应用创建选项
    if out_ds is None:                                   # 创建失败
        raise RuntimeError("无法创建 GeoTIFF: %s" % tmp_path)  # 报错
    # 与输入正射辐亮度严格同网格：直接复制仿射变换与投影（不重采样、不写 RPC）
    out_ds.SetGeoTransform(src_ds.GetGeoTransform())     # 复制仿射变换
    out_ds.SetProjection(src_ds.GetProjection())         # 复制坐标系

    clip_low, clip_high = (REFLECTANCE_CLIP if REFLECTANCE_CLIP else (None, None))  # 裁剪区间
    stats = []                                           # 逐波段统计结果
    masked_pixels = 0                                    # 累计无效像元数
    total_pixels = int(n_rows) * int(n_cols)             # 单波段像元总数

    for index, (wavelength, coef) in enumerate(          # 逐波段处理
            zip(wavelengths, coefficients), start=1):    # 波段号从1开始
        source_band = src_ds.GetRasterBand(index)        # 取输入波段
        out_band = out_ds.GetRasterBand(index)           # 取输出波段
        out_band.SetNoDataValue(NODATA_OUT)              # 设置无效值
        out_band.SetDescription("%.4f um" % float(wavelength))  # 波段描述写波长

        acc = {"count": 0, "sum": 0.0, "sumsq": 0.0,
               "min": None, "max": None, "samples": []}
        blocks_x = (n_cols + 1023) // 1024
        blocks_y = (n_rows + 1023) // 1024
        sample_per_block = max(1, 200000 // max(blocks_x * blocks_y, 1))
        for yoff in range(0, n_rows, 1024):
            height = min(1024, n_rows - yoff)
            for xoff in range(0, n_cols, 1024):
                width = min(1024, n_cols - xoff)
                radiance = source_band.ReadAsArray(
                    xoff=xoff, yoff=yoff, win_xsize=width, win_ysize=height
                ).astype(np.float64)
                invalid = ~np.isfinite(radiance)
                if MASK_INVALID_RADIANCE:
                    invalid |= radiance <= float(RADIANCE_MIN_VALID)

                # v2.2.0：高程分块——逐像元按海拔档取对应 6S 系数；否则单海拔统一系数
                if elev_band_coeffs is not None and elev_band_map is not None:
                    band_blk = elev_band_map[yoff:yoff + height, xoff:xoff + width]
                    coef_blk = elev_band_coeffs[band_blk, index - 1, :]  # (h, w, 3)
                    xa = coef_blk[..., 0]; xb = coef_blk[..., 1]; xc = coef_blk[..., 2]
                else:
                    xa = coef[0]; xb = coef[1]; xc = coef[2]
                y_value = xa * radiance - xb
                denominator = 1.0 + xc * y_value
                with np.errstate(divide="ignore", invalid="ignore"):
                    reflectance = y_value / denominator
                invalid |= ~np.isfinite(reflectance) | (denominator <= 0)
                if clip_low is not None:
                    np.clip(reflectance, clip_low, clip_high, out=reflectance)

                invalid_count = int(invalid.sum())
                masked_pixels += invalid_count
                valid_values = reflectance[~invalid]
                if valid_values.size:
                    acc["count"] += int(valid_values.size)
                    acc["sum"] += float(valid_values.sum(dtype=np.float64))
                    acc["sumsq"] += float(np.square(valid_values, dtype=np.float64).sum())
                    block_min = float(valid_values.min())
                    block_max = float(valid_values.max())
                    acc["min"] = block_min if acc["min"] is None else min(acc["min"], block_min)
                    acc["max"] = block_max if acc["max"] is None else max(acc["max"], block_max)
                    sample = valid_values[::4]
                    if sample.size > sample_per_block:
                        positions = np.linspace(
                            0, sample.size - 1, sample_per_block, dtype=np.int64
                        )
                        sample = sample[positions]
                    acc["samples"].append(sample)

                reflectance[invalid] = NODATA_OUT
                out_band.WriteArray(reflectance.astype(np.float32), xoff=xoff, yoff=yoff)

        if acc["count"]:
            mean = acc["sum"] / acc["count"]
            variance = max(0.0, acc["sumsq"] / acc["count"] - mean * mean)
            sample = np.concatenate(acc["samples"]) if acc["samples"] else np.empty(0)
            p1, p99 = np.percentile(sample, [1, 99]) if sample.size else (mean, mean)
            stats.append({                               # 记录本波段统计
                "band": index,                           # 波段号
                "wavelength": float(wavelength),         # 中心波长
                "min": float(acc["min"]),                # 最小值
                "max": float(acc["max"]),                # 最大值
                "mean": float(mean),                     # 均值
                "std": float(np.sqrt(variance)),         # 标准差
                "p1": float(p1),                         # 1%分位
                "p99": float(p99),                       # 99%分位
                "valid": int(acc["count"]),              # 有效像元数
            })
        else:                                            # 没有有效像元
            logger.warning("  波段 %d 没有有效像元", index)  # 记录告警
            stats.append({"band": index, "wavelength": float(wavelength),  # 补零统计
                          "min": 0.0, "max": 0.0, "mean": 0.0, "std": 0.0,  # 统计值
                          "p1": 0.0, "p99": 0.0, "valid": 0})  # 有效数为0

        if index == 1 or index % 10 == 0 or index == n_bands:  # 按进度打印
            logger.info("  反射率处理进度 %3d/%d 波段（%.1f%%），当前波段有效像元 %.1f%%",  # 进度
                        index, n_bands, 100.0 * index / n_bands,  # 序号与百分比
                        100.0 * stats[-1]["valid"] / max(total_pixels, 1))  # 有效占比

    out_ds.FlushCache()                                  # 刷写缓存到磁盘
    out_ds = None                                        # 释放数据集关闭文件
    if MASK_INVALID_RADIANCE and masked_pixels:          # 有掩膜且确实掩过像元
        logger.info("  辐亮度无效像元累计: %d 个", masked_pixels)  # 汇总无效像元
    return stats                                         # 返回逐波段统计


def detect_uninformative_bands(stats, wavelengths, georef, logger):
    """v2.2.0：大气校正后判定需剔除波段（确定性剔除，保证跨景波段数一致）。

    类别一 water_vapor（确定性剔除）：落在水汽强吸收带的波段按波长直接剔除——
       吸收带内透过率极低、反射率不可信，会污染 SAM/MTMF 与 ML 特征（对齐 FLAASH
       bad bands 惯例）。
    类别二 dn_zero / low_snr（仅评估留痕，不剔除）：来自 Step1 georef.deferred_bad_bands，
       在反射率域用实测统计（有效比例/标准差/动态范围）写入报告供人工复核，但不剔除——
       这些统计在阈值附近随景波动（如 SWIR 尾段 165/166 动态范围略低于阈值），
       若按统计剔除会导致不同景波段数不一致，破坏跨景波段一致性。

    返回 (drop_0based:set, report_rows:list)。report_rows 每行含判定依据供人工复核。
    """
    deferred = (georef or {}).get("deferred_bad_bands", {}) or {}
    dn_zero_0b = set(int(b) for b in deferred.get("dn_zero_bands_0based", []))
    low_snr_0b = set(int(b) for b in deferred.get("low_snr_bands_0based", []))
    stat_candidates = sorted(dn_zero_0b | low_snr_0b)     # 统计判定候选

    # v2.1.2：按波长识别水汽吸收带（不依赖 georef，直接查波长）
    water_vapor_0b = set()
    if DROP_WATER_VAPOR_BANDS:
        for idx, wl in enumerate(wavelengths):
            wl = float(wl)
            for lo, hi in WATER_VAPOR_BANDS_UM:
                if lo <= wl <= hi:
                    water_vapor_0b.add(idx)
                    break

    total_pixels = max((int(item.get("valid", 0)) for item in stats), default=1) or 1
    report_rows = []
    drop = set()
    header = ["band_1based", "wavelength_um", "candidate_group",
              "valid_ratio", "std", "dynamic_range_p99_p1",
              "thr_valid_ratio", "thr_std", "thr_dynamic_range",
              "decision", "reason"]
    report_rows.append(header)

    # --- 类别一：水汽吸收带（按波长剔除，优先标注）---
    for band0 in sorted(water_vapor_0b):
        if band0 >= len(stats):
            continue
        s = stats[band0]
        valid = int(s.get("valid", 0))
        std = float(s.get("std", 0.0))
        dyn = float(s.get("p99", 0.0)) - float(s.get("p1", 0.0))
        group = "water_vapor"
        if band0 in dn_zero_0b:
            group = "water_vapor+dn_zero"
        drop.add(band0)
        report_rows.append([
            band0 + 1, "%.6f" % float(wavelengths[band0]), group,
            "%.6f" % (valid / max(total_pixels, 1)), "%.6f" % std, "%.6f" % dyn,
            "-", "-", "-", "DROP",
            "水汽强吸收带(%.3f-%.3fμm)，透过率极低反射率不可信，按FLAASH bad bands剔除"
            % next((lo, hi) for lo, hi in WATER_VAPOR_BANDS_UM
                   if lo <= float(wavelengths[band0]) <= hi),
        ])

    # --- 类别二：Step1 候选（反射率域统计判定；已在水汽带的不重复）---
    for band0 in stat_candidates:                        # 逐个候选波段
        if band0 < 0 or band0 >= len(stats):             # 越界保护
            logger.warning("候选波段下标 %d 越界（共 %d 波段），跳过", band0, len(stats))
            continue
        if band0 in water_vapor_0b:                      # 已按水汽带剔除，不重复统计行
            continue
        s = stats[band0]                                 # 该波段反射率统计
        valid = int(s.get("valid", 0))                   # 有效像元数
        std = float(s.get("std", 0.0))                   # 标准差
        dyn = float(s.get("p99", 0.0)) - float(s.get("p1", 0.0))  # 动态范围
        valid_ratio = valid / max(total_pixels, 1)

        group = ("dn_zero" if band0 in dn_zero_0b
                 else "low_snr" if band0 in low_snr_0b else "both")
        reasons = []
        if valid_ratio < BAND_DETECT_MIN_VALID_RATIO:
            reasons.append("有效比例<%.0f%%" % (BAND_DETECT_MIN_VALID_RATIO * 100))
        if std < BAND_DETECT_MIN_STD:
            reasons.append("标准差<%g" % BAND_DETECT_MIN_STD)
        if dyn < BAND_DETECT_MIN_DYNAMIC_RANGE:
            reasons.append("动态范围<%g" % BAND_DETECT_MIN_DYNAMIC_RANGE)
        # v2.2.0：仅评估留痕，不剔除，保证跨景波段数一致
        report_rows.append([
            band0 + 1, "%.6f" % float(wavelengths[band0]), group,
            "%.6f" % valid_ratio, "%.6f" % std, "%.6f" % dyn,
            "%.4f" % BAND_DETECT_MIN_VALID_RATIO, "%.6f" % BAND_DETECT_MIN_STD,
            "%.6f" % BAND_DETECT_MIN_DYNAMIC_RANGE,
            "KEEP", ("评估记录，保留" + ("；" + ";".join(reasons) if reasons else "，三项均达标")),
        ])

    logger.info("反射率域波段检测：水汽带 %d 个（按波长剔），Step1候选 %d 个；"
                "合计剔除 %d 个（水汽 %d / DN=0 %d / 低信噪 %d）",
                len(water_vapor_0b), len(stat_candidates), len(drop),
                len(drop & water_vapor_0b), len(drop & dn_zero_0b),
                len(drop & low_snr_0b))
    if drop:
        logger.info("  剔除波段(1-based): %s",
                    ", ".join(str(b + 1) for b in sorted(drop)))
    return drop, report_rows


def write_band_subset(src_tif, keep_0based, wavelengths, dst_tif, logger):
    """v2.1.0：从全波段反射率 tif 抽取保留波段，写出子集 tif（沿用同网格/投影/NoData）。"""
    src = gdal.Open(src_tif)
    if src is None:
        raise RuntimeError("无法打开反射率 tif 以抽取波段: %s" % src_tif)
    keep = sorted(int(i) for i in keep_0based)           # 保持原波段顺序
    if not keep:
        raise ValueError("保留波段为空，拒绝写出空成果")
    driver = gdal.GetDriverByName("GTiff")
    if os.path.exists(dst_tif):
        os.remove(dst_tif)
    out = driver.Create(dst_tif, src.RasterXSize, src.RasterYSize, len(keep),
                        gdal.GDT_Float32, options=build_gtiff_creation_options())
    if out is None:
        raise RuntimeError("无法创建波段子集 tif: %s" % dst_tif)
    out.SetGeoTransform(src.GetGeoTransform())
    out.SetProjection(src.GetProjection())
    metacol = src.GetMetadata("RPC")
    if metacol:
        out.SetMetadata(metacol, "RPC")
    for out_i, src_i in enumerate(keep, start=1):        # 逐波段复制
        src_band = src.GetRasterBand(src_i + 1)
        dst_band = out.GetRasterBand(out_i)
        nodata = src_band.GetNoDataValue()
        if nodata is not None:
            dst_band.SetNoDataValue(nodata)
        desc = src_band.GetDescription() or ("%.4f um" % float(wavelengths[src_i]))
        dst_band.SetDescription(desc)
        arr = src_band.ReadAsArray()
        dst_band.WriteArray(arr)
    out.FlushCache()
    out = None
    src = None
    logger.info("波段子集写出: %s（%d 波段）", dst_tif, len(keep))


def write_detection_report(report_rows, output_path):
    """把波段检测明细写成 CSV。"""
    tmp_path = output_path + ".tmp"
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:
        csv.writer(handle).writerows(report_rows)
    os.replace(tmp_path, output_path)


def write_band_stats_csv(stats, output_path):
    """把逐波段反射率统计写成 CSV。"""
    tmp_path = output_path + ".tmp"                      # 先写临时文件
    with open(tmp_path, "w", newline="", encoding="utf-8") as handle:  # 打开CSV临时文件
        writer = csv.writer(handle)                      # 创建CSV写对象
        writer.writerow(["band", "wavelength_um", "min", "max", "mean",  # 写表头
                         "std", "p1", "p99", "valid_pixels"])  # 其余表头列
        for item in stats:                               # 逐波段写一行
            writer.writerow([                            # 写入该波段统计
                item["band"],                            # 波段号
                "%.6f" % item["wavelength"],             # 波长
                "%.6f" % item["min"],                    # 最小值
                "%.6f" % item["max"],                    # 最大值
                "%.6f" % item["mean"],                   # 均值
                "%.6f" % item["std"],                    # 标准差
                "%.6f" % item["p1"],                     # 1%分位
                "%.6f" % item["p99"],                    # 99%分位
                item["valid"],                           # 有效像元数
            ])
    os.replace(tmp_path, output_path)                    # 原子替换为正式文件


def read_band_stats_csv(path):
    """读取已有 band_stats_reflectance.csv，用于跳过场景时重绘质检图。"""
    stats = []                                           # 解析出的统计列表
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:  # 兼容BOM读取
        for row in csv.DictReader(handle):               # 按字典逐行解析
            try:                                         # 字段可能缺失或非法
                stats.append({                           # 还原统计字典
                    "band": int(row["band"]),            # 波段号
                    "wavelength": float(row["wavelength_um"]),  # 波长
                    "min": float(row["min"]),            # 最小值
                    "max": float(row["max"]),            # 最大值
                    "mean": float(row["mean"]),          # 均值
                    "std": float(row["std"]),            # 标准差
                    "p1": float(row["p1"]),              # 1%分位
                    "p99": float(row["p99"]),            # 99%分位
                    "valid": int(row["valid_pixels"]),   # 有效像元数
                })
            except (KeyError, ValueError):               # 该行为坏数据
                continue                                 # 跳过该行
    return stats                                         # 返回统计列表


# =================================================================
# 质检：带坐标预览与成果图
# =================================================================

def select_qc_bands(wavelengths):
    """选择质检用的 5 个波段（真彩色 / 假彩色 / 短波红外）。"""
    band_map = {                                         # 各用途选取的波段
        "R": _nearest_band(wavelengths, 0.66),           # 红波段
        "G": _nearest_band(wavelengths, 0.55),           # 绿波段
        "B": _nearest_band(wavelengths, 0.47),           # 蓝波段
        "NIR": _nearest_band(wavelengths, 0.86),         # 近红外
        "SWIR": _nearest_band(wavelengths, 2.20),        # 短波红外
    }
    unique = sorted(set(band_map.values()))              # 去重后的波段列表
    return band_map, unique                              # 返回映射与去重列表


def _stretch(array, valid):
    values = array[valid]                                # 取出有效像元
    if values.size == 0:                                 # 没有有效像元
        return np.zeros_like(array)                      # 返回全零拉伸图
    low, high = np.percentile(values, [2, 98])           # 取2%和98%分位
    if high <= low:                                      # 分位异常退化
        high = low + 1e-6                                # 避免除零
    return np.clip((array - low) / (high - low), 0.0, 1.0)  # 线性拉伸到0~1


def prepare_qc_preview(reflectance_path, wavelengths, preview_tif, logger):
    """生成降采样 RGB 预览 tif，供质检图叠加坐标。

    v2.1.1：输入反射率是 Step 1 的 native 网格（无投影、单位阵、仅 RPC），
    gdal.Translate 按 QC_PREVIEW_SCALE 降采样一个近似预览，不生成真实地理坐标
    （真正正射在 Step 3）。_prepare_figure_data 会据投影是否为空自动回退行列轴。
    """
    band_map, unique = select_qc_bands(wavelengths)      # 选定质检波段
    source = gdal.Open(reflectance_path)                 # 打开反射率文件
    if source is None:                                   # 打开失败
        raise RuntimeError("无法打开反射率: %s" % reflectance_path)  # 报错
    try:                                                 # 确保退出时释放句柄
        geo_transform = source.GetGeoTransform()         # 读取仿射变换
        if not geo_transform or abs(geo_transform[1]) <= 0.0:  # 变换无效
            raise RuntimeError("反射率缺少有效仿射变换，无法生成带坐标预览")  # 报错
        preview_res = abs(geo_transform[1]) * QC_PREVIEW_SCALE  # 预览分辨率
        options = gdal.TranslateOptions(                 # 构造Translate参数
            format="GTiff",                              # 输出GeoTIFF
            bandList=[band + 1 for band in unique],      # 选择预览波段
            xRes=preview_res,                            # X方向分辨率
            yRes=preview_res,                            # Y方向分辨率
            resampleAlg="bilinear",                      # 双线性重采样
            noData=NODATA_OUT,                           # 无效值
            outputType=gdal.GDT_Float32,                 # 输出类型
            creationOptions=build_gtiff_creation_options(),  # GeoTIFF选项
        )
        tmp_path = preview_tif + ".tmp.tif"              # 临时预览路径
        if os.path.exists(tmp_path):                     # 临时文件已存在
            os.remove(tmp_path)                          # 先删除
        result = gdal.Translate(tmp_path, source, options=options)  # 执行降采样
        if result is None:                               # 转换失败
            raise RuntimeError("预览降采样失败: %s" % preview_tif)  # 报错
        result = None                                    # 释放结果数据集
        update = gdal.Open(tmp_path, gdal.GA_Update)     # 以更新方式重开
        for position, band_index in enumerate(unique, start=1):  # 逐个预览波段
            band = update.GetRasterBand(position)        # 取该波段
            band.SetNoDataValue(NODATA_OUT)              # 设置无效值
            band.SetDescription("%.4f um" % wavelengths[band_index])  # 写波长描述
        update = None                                    # 关闭更新句柄
        os.replace(tmp_path, preview_tif)                # 原子替换为正式预览
        logger.info("  降采样预览: %s (%d 波段，尺度×%g，源自 native 反射率网格)",
                    preview_tif, len(unique), QC_PREVIEW_SCALE)
        return band_map, unique                          # 返回波段选择结果
    finally:                                             # 无论如何都释放源
        source = None                                    # 关闭源数据集


def _prepare_figure_data(tif_path, preview_tif, band_map, unique,
                         wavelengths, logger):
    """返回 (arrays, extent, xlabel, ylabel, crs_text, size_text, georeferenced)。"""
    if preview_tif and os.path.isfile(preview_tif):      # 带坐标预览存在
        dataset = gdal.Open(preview_tif)                 # 打开预览文件
        if dataset is not None:                          # 打开成功
            arrays = {}                                  # 各用途波段数组
            for name, band_index in band_map.items():    # 逐个用途取数组
                position = unique.index(band_index) + 1  # 预览中对应的波段号
                arrays[name] = dataset.GetRasterBand(position).ReadAsArray(  # 读该波段
                ).astype(np.float64)                     # 转浮点
            width, height = dataset.RasterXSize, dataset.RasterYSize  # 像素宽高
            geo_transform = dataset.GetGeoTransform()    # 读取仿射变换
            projection_wkt = dataset.GetProjection() or ""  # 投影WKT（native 可为空）
            if projection_wkt.strip():                   # 有投影：按米制坐标画 km 轴
                srs = osr.SpatialReference()             # 新建空间参考
                srs.ImportFromWkt(projection_wkt)        # 导入投影（非空，安全）
                epsg = srs.GetAuthorityCode(None) or "?"  # 取EPSG编号
                extent = [geo_transform[0] / 1000.0,     # 西边界（km）
                          (geo_transform[0] + geo_transform[1] * width) / 1000.0,  # 东边界
                          (geo_transform[3] + geo_transform[5] * height) / 1000.0,  # 南边界
                          geo_transform[3] / 1000.0]      # 北边界
                crs_text = "EPSG:%s @ %g m preview" % (epsg, abs(geo_transform[1]))
                dataset = None                           # 关闭数据集
                return (arrays, extent, "Easting (km)", "Northing (km)",  # 坐标轴与范围
                        crs_text, "%d×%d" % (width, height), True)  # 带坐标标志
            # v2.1.1：native 网格无投影——用行列坐标画，避免 ImportFromWkt 空串崩溃
            dataset = None
            crs_text = "native grid (no projection; ortho in Step 3)"
            return (arrays, [0, width, height, 0], "Column", "Row",
                    crs_text, "%d×%d" % (width, height), False)  # 无地理坐标标志

    # 回退：无带坐标预览时直接读正射反射率（降采样），坐标轴只画行列号
    dataset = gdal.Open(tif_path)                        # 打开原始反射率
    if dataset is None:                                  # 打开失败
        return None                                      # 返回空
    width, height = dataset.RasterXSize, dataset.RasterYSize  # 像素宽高
    scale = max(1, int(np.ceil(max(width, height) / float(QC_MAX_SIZE))))  # 降采样倍数
    buf_w, buf_h = max(1, width // scale), max(1, height // scale)  # 降采样后尺寸
    arrays = {}                                          # 各用途波段数组
    for name, band_index in band_map.items():            # 逐个用途取数组
        arrays[name] = dataset.GetRasterBand(band_index + 1).ReadAsArray(  # 读该波段
            buf_xsize=buf_w, buf_ysize=buf_h).astype(np.float64)  # 按降采样尺寸读
    dataset = None                                       # 关闭数据集
    return (arrays, [0, width, height, 0], "Column", "Row",  # 行列坐标与范围
            "no georeference", "%d×%d" % (width, height), False)  # 无坐标标志


def write_quality_figure(tif_path, wavelengths, output_png, scene,
                         params_summary, logger, band_stats=None,
                         preview_tif=None):
    """生成质检成果图：真彩色 / 假彩色 / NDVI / 均值光谱。"""
    try:                                                 # matplotlib 为可选依赖
        import matplotlib                              # 导入绘图库
        matplotlib.use("Agg")                          # 使用无界面后端
        import matplotlib.pyplot as plt                # 导入绘图接口
    except Exception as exc:                           # 导入失败
        logger.warning("matplotlib 不可用，跳过质检图: %s", exc)  # 告警并跳过
        return False                                   # 表示未生成

    wavelengths = np.asarray(wavelengths, dtype=np.float64)  # 波长转浮点数组
    band_map, unique = select_qc_bands(wavelengths)    # 选定质检波段
    prepared = _prepare_figure_data(tif_path, preview_tif, band_map, unique,  # 准备绘图数据
                                    wavelengths, logger)  # 传入波长与日志器
    if prepared is None:                               # 数据准备失败
        logger.warning("无法读取成果影像，跳过质检图: %s", tif_path)  # 告警并跳过
        return False                                   # 表示未生成
    arrays, extent, xlabel, ylabel, crs_text, size_text, georeferenced = prepared  # 解包

    def combine(names):                                # 把三波段合成为RGB
        channels = [arrays[name] for name in names]    # 取出三个通道
        valid = np.ones_like(channels[0], dtype=bool)  # 有效掩膜初值
        for channel in channels:                       # 逐个通道求交
            valid &= np.isfinite(channel) & (channel != NODATA_OUT)  # 有效像元
        rgb = np.ones(channels[0].shape + (3,), dtype=np.float64)  # RGB初值为白
        for index, channel in enumerate(channels):     # 逐通道拉伸
            rgb[..., index] = _stretch(channel, valid)  # 写入拉伸后的通道
        rgb[~valid] = 1.0                              # 无效像元置白
        return np.clip(rgb, 0.0, 1.0), valid           # 返回RGB与有效掩膜

    true_rgb, valid = combine(("R", "G", "B"))         # 合成真彩色
    false_rgb, _ = combine(("NIR", "R", "G"))          # 合成假彩色

    nir = arrays["NIR"]                                # 近红外通道
    red = arrays["R"]                                  # 红通道
    ndvi_valid = (np.isfinite(nir) & np.isfinite(red)  # 两通道均有效
                  & (nir != NODATA_OUT) & (red != NODATA_OUT)  # 均非无效值
                  & ((nir + red) > 0))                 # 分母为正
    with np.errstate(divide="ignore", invalid="ignore"):  # 忽略除零告警
        ndvi = (nir - red) / (nir + red)               # 计算NDVI
    ndvi_map = np.ma.masked_invalid(np.where(ndvi_valid, ndvi, np.nan))  # 掩膜无效值

    if band_stats:                                     # 有逐波段统计
        figure, axes = plt.subplots(2, 2, figsize=(13, 10))  # 2×2画布
        ax_true, ax_false, ax_ndvi, ax_spec = axes.ravel()  # 解包四个子图
    else:                                              # 无统计只画3幅
        figure, axes = plt.subplots(1, 3, figsize=(15, 5.5))  # 1×3画布
        ax_true, ax_false, ax_ndvi = axes.ravel()      # 解包三个子图
        ax_spec = None                                 # 无光谱子图

    wl_r = wavelengths[band_map["R"]]                  # 红波段中心波长
    wl_g = wavelengths[band_map["G"]]                  # 绿波段中心波长
    wl_b = wavelengths[band_map["B"]]                  # 蓝波段中心波长
    wl_nir = wavelengths[band_map["NIR"]]              # 近红外中心波长

    ax_true.imshow(true_rgb, extent=extent, origin="upper",  # 绘制真彩色
                   interpolation="nearest")            # 最近邻插值
    ax_true.set_title("True color  (R %.2f / G %.2f / B %.2f μm)"  # 真彩色标题
                      % (wl_r, wl_g, wl_b), fontsize=10)  # 标题参数
    ax_false.imshow(false_rgb, extent=extent, origin="upper",  # 绘制假彩色
                    interpolation="nearest")           # 最近邻插值
    ax_false.set_title("False color  (NIR %.2f / R %.2f / G %.2f μm)"  # 假彩色标题
                       % (wl_nir, wl_r, wl_g), fontsize=10)  # 标题参数
    ax_ndvi.imshow(ndvi_map, extent=extent, origin="upper", cmap="RdYlGn",  # 绘制NDVI
                   vmin=-0.2, vmax=0.9, interpolation="nearest")  # 色标范围
    ax_ndvi.set_title("NDVI  (NIR %.2f, R %.2f μm)" % (wl_nir, wl_r),  # NDVI标题
                      fontsize=10)                     # 字号
    for axis in (ax_true, ax_false, ax_ndvi):          # 三幅底图统一设置
        axis.set_xlabel(xlabel, fontsize=9)            # X轴标签
        axis.set_ylabel(ylabel, fontsize=9)            # Y轴标签
        axis.tick_params(labelsize=8)                  # 刻度字号

    if ax_spec is not None:                            # 存在光谱子图
        spectrum = sorted(band_stats, key=lambda item: item["wavelength"])  # 按波长排序
        spectrum_wl = np.array([item["wavelength"] for item in spectrum])  # 波长数组
        spectrum_mean = np.array([item["mean"] for item in spectrum])  # 均值数组
        spectrum_std = np.array([item["std"] for item in spectrum])  # 标准差数组
        ax_spec.plot(spectrum_wl, spectrum_mean, color="black", linewidth=1.0)  # 画均值曲线
        ax_spec.fill_between(spectrum_wl,              # 填充标准差带
                             spectrum_mean - spectrum_std,  # 下边界
                             spectrum_mean + spectrum_std,  # 上边界
                             color="0.75", alpha=0.4, linewidth=0)  # 填充样式
        ax_spec.set_xlabel("Wavelength (μm)", fontsize=9)  # X轴标签
        ax_spec.set_ylabel("Mean reflectance", fontsize=9)  # Y轴标签
        ax_spec.set_ylim(bottom=0.0)                   # 纵轴从0开始
        ax_spec.grid(alpha=0.3)                        # 网格淡化
        ax_spec.set_title("Mean reflectance spectrum (full-res stats)",  # 子图标题
                          fontsize=10)                 # 字号

    valid_ratio = float(valid.mean()) if valid.size else 0.0  # 有效像元占比
    figure.suptitle("%s | %d bands | %s | %s\n%s | valid %.1f%%"  # 总标题
                    % (scene, len(wavelengths), size_text, crs_text,  # 标题参数
                       params_summary, 100.0 * valid_ratio), fontsize=10)  # 参数摘要
    figure.tight_layout(rect=(0, 0, 1, 0.94))          # 调整布局留出标题
    figure.savefig(output_png, dpi=QC_DPI, bbox_inches="tight")  # 保存PNG
    plt.close(figure)                                  # 关闭画布释放内存
    logger.info("质检成果图: %s", output_png)            # 记录输出路径
    return True                                        # 表示已生成


def params_summary_text(payload):
    """把参数 JSON 浓缩成质检图上的多行摘要文本。"""
    return ("%s / %s / AOT550=%.3f / elev=%.3f km"  # 摘要格式
            % (payload.get("atmos_profile", "?"),    # 大气模式
               payload.get("aero_profile", "?"),     # 气溶胶类型
               float(payload.get("aot550", 0.0)),    # AOT550
               float(payload.get("elevation_km", 0.0))))  # 海拔


def write_parameters(payload, json_path, txt_path, logger):
    """写 used_parameters.json 与可读文本摘要，记录本次所有参数来源。"""
    atomic_write_json(payload, json_path)                # 原子写出参数JSON
    lines = [                                            # 文本摘要逐行内容
        "# 步骤2 大气校正（6S）参数记录",                  # 标题行
        "场景: %s" % payload.get("scene", ""),            # 场景名
        "正射辐亮度输入: %s" % payload.get("source_radiance", ""),  # 输入文件
        "原始压缩包: %s" % (payload.get("source_tar") or "缺失"),  # 原始压缩包
        "波段数: %d" % payload.get("band_count", 0),      # 波段数
        "输出网格: EPSG:%s @ %s m（与输入正射辐亮度一致，无重采样、无 RPC）"  # 输出网格
        % (payload.get("output_grid_crs_epsg", "?"),     # 网格EPSG
           payload.get("output_grid_res", "?")),         # 网格分辨率
        "成像日期: %s" % (payload.get("date") or "未知"),  # 成像日期
        "传感器: %s / %s" % (payload.get("satellite", "?"),  # 卫星
                             payload.get("sensor", "?")),  # 传感器
        "太阳天顶角: %s" % payload.get("solar_zenith"),   # 太阳天顶角
        "太阳方位角: %s" % payload.get("solar_azimuth"),  # 太阳方位角
        "卫星天顶角: %s" % payload.get("view_zenith"),    # 卫星天顶角
        "卫星方位角: %s" % payload.get("view_azimuth"),   # 卫星方位角
        "中心纬度: %s" % payload.get("center_lat"),       # 中心纬度
        "中心经度: %s" % payload.get("center_lon"),       # 中心经度
        "大气模式: %s（%s）" % (payload.get("atmos_profile"),  # 大气模式
                              payload.get("atmos_source")),  # 模式来源
        "气溶胶类型: %s（%s）" % (payload.get("aero_profile"),  # 气溶胶类型
                               payload.get("aero_source")),  # 类型来源
        "AOT550: %.4f（%s）" % (float(payload.get("aot550", 0.0)),  # AOT550
                              payload.get("aot_source")),  # AOT来源
        "海拔(km): %.4f（%s）" % (float(payload.get("elevation_km", 0.0)),  # 海拔
                                payload.get("elevation_source")),  # 海拔来源
        "水汽/臭氧: %s" % (payload.get("water_vapor") or "6S 大气模式自带廓线"),  # 水汽
        "6S 可执行文件: %s" % payload.get("sixs_path"),   # 6S程序路径
        "6S 版本: %s（Py6S %s）" % (payload.get("sixs_version"),  # 6S版本
                                  payload.get("py6s_version")),  # Py6S版本
        "6S 波长重试: %s" % (                             # 波长重试说明
            ", ".join("波段%d(%+.1f nm)" % (item["band"], item["offset_nm"])  # 拼重试项
                      for item in payload.get("sixs_perturbed_bands") or [])  # 遍历重试
            or "无"),                                    # 无重试时显示无
        "校正公式: %s" % payload.get("formula"),           # 校正公式
        "反射率裁剪: %s" % (payload.get("reflectance_clip") or "不裁剪"),  # 裁剪设置
        "生成时间: %s" % payload.get("created"),          # 生成时间
    ]
    atomic_write_text("\n".join(lines) + "\n", txt_path, )  # 原子写出文本摘要
    logger.info("参数记录: %s / %s", json_path, txt_path)  # 记录两个输出路径


def regenerate_qc(scene_dir, scene, logger):
    """场景已存在时仅补画缺失的 preview_quality.png。"""
    """已有主干成果时，仅补齐质检图/预览。"""
    reflectance = os.path.join(scene_dir, scene + OUT_REFLECTANCE_SUFFIX)  # 反射率文件
    wavelengths_path = os.path.join(scene_dir, OUT_WAVELENGTHS_NAME)  # 波长文件
    stats_path = os.path.join(scene_dir, OUT_STATS_NAME)  # 统计CSV
    qc_path = os.path.join(scene_dir, OUT_QC_NAME)        # 质检图路径
    preview_path = os.path.join(scene_dir, OUT_PREVIEW_GEO_NAME)  # 带坐标预览
    try:                                                 # 补齐失败不影响主流程
        wavelengths = np.load(wavelengths_path)          # 读取波长数组
        band_stats = read_band_stats_csv(stats_path) if os.path.isfile(stats_path) else None  # 读统计
        params = {}                                      # 参数字典初值
        params_path = os.path.join(scene_dir, OUT_PARAMS_JSON)  # 参数JSON路径
        if os.path.isfile(params_path):                  # 参数文件存在
            params = load_georef_json(params_path)       # 载入参数JSON
        preview_tif = None                               # 预览路径初值
        if QC_GEOREF_PREVIEW:                            # 启用了带坐标预览
            try:                                         # 预览生成可能失败
                prepare_qc_preview(reflectance, wavelengths, preview_path, logger)  # 生成预览
                preview_tif = preview_path               # 记录预览路径
            except Exception:                            # 预览失败降级
                logger.warning("补齐带坐标预览失败，质检图回退无坐标降采样",  # 记录告警
                               exc_info=True)            # 附带堆栈
        write_quality_figure(reflectance, wavelengths, qc_path, scene,  # 生成质检图
                             params_summary_text(params), logger,  # 参数摘要
                             band_stats=band_stats, preview_tif=preview_tif)  # 统计与预览
    except Exception as exc:                             # 质检图必需，补齐失败即停（v3）
        logger.error("补齐质检图失败: %s", scene_dir, exc_info=True)  # 记录错误
        raise RuntimeError("质检图补齐失败，用户需要质检图，拒绝静默跳过: %s" % scene_dir) from exc


# =================================================================
# 单景处理
# =================================================================

def resolve_py6s_version():
    try:                                                 # 导入失败视为未知版本
        import Py6S                                     # 导入Py6S包
        return getattr(Py6S, "__version__", None) or \
            getattr(getattr(Py6S, "SixS", None), "__version__", None)  # 包级优先，否则类级
    except Exception:                                    # 导入异常
        return None                                      # 返回未知


def process_scene(scene_dir, out_root, args, logger):
    """处理一景：读步骤1 成果 → 6S 系数 → 逐波段校正 → 参数/统计/质检。"""
    scene = os.path.basename(os.path.normpath(scene_dir))  # 由目录名取场景名
    in_radiance = resolve_radiance_path(scene_dir)       # 定位步骤1辐亮度
    in_wavelengths = os.path.join(scene_dir, IN_WAVELENGTHS_NAME)  # 波长文件路径
    in_georef = os.path.join(scene_dir, IN_GEOREF_NAME)  # 几何元数据路径

    out_dir = os.path.join(out_root, scene)              # 本景输出目录
    final_tif = os.path.join(out_dir, scene + OUT_REFLECTANCE_SUFFIX)  # 反射率成果
    legacy_tifs = [os.path.join(out_dir, name)           # 旧命名成果路径
                   for name in LEGACY_REFLECTANCE_NAMES]  # 遍历旧名列表
    final_wl = os.path.join(out_dir, OUT_WAVELENGTHS_NAME)  # 输出波长文件
    final_georef = os.path.join(out_dir, OUT_GEOREF_NAME)  # 输出元数据文件
    final_coef = os.path.join(out_dir, OUT_COEF_NAME)    # 输出系数文件
    final_stats = os.path.join(out_dir, OUT_STATS_NAME)  # 输出统计CSV
    final_params_json = os.path.join(out_dir, OUT_PARAMS_JSON)  # 输出参数JSON
    final_params_txt = os.path.join(out_dir, OUT_PARAMS_TXT)  # 输出参数文本
    final_qc = os.path.join(out_dir, OUT_QC_NAME)        # 输出质检图
    final_preview = os.path.join(out_dir, OUT_PREVIEW_GEO_NAME)  # 输出预览tif

    # v2.2.0：不再"旧文件名→新文件名"改名复用（会把旧代 142/158 波段成果冒充新成果）。
    # 仅当最终成果自洽（tif 波段数 == 参数里的 band_count）且其来源波段数等于当前输入
    # 波段数时才跳过；否则强制重算，防止旧代成果冒充新成果。
    have_main = all(os.path.isfile(path)                 # 主成果是否齐全
                    for path in (final_tif, final_coef, final_params_json))  # 三项检查
    if not args.overwrite and have_main:
        try:
            _input_wl = np.load(in_wavelengths, allow_pickle=False).ravel()
            _ds_final = gdal.Open(final_tif)
            _params = None
            try:
                with open(final_params_json, "r", encoding="utf-8") as _fh:
                    _params = json.load(_fh)
            except Exception:
                _params = None
            _final_bands = _ds_final.RasterCount if _ds_final is not None else -1
            _rec_bands = int(_params.get("band_count", -1)) if _params else -1
            _rec_src = int(_params.get("source_band_count", -1)) if _params else -1
            if _final_bands != _rec_bands or _rec_src != int(_input_wl.size):
                logger.warning("已有成果与当前输入不一致（tif=%s, 参数 band=%s src=%s, 输入=%d），强制重算",
                               _final_bands, _rec_bands, _rec_src, int(_input_wl.size))
                have_main = False
            _ds_final = None
        except Exception as _exc:
            logger.warning("校验已有成果波段一致性失败，按需重算: %s", _exc)
            have_main = False
    if not args.overwrite and have_main:                 # 无需重做且成果完整
        # 补齐从步骤1 复制的光谱/几何小文件，避免重复计算
        try:                                             # 复制失败不影响跳过
            if not os.path.isfile(final_wl) and os.path.isfile(in_wavelengths):  # 缺波长
                copy_file_atomic(in_wavelengths, final_wl)  # 补复制波长文件
            if not os.path.isfile(final_georef) and os.path.isfile(in_georef):  # 缺元数据
                copy_file_atomic(in_georef, final_georef)  # 补复制元数据
        except OSError:                                  # 复制异常
            logger.warning("补齐 wavelengths_clean.npy/georef.json 失败", exc_info=True)  # 告警
        if MAKE_QC and not args.no_qc and not os.path.isfile(final_qc):  # 缺质检图
            regenerate_qc(out_dir, scene, logger)        # 仅补画质检图
        logger.info("跳过（成果已存在）: %s", out_dir)     # 记录跳过
        return {"scene": scene, "status": "跳过", "bands": None,  # 返回跳过结果
                "elapsed": 0.0, "message": "成果已存在"}  # 耗时与说明

    started = time.time()                                # 记录本景起始时刻
    current_stage = 1                                    # 当前阶段号
    stage_description = "校验步骤1输入与正射网格"          # 阶段说明
    created_dir = not os.path.isdir(out_dir)             # 输出目录是否需新建
    os.makedirs(out_dir, exist_ok=True)                  # 确保输出目录存在
    tmp_tif = os.path.join(out_dir, scene + "_cal_badband_ortho_atmcorr.tmp.tif")  # 临时反射率
    tmp_subset = os.path.join(out_dir, scene + "_band_subset.tmp.tif")  # 波段剔除后的临时反射率
    tmp_coef = os.path.join(out_dir, "py6s_coefficients.tmp.npy")  # 临时系数
    tmp_wl = os.path.join(out_dir, "wavelengths_clean.tmp.npy")  # 临时波长
    final_band_report = os.path.join(out_dir, OUT_BAND_DETECT_REPORT)  # 波段检测报告
    final_repair_log = os.path.join(out_dir, OUT_COEF_REPAIR_LOG)  # NaN 修复台账
    written_targets = []                                 # 已写出的目标文件
    main_complete = False                                # 主成果是否完成
    source_ds = None                                     # 输入数据集句柄
    try:                                                 # 单景异常统一收尾
        # ① 校验步骤1 输入（辐亮度/波长/波段数），输出沿用同一正射网格、不重采样
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        if in_radiance is None:                          # 找不到辐亮度
            raise FileNotFoundError(                     # 报错
                "找不到步骤1正射辐亮度（*_cal_badband_ortho.tif 或旧格式 %s）: %s"  # 提示
                % (LEGACY_RADIANCE_NAME, scene_dir))     # 旧名与目录
        if not os.path.isfile(in_wavelengths):           # 缺少波长文件
            raise FileNotFoundError("找不到步骤1 波长文件: %s" % in_wavelengths)  # 报错

        wavelengths = np.load(in_wavelengths, allow_pickle=False)  # 读取波长数组
        wavelengths = np.asarray(wavelengths, dtype=np.float64).ravel()  # 转浮点一维
        if wavelengths.size == 0 or not np.all(np.isfinite(wavelengths)) \
                or np.any(wavelengths <= 0):             # 空/非有限/非正波长均异常
            raise ValueError("波长数据异常: %s" % in_wavelengths)  # 报错
        original_wavelengths_full = wavelengths.copy()   # v2.2.0：剔除前全波段波长，供修复台账/georef 追溯

        # georef.json 仅作元数据来源（足迹/源 tar/目标 CRS/分辨率），本步不做几何校正
        georef = load_georef_json(in_georef) if os.path.isfile(in_georef) else None  # 读元数据
        if georef is None:                               # 缺元数据
            logger.warning("缺少 georef.json，足迹将回退产品 XML；"  # 告警前半句
                           "不影响正射反射率网格（以栅格自带地理参考为准）")  # 告警后半句

        source_ds = gdal.Open(in_radiance)               # 打开正射辐亮度
        if source_ds is None:                            # 打开失败
            raise RuntimeError("无法打开正射辐亮度: %s" % in_radiance)  # 报错
        if source_ds.RasterCount != wavelengths.size:    # 波段数与波长数不符
            raise ValueError("正射辐亮度波段数 %d != 波长数 %d"  # 报错
                             % (source_ds.RasterCount, wavelengths.size))  # 拼参数
        _raster_size = (georef or {}).get("raster_size") or (georef or {}).get("native_grid", {}).get("raster_size")
        if _raster_size:                                 # 元数据记录了尺寸
            expected = [int(_raster_size[0]), int(_raster_size[1])]  # 记录尺寸
            actual = [source_ds.RasterXSize, source_ds.RasterYSize]  # 实际尺寸
            if expected != actual:                       # 尺寸不一致
                logger.warning("正射辐亮度尺寸 %s 与 georef.json 记录 %s 不一致",  # 告警
                               actual, expected)         # 实际与记录

        # v2.1.1 修复：Step 1 输出 native 网格（无投影、单位阵 GeoTransform、仅 RPC），
        # 正射按设计推迟到 Step 3。因此输入无投影是设计预期，不是错误；直接继承该网格，
        # 不对空投影调用 ImportFromWkt（否则 OGR Error: Corrupt data），也不做 CRS/分辨率校验。
        source_gt = source_ds.GetGeoTransform()          # 读取仿射变换
        projection_wkt = source_ds.GetProjection() or ""  # 投影WKT（native 网格为空）
        is_native_grid = (not projection_wkt.strip())    # 空投影 → native 网格
        source_epsg = None                               # EPSG 初值
        source_res = None                                # 分辨率(m)初值
        if is_native_grid:                               # native 网格：设计预期，跳过校验
            gtf_identity = (source_gt and abs(source_gt[1] - 1.0) < 1e-9
                            and abs(source_gt[5] - 1.0) < 1e-9)
            if not gtf_identity:                         # 空投影却非单位阵，来源异常即停
                raise ValueError("辐亮度无投影但 GeoTransform 非单位阵（%s），"
                                 "来源不符合 native 网格约定" % source_gt)
            logger.info("输入为 native 网格（无投影/单位阵/仅 RPC，符合 v2 设计；"
                        "正射在 Step 3 统一完成），本步原样继承、不重采样")  # 记录
        else:                                            # 有投影：正常读取并校验
            source_srs = osr.SpatialReference()          # 新建空间参考
            source_srs.ImportFromWkt(projection_wkt)     # 导入投影（非空，安全）
            source_epsg = source_srs.GetAuthorityCode(None)  # 取EPSG编号
            source_res = abs(source_gt[1]) if source_gt and source_gt[1] else None  # 像元大小
            if source_epsg is None or source_res is None:  # 有投影但缺编号/分辨率，异常即停
                raise ValueError("辐亮度有投影但缺 EPSG 编号或有效分辨率，"
                                 "无法确认网格（EPSG=%s, res=%s）" % (source_epsg, source_res))
            logger.info("正射网格: EPSG:%s @ %g m, %d×%d（不做重采样）",  # 打印网格
                        source_epsg, source_res,         # EPSG与分辨率
                        source_ds.RasterXSize, source_ds.RasterYSize)  # 行列数
            if georef and georef.get("target_crs_epsg") is not None:  # 元数据有EPSG
                try:                                     # 数值转换可能失败
                    if int(georef["target_crs_epsg"]) != int(source_epsg):  # EPSG不一致
                        logger.warning("辐亮度 CRS (EPSG:%s) 与 georef.json "  # 告警
                                       "target_crs_epsg=%s 不一致，以栅格为准",  # 说明
                                       source_epsg, georef["target_crs_epsg"])  # 两值
                except (TypeError, ValueError):          # 非法值
                    pass                                 # 忽略
            if georef and georef.get("target_res") is not None:  # 元数据有分辨率
                try:                                     # 数值转换可能失败
                    if abs(float(georef["target_res"]) - source_res) > 1e-6:  # 分辨率不一致
                        logger.warning("辐亮度分辨率 %g m 与 georef.json "  # 告警
                                       "target_res=%s 不一致，以栅格为准",  # 说明
                                       source_res, georef["target_res"])  # 两值
                except (TypeError, ValueError):          # 非法值
                    pass                                 # 忽略

        # 元数据展示用：native 网格无自身 CRS/res，用 georef 目标值标注（说明正射待 Step 3）
        meta_epsg = source_epsg
        meta_res = source_res
        if is_native_grid and georef:                    # native 网格用目标值做展示
            meta_epsg = georef.get("target_crs_epsg")    # 展示目标EPSG
            meta_res = georef.get("target_res")          # 展示目标分辨率

        logger.info("正射辐亮度: %d 波段, %d×%d, 波长 %.4f~%.4f μm",  # 打印输入信息
                    source_ds.RasterCount, source_ds.RasterXSize,  # 波段与列
                    source_ds.RasterYSize, wavelengths.min(), wavelengths.max())  # 行与波长
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ② 读取观测几何、日期、中心位置和足迹
        current_stage = 2                                # 进入第2阶段
        stage_description = "解析成像时间、观测几何和场景足迹"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        meta = resolve_metadata(scene, scene_dir, georef, getattr(args, "tar_dir", TAR_DIR), logger)  # 汇总元数据
        ring = footprint_ring(meta["footprint"])         # 生成足迹闭合环
        logger.info("观测几何: 太阳 %.2f°/%.2f°, 卫星 %.2f°/%.2f°, 日期 %s",  # 打印几何
                    meta["solar_z"], meta["solar_a"], meta["view_z"],  # 太阳与观测角
                    meta["view_a"], meta["date"].strftime("%Y-%m-%d"))  # 方位角与日期
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ③ 确定 AOT550、海拔、大气模式和气溶胶类型
        current_stage = 3                                # 进入第3阶段
        stage_description = "确定 AOT550、DEM 海拔及大气/气溶胶模式"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        elevation_km, elevation_detail, elev_map = resolve_elevation(meta, ring, args, logger)  # 确定海拔
        atmos_profile, atmos_source = resolve_atmos_profile(meta, args, logger)  # 确定大气模式
        aero_profile, aero_source = resolve_aero_profile(meta, args, logger)  # 确定气溶胶
        aot_method = args.aot_method or AOT_METHOD
        # v2.1.2：规范化方法名（兼容 CLI 显示值 "2-Band (K-T)" 与配置内部值 "KT_DARK_PIXEL"）
        _method_key = {"2-Band (K-T)": "KT_DARK_PIXEL",
                       "2-Band Over Water": "KT_OVER_WATER"}.get(aot_method, aot_method)
        if args.aot is not None or _method_key == "None":
            aot_wavelengths = np.asarray([], dtype=np.float64)
            radiance_cube = None
        else:
            # K-T 需 3 波段（SWIR+Red+NIR 阴影判据）；Over-Water 用 NIR+Red
            target_wavelengths = (
                (2.175, 0.660, KT_NIR_TARGET_UM) if _method_key == "KT_DARK_PIXEL"
                else (0.865, 0.660)
            )
            aot_band_indices = sorted({
                _nearest_band(wavelengths, target) for target in target_wavelengths
            })
            aot_wavelengths = wavelengths[aot_band_indices]
            radiance_cube = np.stack([
                source_ds.GetRasterBand(index + 1).ReadAsArray().astype(np.float32)
                for index in aot_band_indices
            ])
            logger.info("AOT 反演读取 %d 个波段: %s", len(aot_band_indices),
                        ", ".join("%.4fμm" % w for w in aot_wavelengths))
        sixs_exe = resolve_sixs_exe(logger)

        def kt_sixs_builder(aot, wavelength):
            sixs = build_sixs(
                sixs_exe, meta, atmos_profile, aero_profile,
                float(aot), elevation_km, logger
            )
            from Py6S import Wavelength
            sixs.wavelength = Wavelength(float(wavelength))
            return sixs

        aot550, aot_detail = resolve_aot550(
            aot_wavelengths, radiance_cube, meta, kt_sixs_builder, args, logger
        )
        logger.info("6S 参数: 大气=%s, 气溶胶=%s, AOT550=%.4f, 海拔=%.4f km",  # 打印参数
                    atmos_profile, aero_profile, aot550, elevation_km)  # 参数值
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ④ 逐波段计算 6S 系数 (xa, xb, xc)
        current_stage = 4                                # 进入第4阶段
        stage_description = "运行 6S 并计算各波段校正系数"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        # v2.2.0：高程分块——用 DEM 得到逐像元海拔档，按档分别算 6S 系数（山区精度提升）
        elev_band_map, band_elevs_km = build_native_elevation_bands(source_ds, elev_map, logger)
        elev_band_coeffs = None
        sixs_version = None
        coef_elapsed = 0.0
        perturbed = {}
        if elev_band_map is not None and band_elevs_km:
            elev_band_coeffs_list = []
            for band_elev in band_elevs_km:
                sixs_eb = build_sixs(sixs_exe, meta, atmos_profile, aero_profile,
                                     aot550, band_elev, logger)
                coef_eb, sixs_version, band_elapsed, band_perturbed = \
                    compute_sixs_coefficients(sixs_eb, wavelengths, logger)
                elev_band_coeffs_list.append(coef_eb)
                coef_elapsed += band_elapsed
                perturbed.update(band_perturbed)
            elev_band_coeffs = np.stack(elev_band_coeffs_list)  # (n_elev, n_bands, 3)
            # 代表系数取最接近平均海拔的那一档，供元数据/子集/统计使用
            nearest = int(np.argmin(np.abs(np.asarray(band_elevs_km) - elevation_km)))
            coefficients = elev_band_coeffs[nearest].copy()
        else:
            sixs = build_sixs(sixs_exe, meta, atmos_profile, aero_profile,  # 构造6S对象
                              aot550, elevation_km, logger)  # 传入AOT与海拔
            coefficients, sixs_version, coef_elapsed, perturbed = \
                compute_sixs_coefficients(sixs, wavelengths, logger)  # 逐波段跑6S取系数
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ⑤ 逐波段反演反射率 rho=(xa*L-xb)/(1+xc*(xa*L-xb))，保持输入网格且不写 RPC
        current_stage = 5                                # 进入第5阶段
        stage_description = "逐波段计算并写出地表反射率"     # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        stats = write_reflectance(source_ds, coefficients, wavelengths,  # 写反射率
                                  tmp_tif, logger,       # 输出到临时文件
                                  elev_band_coeffs=elev_band_coeffs,
                                  elev_band_map=elev_band_map)  # 高程分块
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ⑤b v2.1.0 核心：大气校正后在反射率域检测候选波段，剔除真正无信息的波段
        current_stage = 6                                # 进入第6阶段
        stage_description = "反射率域检测并剔除无用波段"   # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时（共8阶段）
                                        stage_description)  # 阶段说明
        drop_0b, report_rows = detect_uninformative_bands(  # 反射率域检测
            stats, wavelengths, georef, logger)          # 返回剔除下标与报告行
        write_detection_report(report_rows, final_band_report)  # 写出检测报告（可复核）
        written_targets.append(final_band_report)         # 纳入失败清理清单
        original_band_count = int(len(wavelengths))       # 记录剔除前波段数

        drop_1based = [int(b) + 1 for b in sorted(drop_0b)]  # 1-based 便于人工核对
        # v2.1.2：从报告行解析水汽带剔除清单（group 含 water_vapor），供元数据/参数分开标注
        water_vapor_dropped_1b = [
            int(row[0]) for row in report_rows[1:]
            if len(row) > 2 and "water_vapor" in str(row[2]) and str(row[9]) == "DROP"
        ]
        stat_dropped_1b = [b for b in drop_1based if b not in water_vapor_dropped_1b]
        if drop_0b:                                      # 有需要剔除的波段
            keep_0b = [i for i in range(original_band_count) if i not in drop_0b]  # 保留下标
            write_band_subset(tmp_tif, keep_0b, wavelengths, tmp_subset, logger)  # 写出子集
            os.replace(tmp_subset, tmp_tif)              # 原子替换为剔除后成果
            wavelengths = wavelengths[keep_0b]           # 同步波长数组
            coefficients = coefficients[keep_0b, :]      # 同步6S系数数组
            stats = [stats[i] for i in keep_0b]          # 同步波段统计
            logger.info("反射率域剔除 %d 个波段（1-based: %s），保留 %d 波段",
                        len(drop_1based),
                        ", ".join(str(b) for b in drop_1based), int(len(wavelengths)))
        else:                                            # 无波段被剔除
            logger.info("反射率域检测：候选波段均保留，成果仍为 %d 波段", original_band_count)
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ⑦ 写入成果元数据、参数记录和波段统计
        current_stage = 7                                # 进入第7阶段
        stage_description = "写入影像元数据、参数记录和波段统计"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        update_ds = gdal.Open(tmp_tif, gdal.GA_Update)   # 以更新方式打开临时成果
        update_ds.SetMetadataItem("SOURCE_RADIANCE", os.path.basename(in_radiance))  # 输入文件名
        update_ds.SetMetadataItem("SOURCE_TAR",          # 原始压缩包名
                                  os.path.basename(meta["tar"])  # 有tar取文件名
                                  if meta["tar"] else "MISSING")  # 无tar标记缺失
        update_ds.SetMetadataItem("BAND_COUNT", str(int(len(wavelengths))))  # 剔除后波段数
        update_ds.SetMetadataItem("SOURCE_BAND_COUNT", str(original_band_count))  # 剔除前波段数
        update_ds.SetMetadataItem("DROPPED_BANDS_1BASED",  # 反射率域剔除的波段
                                  ",".join(str(b) for b in drop_1based))
        update_ds.SetMetadataItem("DROPPED_WATER_VAPOR_1BASED",  # 水汽带剔除（v2.1.2）
                                  ",".join(str(b) for b in water_vapor_dropped_1b))
        update_ds.SetMetadataItem(                       # 写入成像日期
            "IMAGING_DATE", meta["date"].strftime("%Y-%m-%d"))  # 日期字符串
        update_ds.SetMetadataItem("ATMOSPHERIC_CORRECTION", "6S (Py6S)")  # 校正方法
        update_ds.SetMetadataItem("SIXS_EXE", os.path.basename(sixs_exe))  # 6S程序名
        update_ds.SetMetadataItem("SIXS_VERSION", str(sixs_version or "unknown"))  # 6S版本
        update_ds.SetMetadataItem("ATMOS_PROFILE", str(atmos_profile))  # 大气模式
        update_ds.SetMetadataItem("AERO_PROFILE", str(aero_profile))  # 气溶胶类型
        update_ds.SetMetadataItem("AOT550", "%.6f" % aot550)  # AOT550
        update_ds.SetMetadataItem("ELEVATION_KM", "%.6f" % elevation_km)  # 海拔
        update_ds.SetMetadataItem("AOT_SOURCE", str(aot_detail.get("source", "")))  # AOT来源
        update_ds.SetMetadataItem("ELEVATION_SOURCE",    # 海拔来源
                                  str(elevation_detail.get("source", "")))  # 来源文本
        update_ds.SetMetadataItem(                       # 写入校正公式
            "REFLECTANCE_FORMULA",                       # 元数据键
            "rho = (xa*L - xb) / (1 + xc*(xa*L - xb))")  # 公式文本
        update_ds.SetMetadataItem(                       # 写入几何说明
            "GEOMETRY",                                  # 元数据键
            "native grid; RPC orthorectification deferred to step 3")  # 说明
        res_text = ("%g" % meta_res) if meta_res else "?"  # 分辨率文本
        grid_desc = ("native grid (no projection; RPC ortho in Step 3; "
                     "target EPSG:%s @ %s m)" % (meta_epsg or "?", res_text) \
                     if is_native_grid else \
                     "EPSG:%s @ %s m (inherited from %s)" \
                     % (source_epsg or "?", res_text, IN_RADIANCE_SUFFIX))
        update_ds.SetMetadataItem("OUTPUT_GRID", grid_desc)  # 写入输出网格描述
        update_ds.SetMetadataItem("PY6S_COEFFICIENTS", OUT_COEF_NAME)  # 系数文件名
        rpc_metadata = source_ds.GetMetadata("RPC")
        if rpc_metadata:
            update_ds.SetMetadata(rpc_metadata, "RPC")
        elif georef is not None and georef.get("rpc_corrected"):
            update_ds.SetMetadata(
                {str(key): str(value)
                 for key, value in georef["rpc_corrected"].items()},
                "RPC"
            )
        if georef is not None:                           # 有步骤1元数据
            # 以下仅为步骤1 元数据留档，本步几何完全由栅格自身地理参考决定
            update_ds.SetMetadataItem("GEOREF_PARAMS", IN_GEOREF_NAME)  # 元数据文件名
            update_ds.SetMetadataItem("GEOREF_TARGET_CRS",  # 步骤1目标CRS
                                      "EPSG:%s" % georef.get("target_crs_epsg", "?"))  # EPSG文本
            update_ds.SetMetadataItem("GEOREF_TARGET_RES",  # 步骤1目标分辨率
                                      "%g" % float(georef.get("target_res", 0.0)))  # 分辨率文本
        update_ds = None                                 # 关闭并落盘元数据

        # ⑥ 续：原子替换主成果并写出参数/统计
        np.save(tmp_coef, coefficients)                  # 系数写入临时npy
        np.save(tmp_wl, wavelengths)                     # 波长写入临时npy
        os.replace(tmp_tif, final_tif)                   # 原子替换反射率成果
        written_targets.append(final_tif)                # 记录已写出文件
        os.replace(tmp_coef, final_coef)                 # 原子替换系数文件
        written_targets.append(final_coef)               # 记录已写出文件
        os.replace(tmp_wl, final_wl)                     # 原子替换波长文件
        written_targets.append(final_wl)                 # 记录已写出文件
        georef_local = None                              # 本地元数据路径初值
        if georef is not None and os.path.isfile(in_georef):  # 有元数据可复制
            # v2.2.0：重写波段口径，与最终输出(剔除后)一致，避免下游按 georef 误读 166 波段
            try:
                with open(in_georef, "r", encoding="utf-8") as _fh:
                    georef_out = json.load(_fh)
            except Exception:
                georef_out = dict(georef)
            georef_out["band_count"] = int(len(wavelengths))
            georef_out["source_band_count"] = int(original_band_count)
            georef_out["kept_bands_1based"] = [int(i) + 1 for i in range(original_band_count) if i not in drop_0b]
            georef_out["kept_wavelengths_um"] = [float(w) for w in wavelengths]
            georef_out["removed_bands_1based"] = [int(b) for b in drop_1based]
            georef_out["removed_bands_0based"] = sorted(int(b) for b in drop_0b)
            georef_out["removed_wavelengths_um"] = [float(original_wavelengths_full[i]) for i in sorted(drop_0b)]
            atomic_write_json(georef_out, final_georef)  # 写重写后的元数据
            written_targets.append(final_georef)         # 记录已写出文件
            georef_local = final_georef                  # 记录本地元数据路径
        write_band_stats_csv(stats, final_stats)         # 写出波段统计CSV

        # v2.1.0：NaN 修复台账（微调平移波段 + 邻带一致性校验结论），可复核追溯
        repair_payload = {
            "scene": scene,
            "strategy": "fine_tune_pm1nm -> neighbor_interpolate -> raise",
            "consistency_thresholds": {
                "xa_rel_max": COEF_CONSISTENCY_MAX_REL_XA,
                "xb_xc_rel_max": COEF_CONSISTENCY_MAX_REL_XBXC,
            },
            "perturbed_bands": [
                {
                    "band_1based": int(band),
                    "original_wavelength_um": round(float(original_wavelengths_full[band - 1]), 6)
                    if band - 1 < len(original_wavelengths_full) else None,
                    "wavelength_used_um": round(float(value[0]), 6),
                    "offset_nm": round(float(value[1]) * 1000.0, 3),
                    "method": "fine_tune_then_neighbor_consistency_check",
                }
                for band, value in sorted(perturbed.items())
            ],
            "note": "如 perturbed_bands 为空，说明本景 6S 系数无 NaN 坏点，未触发修复。",
        }
        atomic_write_json(repair_payload, final_repair_log)  # 原子写出修复台账
        written_targets.append(final_repair_log)         # 纳入失败清理清单

        payload = {                                      # 参数记录字典
            "scene": scene,                              # 场景名
            "source_radiance": in_radiance,              # 输入辐亮度路径
            "source_tar": meta["tar"],                   # 原始压缩包
            "georef_json": in_georef if georef is not None else None,  # 步骤1元数据
            "georef_json_local": georef_local,           # 输出目录元数据
            "band_count": int(len(wavelengths)),         # 剔除后输出波段数
            "source_band_count": int(original_band_count),  # 剔除前/6S参与计算波段数
            "band_detection": {                          # 反射率域波段检测与剔除记录
                "dropped_1based": drop_1based,           # 被剔除波段（1-based）
                "dropped_count": len(drop_1based),       # 剔除波段数
                "dropped_water_vapor_1based": water_vapor_dropped_1b,  # 水汽带（v2.1.2）
                "dropped_by_statistics_1based": stat_dropped_1b,       # 统计判定剔除
                "water_vapor_bands_um": WATER_VAPOR_BANDS_UM,          # 水汽带定义
                "kept_count": int(len(wavelengths)),     # 保留波段数
                "report": OUT_BAND_DETECT_REPORT,        # 检测报告文件名
            },
            "raster_size": [int(source_ds.RasterXSize),  # 影像列数
                            int(source_ds.RasterYSize)],  # 影像行数
            "date": meta["date"].strftime("%Y-%m-%d %H:%M:%S")  # 成像时间
            if meta.get("date") else None,               # 无日期时为None
            "satellite": meta.get("satellite"),          # 卫星名
            "sensor": meta.get("sensor"),                # 传感器名
            "solar_zenith": meta["solar_z"],             # 太阳天顶角
            "solar_azimuth": meta["solar_a"],            # 太阳方位角
            "view_zenith": meta["view_z"],               # 观测天顶角
            "view_azimuth": meta["view_a"],              # 观测方位角
            "center_lat": meta.get("center_lat"),        # 中心纬度
            "center_lon": meta.get("center_lon"),        # 中心经度
            "footprint_wgs84": {key: [float(value[0]), float(value[1])]  # 足迹经纬度
                                for key, value in meta["footprint"].items()},  # 遍历四角
            "atmos_profile": str(atmos_profile),         # 大气模式
            "atmos_source": atmos_source,                # 大气模式来源
            "aero_profile": str(aero_profile),           # 气溶胶类型
            "aero_source": aero_source,                  # 气溶胶来源
            "aot550": aot550,                            # AOT550
            "aot_source": aot_detail.get("source", ""),   # AOT来源
            "aot_detail": aot_detail,                    # AOT明细
            "elevation_km": elevation_km,                # 海拔
            "elevation_bands_enabled": bool(elev_band_map is not None),  # 是否启用高程分块
            "elevation_bands_km": [float(v) for v in band_elevs_km] if (elev_band_map is not None and band_elevs_km) else None,  # 各档海拔
            "elevation_source": elevation_detail.get("source", ""),  # 海拔来源
            "elevation_detail": elevation_detail,        # 海拔明细
            "water_vapor": ("water=%.3f g/cm2, ozone=%.3f cm-atm"  # 水汽臭氧文本
                            % (float(WATER_VAPOR_CM),    # 水汽含量
                               float(OZONE_CM if OZONE_CM is not None  # 臭氧取值
                                     else OZONE_CM_DEFAULT)))  # 默认臭氧
            if WATER_VAPOR_CM is not None else None,     # 未配置水汽则为None
            "sixs_path": sixs_exe,                       # 6S程序路径
            "sixs_version": sixs_version,                # 6S版本
            "py6s_version": resolve_py6s_version(),      # Py6S版本
            "sixs_bands": int(original_band_count),      # 参与6S计算的波段数（剔除前全量）
            "sixs_elapsed_s": round(float(coef_elapsed), 2),  # 系数计算耗时
            "sixs_perturbed_bands": [                    # 波长平移重试记录
                {"band": int(band),                      # 波段号
                 "wavelength_used_um": round(float(value[0]), 6),  # 实际波长
                 "offset_nm": round(float(value[1]) * 1000.0, 3)}  # 偏移量(nm)
                for band, value in sorted(perturbed.items())  # 按波段排序
            ],
            "formula": "rho = (xa*L - xb) / (1 + xc*(xa*L - xb))",  # 校正公式
            "reflectance_clip": list(REFLECTANCE_CLIP) if REFLECTANCE_CLIP else None,  # 裁剪区间
            "nodata": NODATA_OUT,                        # 无效值
            # 输出网格直接来自正射辐亮度栅格；georef.json 仅作元数据参考
            "output_grid_crs_epsg": source_epsg,         # 输出网格EPSG（native 网格为 None）
            "output_grid_res": source_res,               # 输出分辨率（native 网格为 None）
            "grid_kind": "native" if is_native_grid else "ortho",  # 网格类型
            "target_crs_epsg": (georef.get("target_crs_epsg")  # 步骤1目标EPSG
                                if georef else source_epsg),  # 无元数据用栅格值
            "target_res": (georef.get("target_res") if georef else source_res),  # 目标分辨率
            "output_grid_size": [int(source_ds.RasterXSize),  # 输出列数
                                 int(source_ds.RasterYSize)],  # 输出行数
            "geometry": "native grid; RPC orthorectification deferred to step 3",  # 几何说明
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),  # 生成时间
        }
        write_parameters(payload, final_params_json, final_params_txt, logger)  # 写参数记录
        main_complete = True                             # 标记主成果已完成
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        # ⑧ 生成带坐标预览和质检图
        current_stage = 8                                # 进入第8阶段
        stage_description = "生成降采样预览和质检图"       # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 8,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        preview_tif = None                               # 预览路径初值
        if MAKE_QC and not args.no_qc:                   # 启用质检
            if QC_GEOREF_PREVIEW:                        # 启用带坐标预览
                try:                                     # 预览生成可能失败
                    prepare_qc_preview(final_tif, wavelengths, final_preview,  # 生成预览
                                       logger)           # 传入日志器
                    preview_tif = final_preview          # 记录预览路径
                except Exception:                        # 预览失败降级
                    logger.warning("带坐标预览生成失败，质检图回退无坐标降采样",  # 告警
                                   exc_info=True)        # 附带堆栈
            # v2.2.0：质检图为用户必需成果，生成失败即停（不再 warning 降级）
            qc_ok = write_quality_figure(final_tif, wavelengths, final_qc, scene,  # 生成质检图
                                         params_summary_text(payload), logger,  # 参数摘要
                                         band_stats=stats, preview_tif=preview_tif)  # 统计与预览
            if not qc_ok:
                raise RuntimeError("质检图生成失败（用户必需成果，不允许静默跳过）: %s" % final_qc)
        else:                                            # 已关闭质检
            logger.info("按配置跳过质检图生成")           # 记录跳过
        log_stage_done(logger, current_stage, 8, stage_description, stage_started)  # 阶段结束

        elapsed = time.time() - started                  # 本景总耗时
        valid_ratio = 100.0 * float(np.mean([item["valid"] for item in stats])) \
            / max(int(source_ds.RasterXSize) * int(source_ds.RasterYSize), 1)  # 除以总像元数
        logger.info("完成: %d 波段，%d×%d，平均有效像元 %.1f%%，"  # 完成信息
                    "6S %s，总计 %.1f s，大小 %s",            # 后续字段
                    len(wavelengths), source_ds.RasterXSize,  # 波段数与列
                    source_ds.RasterYSize, valid_ratio,   # 行数与有效率
                    fmt_elapsed(coef_elapsed), elapsed,   # 系数耗时与总耗时
                    human_size(os.path.getsize(final_tif)))  # 成果文件大小
        logger.info("  成果目录: %s", out_dir)             # 记录成果目录
        return {"scene": scene, "status": "成功", "bands": int(len(wavelengths)),  # 返回成功
                "elapsed": elapsed,                       # 耗时
                "message": "%d 波段，AOT=%.3f，海拔=%.3f km，%.1f s"  # 摘要文本
                           % (len(wavelengths), aot550, elevation_km, elapsed)}  # 摘要参数
    except Exception:                                    # 本景处理失败
        for path in (tmp_tif, tmp_subset, tmp_coef, tmp_wl):  # 清理临时文件（含波段剔除中间产品）
            try:                                         # 删除可能失败
                if os.path.exists(path):                 # 文件存在
                    os.remove(path)                      # 删除临时文件
            except OSError:                              # 忽略删除错误
                pass                                     # 继续清理
        # 本次未完整写出参数记录时，清掉已替换的半成品，避免下次误判为已完成
        if not main_complete:                            # 主成果未完成
            for path in written_targets:                 # 遍历已写出文件
                try:                                     # 删除可能失败
                    if os.path.exists(path):             # 文件存在
                        os.remove(path)                  # 删除半成品
                except OSError:                          # 忽略删除错误
                    pass                                 # 继续清理
        if created_dir:                                  # 目录是本景新建的
            try:                                         # 删除可能失败
                if not os.listdir(out_dir):              # 目录已空
                    os.rmdir(out_dir)                    # 删除空目录
            except OSError:                              # 忽略删除错误
                pass                                     # 继续向上抛
        logger.error("场景处理在步骤 %d/8 失败：%s", current_stage,  # 记录失败阶段
                     stage_description)                  # 阶段说明
        raise                                            # 向上抛出原异常
    finally:                                             # 无论如何都释放输入
        source_ds = None                                 # 关闭输入数据集


# =================================================================
# 批处理主流程
# =================================================================

def main(argv=None):
    args = parse_args(argv)                              # 解析命令行参数
    out_root = os.path.abspath(args.out_root)            # 输出根目录绝对路径
    os.makedirs(out_root, exist_ok=True)                 # 确保输出根目录存在
    log_path = args.log or os.path.join(out_root, "step2_batch.log")  # 日志路径
    logger = setup_logging(log_path, args.verbose)       # 初始化日志
    gdal.SetCacheMax(GDAL_CACHE_MB * 1024 * 1024)        # 设置GDAL缓存上限

    in_root = os.path.abspath(args.in_root)              # 输入根目录绝对路径
    logger.info("开始扫描步骤1成果目录：%s", in_root)     # 记录扫描目录
    try:                                                 # 扫描失败可返回错误码
        scene_dirs = discover_scene_dirs(in_root)        # 扫描步骤1成果
    except Exception as exc:                             # 目录不存在等异常
        logger.error("扫描步骤1 成果失败: %s", exc)       # 记录错误
        return 2                                         # 返回错误码

    if args.scene:                                       # 指定了场景过滤
        patterns = list(args.scene)                      # 复制过滤模式
        scene_dirs = [path for path in scene_dirs        # 逐个匹配保留
                      if any(fnmatch.fnmatch(os.path.basename(path), pattern)  # 按通配匹配
                             for pattern in patterns)]   # 遍历模式
    if args.limit > 0:                                   # 限制处理数量
        scene_dirs = scene_dirs[: args.limit]            # 截取前N景
    if not scene_dirs:                                   # 没有待处理场景
        if args.scene:                                   # 是过滤导致为空
            logger.error("在 %s 下没有匹配 --scene %s 的步骤1 成果（%s）",  # 记录过滤
                         in_root, ", ".join(args.scene),  # 目录与模式
                         "*" + IN_RADIANCE_SUFFIX)       # 需要的文件名后缀
        else:                                            # 目录本身无成果
            logger.error("在 %s 下未找到可处理的步骤1 成果（%s）",  # 记录原因
                         in_root, "*" + IN_RADIANCE_SUFFIX)  # 目录与后缀
        return 2                                         # 返回错误码

    logger.info("=" * 72)                                # 打印分隔线
    logger.info(                                         # 打印程序标题
        "ZY1E AHSI 步骤2：6S 大气校正（正射格网，输入 *%s）",  # 标题文本
        IN_RADIANCE_SUFFIX)                              # 输入后缀
    logger.info("步骤1 成果: %s", in_root)                # 输入目录
    logger.info("输出根目录: %s", out_root)               # 输出目录
    logger.info("待处理    : %d 景", len(scene_dirs))     # 待处理景数
    logger.info("辅助数据  : DEM=%s", DEM_DIR)            # 辅助数据目录
    logger.info("处理选项  : AOT方法=%s, 覆盖已有成果=%s, 质检图=%s",  # 处理选项
                args.aot_method or AOT_METHOD, args.overwrite,  # AOT方法与覆盖
                MAKE_QC and not args.no_qc)              # 是否出质检图
    logger.info("日志文件  : %s", log_path)              # 日志路径
    logger.info("=" * 72)                                # 打印分隔线

    if args.dry_run:                                     # 仅演练不执行
        for index, scene_dir in enumerate(scene_dirs, start=1):  # 逐景打印计划
            scene = os.path.basename(os.path.normpath(scene_dir))  # 取场景名
            logger.info("[%d/%d] %s", index, len(scene_dirs), scene)  # 打印序号
            logger.info("        <- %s", resolve_radiance_path(scene_dir))  # 输入文件
            logger.info("        -> %s", os.path.join(           # 输出文件
                out_root, scene, scene + OUT_REFLECTANCE_SUFFIX))  # 拼输出路径
        logger.info("dry-run 结束，未执行任何处理")        # 提示演练结束
        return 0                                         # 正常退出

    results = []                                         # 逐景处理结果
    started = time.time()                                # 批处理起始时刻
    for index, scene_dir in enumerate(scene_dirs, start=1):  # 逐景处理
        scene = os.path.basename(os.path.normpath(scene_dir))  # 取场景名
        logger.info("-" * 72)                            # 打印分隔线
        logger.info("开始处理场景 (%d/%d)：%s", index, len(scene_dirs), scene)  # 记录开始
        scene_started = time.time()                      # 本景起始时刻
        try:                                             # 单景失败不中断批处理
            result = process_scene(scene_dir, out_root, args, logger)  # 处理一景
            results.append(result)                       # 收集结果
            logger.info("场景结束 (%d/%d)：%s，状态=%s，用时 %s",  # 记录结束
                        index, len(scene_dirs), scene, result["status"],  # 序号与状态
                        fmt_elapsed(time.time() - scene_started))  # 打印耗时
        except Exception as exc:                         # 单景异常
            logger.exception("处理失败: %s", scene)       # 记录堆栈
            logger.error("场景结束 (%d/%d)：%s，状态=失败，用时 %s",  # 记录失败
                         index, len(scene_dirs), scene,  # 序号与场景
                         fmt_elapsed(time.time() - scene_started))  # 打印耗时
            results.append({"scene": scene, "status": "失败", "bands": None,  # 失败结果
                            "elapsed": 0.0, "message": str(exc)})  # 异常信息
            if not args.continue_on_error:               # v2.1.3：默认出错即停
                logger.error("出错即停（如需继续请加 --continue-on-error）: %s", scene)
                raise                                    # 非零退出，不留半截成果

    elapsed_total = time.time() - started                # 批处理总耗时
    success = [item for item in results if item["status"] == "成功"]  # 成功景
    skipped = [item for item in results if item["status"].startswith("跳过")]  # 跳过景
    failed = [item for item in results if item["status"] == "失败"]  # 失败景

    logger.info("=" * 72)                                # 打印分隔线
    logger.info("批处理结束: 成功 %d，跳过 %d，失败 %d，总用时 %s",  # 汇总统计
                len(success), len(skipped), len(failed),  # 各状态数量
                fmt_elapsed(elapsed_total))              # 总耗时
    for result in results:                               # 逐景输出结果
        note = result.get("message") or ""               # 取附加说明
        logger.info("  [%s] %s%s", result["status"], result["scene"],  # 打印状态与场景
                    (" - " + note) if note else "")      # 有说明才拼接

    if failed:                                           # 有失败场景
        logger.error("有 %d 景处理失败，详见日志: %s", len(failed), log_path)  # 记录错误
        return 1                                         # 返回失败码
    logger.info("全部完成，成果位于: %s", out_root)        # 提示全部完成
    return 0                                             # 返回成功码


if __name__ == "__main__":                               # 脚本直接运行时
    sys.exit(main())                                     # 用返回码退出
