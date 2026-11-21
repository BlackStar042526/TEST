# -*- coding: utf-8 -*-
"""
02_atmospheric_correction.py
ZY1E AHSI Step 2：K-T 气溶胶反演 + 6S/Py6S 大气校正（v2.0.0）

v2.0.0 重构要点（GitHub 比对 + 用户需求）：
  ✅ 废弃 MODIS/AERONET 外部气溶胶产品 → K-T 暗像元 2-Band 反演 (SWIR 2100-2250 nm + Red 640-680 nm)
  ✅ K-T 失败回退：Koschmieder 能见度法（默认 60 km 晴朗 → AOT550 ≈ 0.065）
  ✅ 6S 轨道高度：set_sensor_custom_altitude(778.0)（ZY1E 实际轨道高度，v1.1.4 错用默认 100km）
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

# =================================================================
# 配置区（按需修改）
# =================================================================

# 输入步骤1 正射成果目录；输出大气校正成果根目录
IN_ROOT_STEP1 = r"E:\ZhouYang\test\ZY1E_output_step1"   # 步骤1正射成果输入根目录
OUT_ROOT_STEP2 = r"E:\ZhouYang\test\ZY1E_output_step2"  # 本步大气校正成果输出根目录

# 原始 TAR 目录；读取观测几何和日期，优先使用 georef.json 的 source_tar
TAR_DIR = r"E:\ZhouYang\test\ZY1E_AHSI"                 # 原始TAR目录，用于读观测几何

# MODIS MOD04/MYD04 AOD 产品目录（按日期及场景足迹选有效像元）

# ASTGTM DEM 目录；可放置 GeoTIFF/IMG 或包含栅格的 ZIP
DEM_DIR = r"E:\ZhouYang\test\DEM"                        # DEM目录，用于算海拔

# 6S 可执行程序路径；留空时依次查环境变量、脚本 _bin 目录和 PATH
SIXS_EXE = ""                                            # 6S程序路径，空则自动查找

# 自动按日期和纬度选择大气模式；关闭时使用手工模式
AUTO_DETECT_ATMOS = True                     # 按日期和纬度自动选大气模式
ATMOS_PROFILE_MANUAL = "MidlatitudeWinter"  # Tropical/MidlatitudeSummer/MidlatitudeWinter/SubarcticSummer/SubarcticWinter/USStandard

# 气溶胶类型（无土地覆盖数据时默认 Continental）
AERO_PROFILE_MANUAL = "Continental"  # Continental/Urban/Maritime/Desert
AUTO_DETECT_AERO = False             # 启用时按经纬度粗略判断

# 暗像元反射比 AOT 模式：None / KT_DARK_PIXEL / 2-BAND_OVER_WATER
AOT_METHOD = "KT_DARK_PIXEL"                 # 直接指定暗像元反射比算法
AOT_METHOD_CHOICES = ("None", "2-Band (K-T)", "2-Band Over Water")
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

    # 手工覆盖 6S 关键参数（便于快速试验，不写死文件）
    parser.add_argument("--aot", type=float, default=None,         # 手工指定AOT550
                        help="直接指定 AOT550（覆盖自动 AOD 检索）")  # 参数帮助文本
    parser.add_argument("--aot-method", choices=AOT_METHOD_CHOICES,  # 暗像元反射比算法
                        default=None,
                        help="AOT 算法：None / 2-Band (K-T) / 2-Band Over Water")
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
    K-T 暗像元反演 AOT550（整景均匀值）。
    
    流程：
      1. 找 SWIR (2100-2250 nm, 中心 2175 μm) 和 Red (640-680 nm, 中心 660 μm) 波段
      2. 取 SWIR 辐亮度最低的 0.1% 像元作为暗像元
      3. 用 6S 正向模拟建立 AOT ∈ {0, 0.05, ..., 0.6} → 红光辐亮度查找表
      4. 在查找表里匹配暗像元红光辐亮度中位数
      5. 回退：Koschmieder 能见度法
    
    参数：
      wavelengths_um: 保留波段的中心波长数组（μm）
      radiance_cube: (n_bands, H, W) Float32 native 网格辐亮度
                     单位 W·m⁻²·sr⁻¹·μm⁻¹（v2.0.0 确认无 ×10 换算）
      meta: XML 观测几何元数据
      sixs_builder_func: 构造 Py6S 实例的函数（传入 atmos/aero/aot 参数）
    
    返回 (aot550, info_dict)
    """
    # --- Step 1: 定位 K-T 波段 ---
    swir_idx = _nearest_band(wavelengths_um, 2.175)  # 2175 μm
    red_idx = _nearest_band(wavelengths_um, 0.660)   # 660 μm
    swir_wl = wavelengths_um[swir_idx]
    red_wl = wavelengths_um[red_idx]

    kt_ok = (2.10 <= swir_wl <= 2.25) and (0.64 <= red_wl <= 0.68)
    logger.info("K-T 波段: SWIR idx=%d (%.4f μm, FLAASH 2.10-2.25 μm ✓), "
                "Red idx=%d (%.4f μm, FLAASH 0.64-0.68 μm %s)",
                swir_idx + 1, swir_wl, red_idx + 1, red_wl,
                "✓" if kt_ok else "⚠ 偏离标准范围")

    if not kt_ok:
        logger.warning("K-T 波段不在 FLAASH 标准范围内，反演可能不准")

    # --- Step 2: 取有效辐亮度 ---
    swir_l = radiance_cube[swir_idx].astype(np.float64)
    red_l = radiance_cube[red_idx].astype(np.float64)
    valid_mask = (np.isfinite(swir_l) & np.isfinite(red_l)
                  & (swir_l > 0.0) & (red_l > 0.0))
    n_valid = int(valid_mask.sum())
    if n_valid < 1000:
        logger.warning("有效像元不足（%d），K-T 反演失效", n_valid)
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "visibility_km": VISIBILITY_MANUAL_KM,
            "reason": "insufficient_valid_pixels",
        }

    # --- Step 3: 6S 正向查找表 ---
    aot_candidates = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5, 0.6]
    red_l_simulated = []
    for aot_hypo in aot_candidates:
        try:
            s = sixs_builder_func(aot=aot_hypo, wavelength=float(red_wl))
            s.run()
            red_l_simulated.append(float(s.outputs.pixel_radiance))
        except Exception as e:
            logger.warning("6S 正向模拟失败 AOT=%.2f: %s", aot_hypo, e)
            red_l_simulated.append(np.nan)

    # 跳过 NaN
    valid_lut = [(a, l) for a, l in zip(aot_candidates, red_l_simulated) if np.isfinite(l)]
    if len(valid_lut) < 3:
        logger.warning("6S 查找表有效点不足，回退能见度法")
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "reason": "lut_too_few_points",
        }

    lut_aot = np.array([v[0] for v in valid_lut])
    lut_red = np.array([v[1] for v in valid_lut])

    # --- Step 4: 暗像元选择 + 匹配 ---
    swir_valid = swir_l[valid_mask]
    percentile_0p1 = np.percentile(swir_valid, 0.1)
    dark_mask = valid_mask & (swir_l <= percentile_0p1)
    n_dark = int(dark_mask.sum())

    if n_dark < 10:
        logger.warning("暗像元不足（%d 个），K-T 反演失效", n_dark)
        return visibility_to_aot550(VISIBILITY_MANUAL_KM), {
            "source": "visibility_fallback", "visibility_km": VISIBILITY_MANUAL_KM,
            "reason": "no_dark_pixels",
        }

    red_measured_dark = float(np.median(red_l[dark_mask]))
    aot_estimated = float(np.clip(np.interp(red_measured_dark, lut_red, lut_aot), 0.01, 1.0))

    logger.info("K-T 反演: AOT550=%.4f, 暗像元 %d 个, SWIR p0.1=%.3f, Red median=%.3f",
                aot_estimated, n_dark, percentile_0p1, red_measured_dark)

    return aot_estimated, {
        "source": "Kaufman-Tanré 2-Band",
        "method": "6S forward LUT + SWIR 0.1% percentile dark pixels",
        "aot550": aot_estimated,
        "visibility_km": float(3.912 / aot_estimated),
        "swir_band_1based": int(swir_idx + 1),
        "swir_wavelength_um": float(swir_wl),
        "red_band_1based": int(red_idx + 1),
        "red_wavelength_um": float(red_wl),
        "dark_pixel_count": n_dark,
        "dark_pixel_swir_p01": float(percentile_0p1),
        "dark_pixel_red_median": red_measured_dark,
        "lut_aot": [round(a, 2) for a in lut_aot.tolist()],
        "lut_red": [round(l, 4) for l in lut_red.tolist()],
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

    if method_key == "KT_DARK_PIXEL":
        try:
            return resolve_aot550_kt(
                wavelengths_um, radiance_cube, meta, sixs_builder_func, logger
            )
        except Exception as exc:
            logger.warning("K-T 反演异常: %s，回退固定能见度", exc)

    if method_key == "KT_OVER_WATER":
        try:
            return resolve_aot550_over_water(
                wavelengths_um, radiance_cube, meta, sixs_builder_func, logger
            )
        except Exception as exc:
            logger.warning("2-Band Over Water 反演异常: %s，回退固定能见度", exc)

    aot = visibility_to_aot550(VISIBILITY_MANUAL_KM)
    logger.info("固定能见度回退 AOT550=%.4f (visibility=%d km)",
                aot, VISIBILITY_MANUAL_KM)
    return aot, {
        "source": "visibility_fallback",
        "method": method,
        "visibility_km": VISIBILITY_MANUAL_KM,
        "aot550": aot,
    }



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
        return None                                      # 返回无结果
    logger.info("DEM 数据: %d 个瓦片", len(paths))        # 打印瓦片数量

    vrt_path = "/vsimem/zy1e_step2_dem.vrt"              # 内存虚拟VRT路径
    vrt = None                                           # VRT句柄初值
    try:                                                 # 构建与裁剪可能失败
        vrt = gdal.BuildVRT(vrt_path, paths)             # 把多个瓦片拼成VRT
        if vrt is None:                                  # 构建失败
            raise RuntimeError("DEM VRT 构建失败")        # 报错
        first = gdal.Open(paths[0])                      # 打开首个瓦片
        if first is None:                                # 打开失败
            raise RuntimeError("无法打开 DEM 瓦片: %s" % paths[0])  # 报错
        srs_wkt = first.GetProjection()                  # 读取DEM投影
        geo_transform = first.GetGeoTransform()          # 读取DEM仿射变换
        first = None                                     # 关闭首个瓦片
        if not srs_wkt:                                  # DEM没有投影
            logger.warning("DEM 缺少投影信息，跳过自动海拔")  # 记录告警
            return None                                  # 返回无结果

        srs_dem = osr.SpatialReference()                 # DEM坐标系对象
        srs_dem.ImportFromWkt(srs_wkt)                   # 由WKT导入
        srs_wgs84 = osr.SpatialReference()               # WGS84坐标系对象
        srs_wgs84.ImportFromEPSG(4326)                   # 设为地理坐标系
        srs_dem.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)  # 传统经纬顺序
        srs_wgs84.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)  # 传统经纬顺序
        transformer = osr.CoordinateTransformation(srs_wgs84, srs_dem)  # 建坐标转换
        points = [transformer.TransformPoint(lon, lat) for lon, lat in ring]  # 足迹转DEM坐标
        min_x = min(point[0] for point in points)         # 西边界
        max_x = max(point[0] for point in points)         # 东边界
        min_y = min(point[1] for point in points)         # 南边界
        max_y = max(point[1] for point in points)         # 北边界

        dataset = gdal.Warp(                             # 内存中裁剪重投影
            "", vrt, format="MEM",                       # 输出到内存
            dstSRS=srs_dem.ExportToWkt(),                # 目标坐标系
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
        logger.info("DEM 足迹矩形网格读取完成：%d×%d @ %g m",  # 打印网格信息
                    dataset.RasterXSize, dataset.RasterYSize, DEM_RES)  # 行列数与分辨率
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
        return None                                      # 返回无结果
    mean_m = float(elevation[valid].mean())              # 有效像元平均高程（m）
    coverage = 100.0 * count / max(total, 1)             # 有效覆盖率（%）
    detail = {                                           # 海拔统计明细
        "source": "ASTGTM DEM (%d tiles)" % len(paths),  # 数据来源与瓦片数
        "tiles": paths,                                  # 使用的瓦片路径
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
    return detail                                        # 返回海拔明细


def resolve_elevation(meta, ring, args, logger):
    """确定目标海拔（km）：命令行 → DEM 平均 → 手工兜底，并记录来源。"""
    """确定目标海拔（km）：命令行 > DEM > 手工兜底。"""
    if args.elev is not None:                            # 命令行直接给了海拔
        logger.info("海拔 = %.3f km（命令行指定）", args.elev)  # 记录来源
        return float(args.elev), {"source": "命令行 --elev",  # 返回海拔与来源
                                  "elevation_km": float(args.elev)}  # 海拔值
    if AUTO_DETECT_DEM:                                  # 启用了DEM自动计算
        detail = compute_mean_elevation(ring, DEM_DIR, logger)  # 用DEM算平均海拔
        if detail is not None:                           # 计算成功
            return detail["elevation_km"], detail        # 返回海拔与明细
    logger.warning("DEM 自动计算失败，使用兜底海拔 %.3f km", ELEVATION_KM_FALLBACK)  # 告警
    return float(ELEVATION_KM_FALLBACK), {"source": "ELEVATION_KM_FALLBACK",  # 返回兜底海拔
                                          "elevation_km": float(ELEVATION_KM_FALLBACK)}  # 海拔值


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
    sixs.altitudes.set_sensor_custom_altitude(SENSOR_ALTITUDE_KM)  # v2.0.0: ZY1E 实际 778 km（旧版错用默认 100 km）          # 传感器按卫星高度
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


def write_reflectance(src_ds, coefficients, wavelengths, tmp_path, logger):
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

                y_value = coef[0] * radiance - coef[1]
                denominator = 1.0 + coef[2] * y_value
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


def _nearest_band(wavelengths, target):
    return int(np.argmin(np.abs(wavelengths - target)))  # 取最接近目标波长的下标


def _stretch(array, valid):
    values = array[valid]                                # 取出有效像元
    if values.size == 0:                                 # 没有有效像元
        return np.zeros_like(array)                      # 返回全零拉伸图
    low, high = np.percentile(values, [2, 98])           # 取2%和98%分位
    if high <= low:                                      # 分位异常退化
        high = low + 1e-6                                # 避免除零
    return np.clip((array - low) / (high - low), 0.0, 1.0)  # 线性拉伸到0~1


def prepare_qc_preview(reflectance_path, wavelengths, preview_tif, logger):
    """生成带地理参考的降采样 RGB 预览 tif，供质检图叠加坐标。"""
    """对正射反射率的 5 个预览波段直接降采样出带坐标预览，不做 RPC/warp。

    反射率已与步骤1 正射辐亮度同网格（EPSG:32650/30 m），这里只用
    gdal.Translate 按 QC_PREVIEW_SCALE 降采样（默认 60 m），保留原 CRS。
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
        logger.info("  带坐标预览: %s (%d 波段 @ %g m，降采样自正射反射率)",  # 记录预览信息
                    preview_tif, len(unique), preview_res)  # 路径、波段数与分辨率
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
            srs = osr.SpatialReference()                 # 新建空间参考
            srs.ImportFromWkt(dataset.GetProjection())   # 导入投影
            epsg = srs.GetAuthorityCode(None) or "?"     # 取EPSG编号
            extent = [geo_transform[0] / 1000.0,         # 西边界（km）
                      (geo_transform[0] + geo_transform[1] * width) / 1000.0,  # 东边界
                      (geo_transform[3] + geo_transform[5] * height) / 1000.0,  # 南边界
                      geo_transform[3] / 1000.0]          # 北边界
            dataset = None                               # 关闭数据集
            crs_text = "EPSG:%s @ %g m preview" % (epsg, abs(geo_transform[1]))  # 坐标系说明
            return (arrays, extent, "Easting (km)", "Northing (km)",  # 坐标轴与范围
                    crs_text, "%d×%d" % (width, height), True)  # 说明文本与带坐标标志

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
    except Exception:                                    # 整体失败
        logger.warning("补齐质检图失败: %s", scene_dir, exc_info=True)  # 告警


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

    legacy_tif = next((path for path in legacy_tifs if os.path.isfile(path)),  # 取旧成果
                      None)                              # 没有则为None
    if (not args.overwrite and not os.path.isfile(final_tif)  # 未要求覆盖且无新成果
            and legacy_tif and os.path.isfile(final_coef)  # 有旧成果与系数
            and os.path.isfile(final_params_json)):      # 参数也齐全
        os.replace(legacy_tif, final_tif)                # 旧成果改名复用
        logger.info("复用已有完整成果并更新反射率文件名：%s", final_tif)  # 记录复用

    have_main = all(os.path.isfile(path)                 # 主成果是否齐全
                    for path in (final_tif, final_coef, final_params_json))  # 三项检查
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
    tmp_coef = os.path.join(out_dir, "py6s_coefficients.tmp.npy")  # 临时系数
    tmp_wl = os.path.join(out_dir, "wavelengths_clean.tmp.npy")  # 临时波长
    written_targets = []                                 # 已写出的目标文件
    main_complete = False                                # 主成果是否完成
    source_ds = None                                     # 输入数据集句柄
    try:                                                 # 单景异常统一收尾
        # ① 校验步骤1 输入（辐亮度/波长/波段数），输出沿用同一正射网格、不重采样
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
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
        if georef and georef.get("raster_size"):         # 元数据记录了尺寸
            expected = [int(georef["raster_size"][0]), int(georef["raster_size"][1])]  # 记录尺寸
            actual = [source_ds.RasterXSize, source_ds.RasterYSize]  # 实际尺寸
            if expected != actual:                       # 尺寸不一致
                logger.warning("正射辐亮度尺寸 %s 与 georef.json 记录 %s 不一致",  # 告警
                               actual, expected)         # 实际与记录

        # 输出网格完全继承输入正射栅格：读取 CRS/分辨率并做一致性校验（不重采样）
        source_gt = source_ds.GetGeoTransform()          # 读取仿射变换
        source_srs = osr.SpatialReference()              # 新建空间参考
        source_srs.ImportFromWkt(source_ds.GetProjection())  # 导入投影
        source_epsg = source_srs.GetAuthorityCode(None)  # 取EPSG编号
        source_res = abs(source_gt[1]) if source_gt and source_gt[1] else None  # 像元大小
        if source_epsg is None or source_res is None:    # 缺CRS或分辨率
            logger.warning("正射辐亮度缺少有效 CRS/仿射变换，输出将沿用该空地理参考")  # 告警
        else:                                            # 地理参考有效
            logger.info("正射网格: EPSG:%s @ %g m, %d×%d（不做重采样）",  # 打印网格
                        source_epsg, source_res,         # EPSG与分辨率
                        source_ds.RasterXSize, source_ds.RasterYSize)  # 行列数
            if georef and georef.get("target_crs_epsg") is not None:  # 元数据有EPSG
                try:                                     # 数值转换可能失败
                    if int(georef["target_crs_epsg"]) != int(source_epsg):  # EPSG不一致
                        logger.warning("正射辐亮度 CRS (EPSG:%s) 与 georef.json "  # 告警
                                       "target_crs_epsg=%s 不一致，以栅格为准",  # 说明
                                       source_epsg, georef["target_crs_epsg"])  # 两值
                except (TypeError, ValueError):          # 非法值
                    pass                                 # 忽略
            if georef and georef.get("target_res") is not None:  # 元数据有分辨率
                try:                                     # 数值转换可能失败
                    if abs(float(georef["target_res"]) - source_res) > 1e-6:  # 分辨率不一致
                        logger.warning("正射辐亮度分辨率 %g m 与 georef.json "  # 告警
                                       "target_res=%s 不一致，以栅格为准",  # 说明
                                       source_res, georef["target_res"])  # 两值
                except (TypeError, ValueError):          # 非法值
                    pass                                 # 忽略

        logger.info("正射辐亮度: %d 波段, %d×%d, 波长 %.4f~%.4f μm",  # 打印输入信息
                    source_ds.RasterCount, source_ds.RasterXSize,  # 波段与列
                    source_ds.RasterYSize, wavelengths.min(), wavelengths.max())  # 行与波长
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

        # ② 读取观测几何、日期、中心位置和足迹
        current_stage = 2                                # 进入第2阶段
        stage_description = "解析成像时间、观测几何和场景足迹"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        meta = resolve_metadata(scene, scene_dir, georef, TAR_DIR, logger)  # 汇总元数据
        ring = footprint_ring(meta["footprint"])         # 生成足迹闭合环
        logger.info("观测几何: 太阳 %.2f°/%.2f°, 卫星 %.2f°/%.2f°, 日期 %s",  # 打印几何
                    meta["solar_z"], meta["solar_a"], meta["view_z"],  # 太阳与观测角
                    meta["view_a"], meta["date"].strftime("%Y-%m-%d"))  # 方位角与日期
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

        # ③ 确定 AOT550、海拔、大气模式和气溶胶类型
        current_stage = 3                                # 进入第3阶段
        stage_description = "确定 AOT550、DEM 海拔及大气/气溶胶模式"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        elevation_km, elevation_detail = resolve_elevation(meta, ring, args, logger)  # 确定海拔
        atmos_profile, atmos_source = resolve_atmos_profile(meta, args, logger)  # 确定大气模式
        aero_profile, aero_source = resolve_aero_profile(meta, args, logger)  # 确定气溶胶
        aot_method = args.aot_method or AOT_METHOD
        if args.aot is not None or aot_method == "None":
            aot_wavelengths = np.asarray([], dtype=np.float64)
            radiance_cube = None
        else:
            target_wavelengths = (
                (2.175, 0.660) if aot_method == "2-Band (K-T)"
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
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

        # ④ 逐波段计算 6S 系数 (xa, xb, xc)
        current_stage = 4                                # 进入第4阶段
        stage_description = "运行 6S 并计算各波段校正系数"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        sixs = build_sixs(sixs_exe, meta, atmos_profile, aero_profile,  # 构造6S对象
                          aot550, elevation_km, logger)  # 传入AOT与海拔
        coefficients, sixs_version, coef_elapsed, perturbed = \
            compute_sixs_coefficients(sixs, wavelengths, logger)  # 逐波段跑6S取系数
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

        # ⑤ 逐波段反演反射率 rho=(xa*L-xb)/(1+xc*(xa*L-xb))，保持输入网格且不写 RPC
        current_stage = 5                                # 进入第5阶段
        stage_description = "逐波段计算并写出地表反射率"     # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        stats = write_reflectance(source_ds, coefficients, wavelengths,  # 写反射率
                                  tmp_tif, logger)       # 输出到临时文件
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

        # ⑥ 写入成果元数据、参数记录和波段统计
        current_stage = 6                                # 进入第6阶段
        stage_description = "写入影像元数据、参数记录和波段统计"  # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
                                        stage_description)  # 阶段说明
        update_ds = gdal.Open(tmp_tif, gdal.GA_Update)   # 以更新方式打开临时成果
        update_ds.SetMetadataItem("SOURCE_RADIANCE", os.path.basename(in_radiance))  # 输入文件名
        update_ds.SetMetadataItem("SOURCE_TAR",          # 原始压缩包名
                                  os.path.basename(meta["tar"])  # 有tar取文件名
                                  if meta["tar"] else "MISSING")  # 无tar标记缺失
        update_ds.SetMetadataItem("BAND_COUNT", str(source_ds.RasterCount))  # 波段数
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
        res_text = ("%g" % source_res) if source_res else "?"  # 分辨率文本
        update_ds.SetMetadataItem(                       # 写入输出网格
            "OUTPUT_GRID", "EPSG:%s @ %s m (inherited from %s)"  # 网格描述
            % (source_epsg or "?", res_text, IN_RADIANCE_SUFFIX))  # 参数
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
            copy_file_atomic(in_georef, final_georef)    # 复制元数据到输出目录
            written_targets.append(final_georef)         # 记录已写出文件
            georef_local = final_georef                  # 记录本地元数据路径
        write_band_stats_csv(stats, final_stats)         # 写出波段统计CSV

        payload = {                                      # 参数记录字典
            "scene": scene,                              # 场景名
            "source_radiance": in_radiance,              # 输入辐亮度路径
            "source_tar": meta["tar"],                   # 原始压缩包
            "georef_json": in_georef if georef is not None else None,  # 步骤1元数据
            "georef_json_local": georef_local,           # 输出目录元数据
            "band_count": int(source_ds.RasterCount),    # 波段数
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
            "sixs_bands": int(len(wavelengths)),         # 参与计算的波段数
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
            "output_grid_crs_epsg": source_epsg,         # 输出网格EPSG
            "output_grid_res": source_res,               # 输出分辨率
            "output_grid_size": [int(source_ds.RasterXSize),  # 输出列数
                                 int(source_ds.RasterYSize)],  # 输出行数
            "geometry": "native grid; RPC orthorectification deferred to step 3",  # 几何说明
            "target_crs_epsg": (georef.get("target_crs_epsg")  # 步骤1目标EPSG
                                if georef else source_epsg),  # 无元数据用栅格值
            "target_res": (georef.get("target_res") if georef else source_res),  # 目标分辨率
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),  # 生成时间
        }
        write_parameters(payload, final_params_json, final_params_txt, logger)  # 写参数记录
        main_complete = True                             # 标记主成果已完成
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

        # ⑦ 生成带坐标预览和质检图
        current_stage = 7                                # 进入第7阶段
        stage_description = "生成降采样预览和质检图"       # 阶段说明
        stage_started = log_stage_start(logger, current_stage, 7,  # 阶段开始计时
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
            try:                                         # 质检图生成可能失败
                write_quality_figure(final_tif, wavelengths, final_qc, scene,  # 生成质检图
                                     params_summary_text(payload), logger,  # 参数摘要
                                     band_stats=stats, preview_tif=preview_tif)  # 统计与预览
            except Exception:                            # 绘图失败
                logger.warning("质检图生成失败: %s", final_qc, exc_info=True)  # 告警
        else:                                            # 已关闭质检
            logger.info("按配置跳过质检图生成")           # 记录跳过
        log_stage_done(logger, current_stage, 7, stage_description, stage_started)  # 阶段结束

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
        for path in (tmp_tif, tmp_coef, tmp_wl):         # 清理临时文件
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
        logger.error("场景处理在步骤 %d/7 失败：%s", current_stage,  # 记录失败阶段
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
