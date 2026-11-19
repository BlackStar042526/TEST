# -*- coding: utf-8 -*-
"""
00_extract_georeferenced_jpg.py
提取快视图 JPG，将 SHP 边界绘制进 JPG，并为图片赋予地理参考（v1.0.2）

处理顺序：定位 JPG/SHP → 提取/复制 → 读 SHP 范围与投影 → 写 .jgw/.prj
→ 绘制边界 → 输出带坐标的 GeoTIFF。

本版在 v1.0.1 基础上补充逐行中文备注，算法、输出文件名与格式完全不变。
"""
import argparse
import glob
import logging
import math
import os
import sys
import time
import tarfile
import shutil

from PIL import Image, ImageDraw
from osgeo import gdal, ogr

gdal.UseExceptions()                          # GDAL 出错时抛异常，而不是返回 None
gdal.PushErrorHandler("CPLQuietErrorHandler")  # 屏蔽 GDAL 非致命错误的标准错误输出
ogr.UseExceptions()                           # OGR 同步开启异常模式

# =================================================================
# 配置区（按需修改）
# =================================================================
SRC_DIR = r"E:\ZhouYang\test\ZY1E_AHSI"           # 默认源目录（L1A 压缩包所在处）
OUT_ROOT = r"E:\ZhouYang\test\ZY1E_Georef_JPG"    # 默认输出根目录

# 支持的压缩包格式
TAR_PATTERNS = ("*.tar", "*.tar.gz", "*.tgz")     # 扫描压缩包时使用的通配模式

# 输出设置
OVERWRITE = False                       # 默认跳过已有成果，加 --overwrite 才重做
LOG_NAME = "zy1e.georef_jpg_shp"        # 日志器名称，避免与第三方库重名
BOUNDARY_COLOR = (255, 0, 0)            # 边界线颜色：纯红
BOUNDARY_WIDTH = 2                      # 边界线宽（像素）
JPEG_QUALITY = 95                       # 输出快视图 JPG 的质量
GEOTIFF_JPEG_QUALITY = 90               # GeoTIFF 内部 JPEG 压缩质量
BOUNDARY_RENDER_MARKER = ".boundary_rendered"  # 边界已绘制的标记文件名

# =================================================================
# 通用工具
# =================================================================
def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(description="ZY1E 利用SHP提取并赋予JPG坐标")  # 创建参数解析器
    parser.add_argument("--src", default=SRC_DIR)          # 源目录
    parser.add_argument("--out", default=OUT_ROOT)         # 输出根目录
    parser.add_argument("--pattern", action="append", default=None)  # 自定义扫描通配符，可多次指定
    parser.add_argument("--scene", action="append", default=None)    # 只处理匹配的景名，支持通配符
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE)  # 覆盖已有成果
    parser.add_argument("--limit", type=int, default=0)    # 仅处理前 N 景，0 表示不限
    parser.add_argument("--dry-run", action="store_true")  # 只打印计划不落盘
    parser.add_argument("--log", default=None)             # 自定义日志文件路径
    parser.add_argument("--verbose", action="store_true")  # 控制台输出 DEBUG 级日志
    return parser.parse_args(argv)                         # 解析并返回命名空间

def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)      # 取得具名日志器
    logger.setLevel(logging.DEBUG)            # 记录器总级别设为 DEBUG
    logger.handlers.clear()                   # 清空旧处理器，防止重复输出
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")  # 统一日志格式
    
    console = logging.StreamHandler(stream=sys.stdout)              # 控制台处理器
    console.setLevel(logging.DEBUG if verbose else logging.INFO)    # 按 --verbose 决定级别
    console.setFormatter(fmt)                                       # 套用统一格式
    logger.addHandler(console)                                      # 挂到日志器
    
    if log_path:                                                    # 指定了日志文件才写文件
        try:
            os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)  # 确保日志目录存在
            file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")  # 追加写，避免覆盖历史
            file_handler.setLevel(logging.DEBUG)                    # 文件始终记录最详细级别
            file_handler.setFormatter(fmt)                          # 套用统一格式
            logger.addHandler(file_handler)                         # 挂到日志器
        except OSError:
            logger.exception("无法创建日志文件: %s", log_path)      # 记录失败堆栈
            raise                                                   # 无法写日志则直接终止
    logger.propagate = False                                        # 不向根日志器传播，避免重复打印
    return logger                                                   # 返回配置好的日志器

def discover_scenes(src_dir, patterns):
    """扫描源目录：返回压缩包路径，以及同时含 .jpg 和 .shp 的目录。"""
    if not os.path.isdir(src_dir):                                  # 源目录必须存在
        raise NotADirectoryError("源目录不存在: %s" % src_dir)      # 不存在直接报错
    found = set()                                                   # 用集合收集，自动去重
    for pattern in patterns:                                        # 逐个通配模式扫描
        found.update(glob.glob(os.path.join(src_dir, pattern), recursive=True))      # 当前层匹配
        found.update(glob.glob(os.path.join(src_dir, "**", pattern), recursive=True))  # 递归子目录匹配
    for current_dir, directories, filenames in os.walk(src_dir):    # 递归遍历全部目录
        directories[:] = [name for name in directories if not name.startswith(".")]  # 跳过隐藏目录
        extensions = {os.path.splitext(name)[1].lower() for name in filenames}       # 当前目录的后缀集合（小写）
        if ".jpg" in extensions and ".shp" in extensions:           # 同时含 JPG 与 SHP 才算候选景
            found.add(current_dir)                                  # 记录该目录
    return sorted(f for f in found if os.path.isfile(f) or os.path.isdir(f))  # 过滤失效路径并按名称排序

def get_scene_name(path):
    """由压缩包或目录路径推导景名（去掉 .tar/.tar.gz/.tgz 等后缀）。"""
    name = os.path.basename(path)           # 取末级文件名/目录名
    if os.path.isdir(path):
        return name                         # 目录场景直接使用目录名
    lower = name.lower()                    # 转小写便于后缀比对
    for ext in (".tar.gz", ".tgz", ".tar"): 
        if lower.endswith(ext): return name[: -len(ext)]  # 命中压缩包后缀则截掉
    return os.path.splitext(name)[0]        # 其余情况仅去掉扩展名

def extract_tar_member(tar, member_name, out_path):
    """从已打开的 tar 包中流式提取单个成员到 out_path。"""
    source = tar.extractfile(member_name)                       # 取得成员的文件对象
    if source is None:                                          # 目录或无效成员会返回 None
        raise OSError("无法从 TAR 包读取文件: %s" % member_name)  # 明确报错便于定位
    with source, open(out_path, "wb") as destination:           # 同时打开源与目标
        shutil.copyfileobj(source, destination)                 # 分块流式复制，省内存


def matching_sidecars(files, shp_name):
    """返回与 .shp 同名、同目录的标准 sidecar（.dbf/.shx/.prj/.cpg）。"""
    shp_stem = os.path.splitext(os.path.basename(shp_name))[0].casefold()  # SHP 主文件名（忽略大小写）
    shp_dir = os.path.dirname(shp_name).casefold()                         # SHP 所在目录（忽略大小写）
    sidecar_extensions = {".dbf", ".shx", ".prj", ".cpg"}                  # 需要一并复制的附属文件后缀
    return [
        path for path in files                                             # 遍历候选文件列表
        if os.path.dirname(path).casefold() == shp_dir                     # 必须与 SHP 同目录
        and os.path.splitext(os.path.basename(path))[0].casefold() == shp_stem  # 主名一致
        and os.path.splitext(path)[1].lower() in sidecar_extensions        # 后缀属于 sidecar 集合
    ]


def burn_shapefile_boundary(jpg_path, shp_path, extent, logger):
    """把 SHP 的线/面边界按 extent 线性映射到 JPG 像素坐标，用红线绘制。

    像素坐标换算（extent 为 SHP 外接矩形 minX/maxX/minY/maxY）：
        px = (x - min_x) / xres
        py = (max_y - y) / yres      # 图像 y 轴向下，故用 max_y 减
    结果先写临时 JPG 再 os.replace 覆盖原图，保证原子性。
    """
    min_x, max_x, min_y, max_y = extent            # 解包 SHP 外接矩形
    image = Image.open(jpg_path).convert("RGB")    # 读入原图并统一为 RGB
    width, height = image.size                     # 记录像素宽高
    xres = (max_x - min_x) / width                 # 每个像素对应的 X 方向地理长度
    yres = (max_y - min_y) / height                # 每个像素对应的 Y 方向地理长度
    draw = ImageDraw.Draw(image)                   # 创建绘图句柄

    vector = ogr.Open(shp_path, 0)                 # 以只读方式打开 SHP
    if vector is None:
        raise RuntimeError("绘制边界时无法重新打开 SHP: %s" % shp_path)  # 打不开则终止
    layer = vector.GetLayer()                      # 取第一个图层
    if layer is None:
        vector = None                              # 先释放数据源
        raise RuntimeError("SHP 数据源中没有可绘制的图层")  # 无图层无法绘制
    feature_count = 0                              # 统计实际绘制的几何部件数
    try:
        for feature in layer:                      # 逐个要素遍历
            geometry = feature.GetGeometryRef()    # 取要素几何
            if geometry is None or geometry.IsEmpty():  # 空几何跳过
                continue
            boundary = geometry.Boundary() if geometry.GetGeometryName().upper() in (
                "POLYGON", "MULTIPOLYGON"
            ) else geometry.Clone()                # 面取边界，线直接克隆
            if boundary is None or boundary.IsEmpty():  # 边界为空也跳过
                continue

            parts = [boundary]                     # 用栈展开多部件几何
            while parts:
                part = parts.pop()                 # 取出一个部件
                if part.GetGeometryCount():        # 仍是集合则继续拆分
                    parts.extend(
                        part.GetGeometryRef(index)
                        for index in range(part.GetGeometryCount())
                    )
                    continue                       # 拆分后不参与绘制
                if part.GetPointCount() < 2:       # 少于两点无法成线
                    continue
                pixels = [                         # 地理坐标 → 像素坐标
                    (
                        (part.GetX(index) - min_x) / xres,      # X 平移后除以分辨率
                        (max_y - part.GetY(index)) / yres,      # Y 翻转（图像向下）后换算
                    )
                    for index in range(part.GetPointCount())    # 遍历该部件全部节点
                ]
                draw.line(
                    pixels,                        # 像素折线
                    fill=BOUNDARY_COLOR,           # 线色（红）
                    width=BOUNDARY_WIDTH,          # 线宽
                    joint="curve",                 # 拐角圆滑连接
                )
                feature_count += 1                 # 绘制计数加一
    finally:
        vector = None                              # 无论成败都释放数据源

    if feature_count == 0:                         # 一条边界都没画出来说明数据有问题
        raise ValueError("SHP 中没有可绘制的线或面边界")

    temp_path = jpg_path + ".boundary.tmp.jpg"     # 临时文件名，避免写坏原图
    try:
        image.save(temp_path, format="JPEG", quality=JPEG_QUALITY)  # 先写临时 JPG
        os.replace(temp_path, jpg_path)            # 原子替换，覆盖原图
    except Exception:
        if os.path.exists(temp_path):              # 失败时清理残留临时文件
            os.remove(temp_path)
        raise                                      # 继续向上抛错
    logger.info(
        "  - 已将 SHP 边界绘制到 JPG（红色，%d 像素宽）",
        BOUNDARY_WIDTH,                            # 记录实际线宽
    )


def create_georeferenced_geotiff(jpg_path, tif_path, extent, projection_wkt,
                                 logger):
    """给已绘制边界的 JPG 写入仿射变换与投影，输出 JPEG 压缩的 GeoTIFF。

    geotransform = (min_x, 像元宽, 0, max_y, 0, -像元高)
    其中 像元宽=(max_x-min_x)/宽，像元高=(max_y-min_y)/高；写出后重新打开
    校验波段数、尺寸、投影与仿射变换，全部通过才原子替换为正式文件。
    """
    min_x, max_x, min_y, max_y = extent                   # 解包 SHP 外接矩形
    source = gdal.Open(jpg_path, gdal.GA_ReadOnly)        # 只读打开已绘制边界的 JPG
    if source is None:
        raise RuntimeError("无法打开已绘制边界的 RGB 预览图: %s" % jpg_path)  # 打不开直接报错
    if source.RasterCount != 3:                           # 必须是三波段彩色图
        source = None                                     # 释放数据源
        raise ValueError("RGB 预览图应为 3 个波段，实际波段数不符")

    width, height = source.RasterXSize, source.RasterYSize  # 记录影像宽高
    geotransform = (
        min_x,                                            # 左上角 X
        (max_x - min_x) / width,                          # 像元宽
        0.0,                                              # 无旋转项
        max_y,                                            # 左上角 Y
        0.0,                                              # 无旋转项
        -(max_y - min_y) / height,                        # 像元高（负值表示向下）
    )
    source.SetGeoTransform(geotransform)                  # 写入仿射变换
    source.SetProjection(projection_wkt)                  # 写入投影 WKT
    temp_path = tif_path + ".tmp.tif"                     # 临时 GeoTIFF 路径
    driver = gdal.GetDriverByName("GTiff")                # 获取 GeoTIFF 驱动
    if driver is None:
        source = None                                     # 释放数据源
        raise RuntimeError("GDAL 未提供 GeoTIFF 驱动")

    try:
        output = driver.CreateCopy(
            temp_path,                                    # 先写临时文件
            source,                                       # 带仿射与投影的内存数据集
            strict=1,                                     # 严格模式，出错即失败
            options=[
                "COMPRESS=JPEG",                          # 采用 JPEG 压缩
                "JPEG_QUALITY=%d" % GEOTIFF_JPEG_QUALITY, # 压缩质量
                "PHOTOMETRIC=YCBCR",                      # 彩色图用 YCbCr 提高压缩率
                "TILED=YES",                              # 分块存储便于局部读取
            ],
        )
        if output is None:
            raise RuntimeError("创建压缩 GeoTIFF 失败: %s" % temp_path)  # 创建失败报错
        output.FlushCache()                               # 强制落盘
        output = None                                     # 关闭写句柄
        source = None                                     # 关闭源数据集

        check = gdal.Open(temp_path, gdal.GA_ReadOnly)    # 重新打开校验
        if check is None:
            raise RuntimeError("生成的 GeoTIFF 无法重新打开")  # 打不开说明文件损坏
        try:
            if (check.RasterCount != 3
                    or check.RasterXSize != width
                    or check.RasterYSize != height):      # 波段数或尺寸必须与 JPG 一致
                raise RuntimeError("GeoTIFF 波段数或图像尺寸与预览 JPG 不一致")
            if not check.GetProjection():                 # 投影必须写入成功
                raise RuntimeError("GeoTIFF 未写入坐标参考系")
            if check.GetGeoTransform() != geotransform:   # 仿射变换必须与预期一致
                raise RuntimeError("GeoTIFF 仿射变换与预期不一致")
        finally:
            check = None                                  # 校验完关闭数据集
        os.replace(temp_path, tif_path)                   # 原子替换为正式成果
    except Exception:
        source = None                                     # 异常时释放源数据集
        if os.path.exists(temp_path):                     # 清理临时文件
            try:
                os.remove(temp_path)
            except OSError:
                logger.exception("清理临时 GeoTIFF 失败: %s", temp_path)  # 清理失败仅记录
        raise                                             # 继续向上抛错

    logger.info(
        "  - 已生成带坐标系的 RGB GeoTIFF: %s (JPEG 压缩)",
        os.path.basename(tif_path),                       # 只打印文件名
    )

# =================================================================
# 单景处理逻辑
# =================================================================
def process_scene(scene_path, out_root, args, logger):
    """处理单景：①定位文件 ②提取/复制 ③读 SHP 范围 ④算分辨率

    ⑤写 .jgw ⑥写 .prj ⑦绘边界并输出 GeoTIFF。
    单景异常在此捕获并记录堆栈，返回状态字典，不中断批处理。
    """
    scene = get_scene_name(scene_path)              # 由路径推导景名
    scene_dir = os.path.join(out_root, scene)       # 该景输出目录

    logger.info("开始处理场景: %s", scene)          # 记录开始
    started = time.time()                           # 整景计时起点
    stage_started = started                         # 分阶段计时起点

    try:
        os.makedirs(scene_dir, exist_ok=True)       # 建立输出目录
        jpg_member = None                           # 压缩包内 JPG 成员名
        shp_member = None                           # 压缩包内 SHP 成员名
        shp_related = []                            # SHP 附属文件列表
        is_tar = not os.path.isdir(scene_path)      # 非目录即视为压缩包
        
        out_jpg_path = None                         # 输出 JPG 路径
        out_shp_path = None                         # 输出 SHP 路径

        # ① 定位文件：压缩包内按后缀查找 JPG/SHP；目录则直接列文件
        if is_tar:
            with tarfile.open(scene_path, "r:*") as tar:          # 自动识别压缩格式
                members = [m for m in tar.getmembers() if m.isfile()]  # 仅保留普通文件
                jpg_members = sorted(m.name for m in members
                                     if m.name.lower().endswith(".jpg"))  # 包内全部 JPG
                shp_members = sorted(m.name for m in members
                                     if m.name.lower().endswith(".shp"))  # 包内全部 SHP
                if not jpg_members:
                    raise FileNotFoundError("包内未找到 JPG")          # 缺 JPG 无法继续
                if not shp_members:
                    raise FileNotFoundError("包内未找到 SHP 矢量边界")  # 缺 SHP 无法继续
                if len(jpg_members) > 1 or len(shp_members) > 1:       # 多份时给出警告
                    logger.warning("包内 JPG/SHP 不止一个，按排序选择首个文件。")
                jpg_member = jpg_members[0]                            # 取排序后的第一个 JPG
                shp_member = shp_members[0]                            # 取排序后的第一个 SHP
                shp_dir = os.path.dirname(shp_member)                  # SHP 所在包内目录
                shp_related = matching_sidecars(
                    [m.name for m in members
                     if os.path.dirname(m.name) == shp_dir],           # 仅同目录成员参与匹配
                    shp_member,
                )

                out_jpg_path = os.path.join(scene_dir, os.path.basename(jpg_member))  # 输出 JPG
                out_shp_path = os.path.join(scene_dir, os.path.basename(shp_member))  # 输出 SHP
                out_tif_path = os.path.splitext(out_jpg_path)[0] + ".tif"             # 输出 GeoTIFF
                jgw_path = os.path.splitext(out_jpg_path)[0] + ".jgw"                 # 输出坐标文件
                jpg_prj_path = os.path.splitext(out_jpg_path)[0] + ".prj"             # 输出投影文件
                expected_outputs = [                                                  # 成果清单
                    out_jpg_path, out_shp_path, out_tif_path, jgw_path, jpg_prj_path,
                    os.path.join(scene_dir, BOUNDARY_RENDER_MARKER),                  # 边界标记
                    *(os.path.join(scene_dir, os.path.basename(path))
                      for path in shp_related),                                       # 附属文件
                ]
                if not args.overwrite and all(os.path.isfile(path)
                                              for path in expected_outputs):          # 成果齐全则跳过
                    logger.info("跳过（成果已存在）: %s", scene)
                    return {"scene": scene, "status": "跳过", "elapsed": 0.0}         # 提前返回

                logger.info("  - 正在提取 JPG 和 SHP 矢量...")
                extract_tar_member(tar, jpg_member, out_jpg_path)                     # 提取 JPG
                extract_tar_member(tar, shp_member, out_shp_path)                     # 提取 SHP
                for rel in shp_related:                                               # 逐个提取附属文件
                    extract_tar_member(tar, rel, os.path.join(scene_dir, os.path.basename(rel)))
        else:
            source_files = [                                                          # 目录内全部文件
                os.path.join(scene_path, name) for name in os.listdir(scene_path)
                if os.path.isfile(os.path.join(scene_path, name))
            ]
            jpg_files = sorted(path for path in source_files
                               if path.lower().endswith(".jpg"))                      # 目录内 JPG
            shp_files = sorted(path for path in source_files
                               if path.lower().endswith(".shp"))                      # 目录内 SHP
            if not jpg_files or not shp_files:                                        # 缺文件直接报错
                raise FileNotFoundError("文件夹内缺少 JPG 或 SHP")
            if len(jpg_files) > 1 or len(shp_files) > 1:                              # 多份时警告
                logger.warning("场景目录内 JPG/SHP 不止一个，按排序选择首个文件。")

            out_jpg_path = os.path.join(scene_dir, os.path.basename(jpg_files[0]))    # 输出 JPG
            out_shp_path = os.path.join(scene_dir, os.path.basename(shp_files[0]))    # 输出 SHP
            out_tif_path = os.path.splitext(out_jpg_path)[0] + ".tif"                 # 输出 GeoTIFF
            shp_related = matching_sidecars(source_files, shp_files[0])               # 匹配附属文件
            jgw_path = os.path.splitext(out_jpg_path)[0] + ".jgw"                     # 输出坐标文件
            jpg_prj_path = os.path.splitext(out_jpg_path)[0] + ".prj"                 # 输出投影文件
            expected_outputs = [                                                      # 成果清单
                out_jpg_path, out_shp_path, out_tif_path, jgw_path, jpg_prj_path,
                os.path.join(scene_dir, BOUNDARY_RENDER_MARKER),                      # 边界标记
                *(os.path.join(scene_dir, os.path.basename(path))
                  for path in shp_related),                                           # 附属文件
            ]
            if not args.overwrite and all(os.path.isfile(path)
                                          for path in expected_outputs):              # 成果齐全则跳过
                logger.info("跳过（成果已存在）: %s", scene)
                return {"scene": scene, "status": "跳过", "elapsed": 0.0}             # 提前返回

            shutil.copy2(jpg_files[0], out_jpg_path)                                  # 复制 JPG（保留时间）
            shutil.copy2(shp_files[0], out_shp_path)                                  # 复制 SHP
            for sidecar in shp_related:                                               # 逐个复制附属文件
                shutil.copy2(sidecar, scene_dir)

        logger.info("  - 文件准备完成，用时 %.1f s", time.time() - stage_started)  # 阶段耗时
        stage_started = time.time()                                               # 重置阶段计时

        # ② 读取 SHP 获取坐标范围和投影
        logger.info("  - 正在解析 SHP 获取空间坐标...")
        shp_ds = ogr.Open(out_shp_path, 0)                    # 只读打开 SHP
        if shp_ds is None: raise RuntimeError(f"无法打开矢量: {out_shp_path}")  # 打不开直接失败
        
        layer = shp_ds.GetLayer()                             # 取第一个图层
        if layer is None:
            raise RuntimeError("SHP 数据源中没有可读取的图层")  # 无图层无法取范围
        extent = layer.GetExtent(force=1) # (minX, maxX, minY, maxY)  # force=1 强制精确计算
        spatial_ref = layer.GetSpatialRef()                   # 取坐标系
        if extent is None or not all(math.isfinite(value) for value in extent):  # 范围必须有效
            raise ValueError("SHP 的空间范围为空或包含无效坐标")
        
        if spatial_ref:
            # 转换为 WKT 以写入 .prj
            prj_wkt = spatial_ref.ExportToWkt()               # 导出投影 WKT
            # 检查是不是地理坐标 (经纬度)
            is_geo = spatial_ref.IsGeographic()               # 是否地理坐标系
            unit_name = "度" if is_geo else "米"               # 根据类型给单位名
            logger.info("  - SHP 坐标系: %s (单位: %s)", spatial_ref.GetName(), unit_name)  # 打印坐标系
        else:
            raise RuntimeError("SHP 文件缺少 .prj 投影信息，无法确定坐标系统")  # 无投影则失败
            
        shp_ds = None # 关闭                                        # 释放矢量数据源

        # ③ 获取 JPG 像素尺寸
        jpg_ds = gdal.Open(out_jpg_path)                      # 打开 JPG 读尺寸
        if jpg_ds is None: raise RuntimeError("无法打开 JPG")  # 打不开直接失败
        width_jpg, height_jpg = jpg_ds.RasterXSize, jpg_ds.RasterYSize  # 记录宽高
        jpg_ds = None                                         # 释放数据集

        # ④ 计算像元分辨率 (根据 SHP 范围计算)
        min_x, max_x, min_y, max_y = extent                   # 解包范围
        if (max_x <= min_x or max_y <= min_y
                or width_jpg <= 0 or height_jpg <= 0):        # 范围或尺寸非法
            raise ValueError(f"SHP 坐标范围异常: {extent}")

        xres_jpg = (max_x - min_x) / width_jpg                # X 方向像元分辨率
        yres_jpg = (max_y - min_y) / height_jpg               # Y 方向像元分辨率
        logger.info(
            "  - 按 SHP 外接矩形计算 JPG 地理参考（假设 JPG 覆盖整个矩形且无旋转）"
        )

        # ⑤ 生成 .jgw (World File)：A=像元宽 D=旋转 B=旋转 E=-像元高 C/F=左上像元中心
        with open(jgw_path, "w") as f:                        # 新建坐标文件
            f.write(f"{xres_jpg}\n")  # A: X 分辨率
            f.write("0.0\n")          # D: 旋转
            f.write("0.0\n")          # B: 旋转
            f.write(f"{-yres_jpg}\n") # E: Y 分辨率 (负值)
            f.write(f"{min_x + xres_jpg / 2.0:.17g}\n") # C: 左上像元中心 X
            f.write(f"{max_y - yres_jpg / 2.0:.17g}\n") # F: 左上像元中心 Y
        logger.info("  - 成功生成坐标: %s", os.path.basename(jgw_path))  # 记录坐标文件

        # ⑥ 生成 .prj (与 JPG 同名)，内容为 SHP 的 WKT 投影
        prj_path = jpg_prj_path                               # 投影文件路径
        with open(prj_path, "w", encoding="utf-8") as f:      # 新建投影文件
            f.write(prj_wkt)                                  # 写入 WKT
        logger.info("  - 成功生成投影: %s", os.path.basename(prj_path))  # 记录投影文件
        logger.info("  - 地理参考写入完成，用时 %.1f s", time.time() - stage_started)  # 阶段耗时

        # ⑦ 绘制 SHP 边界并输出带地理参考的 RGB GeoTIFF
        burn_shapefile_boundary(out_jpg_path, out_shp_path, extent, logger)  # 先画红线边界
        create_georeferenced_geotiff(
            out_jpg_path, out_tif_path, extent, prj_wkt, logger              # 再写坐标并转 GeoTIFF
        )
        with open(os.path.join(scene_dir, BOUNDARY_RENDER_MARKER), "w",
                  encoding="utf-8") as marker:                               # 写边界标记文件
            marker.write("SHP boundary burned into JPG\n")

        logger.info("场景处理完成: %s (耗时 %.1fs)", scene, time.time() - started)  # 整景耗时
        return {"scene": scene, "status": "成功", "elapsed": time.time() - started}  # 成功状态

    except Exception as e:
        logger.exception("处理失败: %s", e)                                            # 记录失败堆栈
        return {"scene": scene, "status": "失败", "elapsed": time.time() - started}    # 失败状态

# =================================================================
# 批处理主流程
# =================================================================
def main(argv=None):
    """批处理入口：扫描 → （dry-run 仅打印计划）→ 逐景处理 → 汇总。"""
    args = parse_args(argv)                                  # 解析命令行
    out_root = os.path.abspath(args.out)                     # 输出根目录绝对路径
    os.makedirs(out_root, exist_ok=True)                     # 确保输出根目录存在
    log_path = args.log or os.path.join(out_root, "step0_georef_jpg_shp.log")  # 默认日志位置
    logger = setup_logging(log_path, args.verbose)           # 初始化日志

    logger.info("=" * 60)
    logger.info("ZY1E AHSI 步骤0：快视图 JPG 边界绘制与地理参考")
    logger.info("源目录  : %s", os.path.abspath(args.src))   # 打印源目录
    logger.info("输出目录: %s", out_root)                    # 打印输出目录
    logger.info("日志文件: %s", log_path)                    # 打印日志路径
    logger.info("=" * 60)

    patterns = tuple(args.pattern) if args.pattern else TAR_PATTERNS  # 命令行优先，否则用默认
    try:
        scenes = discover_scenes(os.path.abspath(args.src), patterns)  # 扫描候选景
    except Exception as e:
        logger.error("扫描失败: %s", e)                       # 扫描失败记录原因
        return 2                                             # 返回错误码

    if args.scene:
        import fnmatch                                        # 局部导入，仅此分支使用
        scenes = [p for p in scenes if any(fnmatch.fnmatch(get_scene_name(p), pat) for pat in args.scene)]  # 按景名过滤
    if args.limit > 0: scenes = scenes[:args.limit]           # 应用数量上限
    if not scenes: 
        logger.error("未找到任何包含 JPG 和 SHP 的文件夹或压缩包")  # 无候选景
        return 2                                             # 返回错误码

    if args.dry_run:
        for idx, p in enumerate(scenes, 1):                   # 仅打印处理计划
            logger.info("[%d/%d] %s -> %s", idx, len(scenes), get_scene_name(p), os.path.join(out_root, get_scene_name(p)))
        logger.info("dry-run 完成：共 %d 景，未处理文件", len(scenes))  # 汇总计划数
        return 0                                             # dry-run 正常退出

    results = []                                             # 收集每景结果
    started = time.time()                                    # 批处理计时起点
    
    for idx, scene_path in enumerate(scenes, start=1):       # 逐景处理
        logger.info("="*50)
        logger.info("正在处理第 %d/%d 景数据", idx, len(scenes))  # 当前进度
        results.append(process_scene(scene_path, out_root, args, logger))  # 单景结果入列

    success = sum(1 for r in results if r["status"] == "成功")              # 成功数
    skipped = sum(1 for r in results if r["status"].startswith("跳过"))     # 跳过数
    failed = sum(1 for r in results if r["status"] == "失败")               # 失败数
    
    logger.info("="*50)
    logger.info("批处理结束: 成功 %d, 跳过 %d, 失败 %d，总用时 %.1f s", success, skipped, failed, time.time() - started)  # 汇总
    for result in results:
        logger.info("  [%s] %s（%.1f s）", result["status"], result["scene"], result["elapsed"])  # 逐景结果
    if failed:
        logger.error("有 %d 景处理失败，详见日志: %s", failed, log_path)  # 有失败则提示
        return 1                                             # 存在失败返回 1
    logger.info("全部完成，成果位于: %s", out_root)           # 全部成功
    return 0                                                 # 成功返回 0

if __name__ == "__main__":
    sys.exit(main())                                         # 以 main 返回值作为进程退出码
