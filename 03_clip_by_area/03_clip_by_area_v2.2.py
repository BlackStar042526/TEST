#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Clip step-2 reflectance rasters by research-area shapefiles (v2.2).

按研究区矢量裁剪步骤2 大气校正影像：不重投影、不重采样，保留源 CRS、
像元大小、波段顺序与像元对齐。流程：扫描影像/矢量 → 求交集像元窗口
→ gdal.Warp 按 cutline 裁剪 → 回写波段元数据 → 校验网格后原子替换。

输出直接写到 <output-root>/<area-name>/<scene-name>_cal_badband_ortho_atmcorr_clip.tif，
不建立中间区域子目录。

Examples:
    python 03_clip_by_area_v2.2.py --dry-run
    python 03_clip_by_area_v2.2.py
    python 03_clip_by_area_v2.2.py --area JLHK --scene "*20250121*"

本版仅补充逐行中文备注，算法与输出不变。
"""

import argparse
import fnmatch
import glob
import logging
import math
import os
import re
import sys
import time

try:
    from osgeo import gdal, ogr, osr
except ImportError as exc:
    raise SystemExit(
        "当前 Python 解释器缺少 GDAL（osgeo）：%s\n"
        "请使用安装了 GDAL 的 Python 解释器运行本脚本。" % exc)

gdal.UseExceptions()   # GDAL 出错即抛异常，便于 try/finally 清理
ogr.UseExceptions()    # OGR 同样开启异常模式

IN_ROOT_STEP2 = r"E:\ZhouYang\test\ZY1E_output_step2"   # 步骤2 大气校正成果根目录
OUT_ROOT_STEP3 = r"E:\ZhouYang\test\ZY1E_output_step3"  # 步骤3 裁剪成果根目录
AREA_ROOT = r"E:\ZhouYang\test\Hyperarea2"              # 研究区矢量根目录

IN_REFLECTANCE_SUFFIX = "_cal_badband_ortho_atmcorr.tif"  # 新格式反射率成果后缀
LEGACY_REFLECTANCE_SUFFIXES = (                           # 需兼容的旧成果后缀
    "_ortho_atmcorr.tif",
    "reflectance_atm.tif",
)
OUTPUT_SUFFIX = "_cal_badband_ortho_atmcorr_clip.tif"     # 裁剪成果文件名后缀
OUTPUT_NODATA = -9999.0                                   # 无 NoData 时的默认值
LOG_NAME = "zy1e.step3.clip"                              # 日志器名称，避免重名
GTIFF_OPTIONS = [                                         # GeoTIFF 写出选项
    "BIGTIFF=IF_SAFER",                                   # 超 4GB 自动切 BigTIFF
    "TILED=YES",                                          # 分块存储便于局部读取
    "INTERLEAVE=BAND",                                    # 按波段交错排列
    "COMPRESS=DEFLATE",                                   # 采用 DEFLATE 无损压缩
    "PREDICTOR=3",                                        # 浮点预测器提升压缩率
    "ZLEVEL=6",                                           # 压缩级别 6，兼顾速度
]


def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(
        description="按研究区矢量裁剪步骤2大气校正影像（不重投影、不重采样）")  # 创建解析器
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP2,
                        help="步骤2成果根目录（默认：%(default)s）")           # 输入根目录
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP3,
                        help="裁剪成果根目录（默认：%(default)s）")           # 输出根目录
    parser.add_argument("--areas", dest="area_root", default=AREA_ROOT,
                        help="研究区矢量根目录（默认：%(default)s）")         # 矢量根目录
    parser.add_argument("--area", action="append", default=None,
                        help="只处理匹配的区域名，可多次指定，支持通配符")    # 区域筛选
    parser.add_argument("--scene", action="append", default=None,
                        help="只处理匹配的场景名，可多次指定，支持通配符")    # 场景筛选
    parser.add_argument("--overwrite", action="store_true",
                        help="覆盖已有裁剪成果（默认跳过）")                  # 覆盖开关
    parser.add_argument("--dry-run", action="store_true",
                        help="只显示处理计划，不写出影像")                    # 仅打印计划
    parser.add_argument("--log", default=None,
                        help="日志文件（默认：<输出根目录>\\clip_batch.log）")  # 日志路径
    parser.add_argument("--verbose", action="store_true",
                        help="输出调试日志")                                  # 调试日志
    return parser.parse_args(argv)                                            # 解析并返回


def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)               # 取得具名日志器
    logger.setLevel(logging.DEBUG)                     # 记录器总级别设为 DEBUG
    for existing_handler in list(logger.handlers):     # 遍历已有处理器
        logger.removeHandler(existing_handler)         # 先摘除
        existing_handler.close()                       # 再关闭，避免重复输出
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")                   # 统一时间格式
    console = logging.StreamHandler(stream=sys.stdout)  # 控制台处理器
    console.setLevel(logging.DEBUG if verbose else logging.INFO)  # 按 --verbose 定级
    console.setFormatter(formatter)                    # 套用统一格式
    logger.addHandler(console)                         # 挂到日志器
    if log_path:                                       # 指定了日志文件才写文件
        parent = os.path.dirname(os.path.abspath(log_path))  # 日志文件所在目录
        if parent:                                     # 目录部分非空才创建
            os.makedirs(parent, exist_ok=True)         # 缺目录自动补齐
        handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")  # 追加写
        handler.setLevel(logging.DEBUG)                # 文件始终记录最详细级别
        handler.setFormatter(formatter)                # 套用统一格式
        logger.addHandler(handler)                     # 挂到日志器
    logger.propagate = False                           # 不向根日志器传播
    return logger                                      # 返回配置好的日志器


def normalize_area_name(name):
    """规范化区域名：去序号前缀、空白与非法字符，作为输出目录名。"""
    value = os.path.basename(os.path.normpath(str(name)))  # 只保留末级名称
    value = re.sub(r"^\s*\d+\s*[\.\-_]?\s*", "", value)    # 去掉“01_”这类序号前缀
    value = re.sub(r"\s+", "", value)                      # 去掉全部空白字符
    value = re.sub(r'[<>:"/\\|?*]', "_", value)            # 路径非法字符换成下划线
    return value.strip(" ._")                              # 去掉首尾的点和下划线


def matches_filters(name, filters):
    """名称是否匹配任一通配符筛选条件；filters 为空表示全部匹配。"""
    if not filters:                                        # 未给筛选条件即全部命中
        return True
    return any(fnmatch.fnmatchcase(name.lower(), pattern.lower())
               for pattern in filters)                     # 忽略大小写逐一比对


def discover_scenes(in_root, filters):
    """扫描步骤2 成果（新后缀优先、兼容旧后缀），返回场景列表。"""
    if not os.path.isdir(in_root):                         # 输入目录必须存在
        raise NotADirectoryError("步骤2成果目录不存在：%s" % in_root)
    paths_by_scene = {}                                    # 场景目录 → 影像路径
    for path in sorted(glob.glob(
            os.path.join(in_root, "**", "*" + IN_REFLECTANCE_SUFFIX),
            recursive=True)):                              # 递归收集新格式成果
        scene_dir = os.path.dirname(path)                  # 影像所在场景目录
        paths_by_scene[scene_dir] = path                   # 记录该场景
    for suffix in LEGACY_REFLECTANCE_SUFFIXES:             # 再补扫兼容的旧后缀
        legacy_name = suffix if suffix.endswith(".tif") else suffix + ".tif"  # 补齐扩展名
        pattern_name = ("*" + suffix if suffix.startswith("_") else legacy_name)  # 拼通配模式
        for path in sorted(glob.glob(
                os.path.join(in_root, "**", pattern_name), recursive=True)):  # 递归扫描
            paths_by_scene.setdefault(os.path.dirname(path), path)  # 新格式优先，旧格式补缺

    scenes = []
    for scene_dir, path in sorted(paths_by_scene.items()):  # 按目录名排序遍历
        scene_name = os.path.basename(scene_dir)            # 场景名取目录名
        if matches_filters(scene_name, filters):            # 通过 --scene 筛选
            scenes.append({"name": scene_name, "path": os.path.abspath(path)})  # 存绝对路径
    if not scenes:                                          # 一个都没扫到
        raise FileNotFoundError(
            "未在步骤2目录中找到 *_cal_badband_ortho_atmcorr.tif "
            "或兼容旧格式成果：%s" % in_root)
    return scenes                                           # 返回场景列表


def discover_areas(area_root, filters, logger):
    """扫描研究区矢量：每区确定唯一 .shp，校验图层、坐标系与要素。"""
    if not os.path.isdir(area_root):                        # 矢量根目录必须存在
        raise NotADirectoryError("研究区目录不存在：%s" % area_root)
    shapefiles = sorted(glob.glob(
        os.path.join(area_root, "**", "*.shp"), recursive=True))  # 递归收集全部 shp
    if not shapefiles:                                      # 没有 shp 无法继续
        raise FileNotFoundError("研究区目录下未找到 .shp 文件：%s" % area_root)

    by_folder = {}
    for shp in shapefiles:                                  # 逐个 shp 归入区域
        folder_name = os.path.basename(os.path.dirname(shp))  # 默认用上级目录名
        if os.path.abspath(os.path.dirname(shp)) == os.path.abspath(area_root):  # 直接放在根目录
            folder_name = os.path.splitext(os.path.basename(shp))[0]  # 改用文件主名
        area_name = normalize_area_name(folder_name)        # 规范化成区域名
        if not area_name or not matches_filters(area_name, filters):  # 空名或被筛掉
            continue
        by_folder.setdefault(area_name, []).append(shp)     # 收进该区域的候选表

    areas = []
    for area_name, candidates in sorted(by_folder.items()):  # 逐个区域定唯一矢量
        matching = [path for path in candidates
                    if normalize_area_name(
                        os.path.splitext(os.path.basename(path))[0]).lower()
                    == area_name.lower()]                    # 与区域同名的候选
        if len(matching) == 1:                               # 唯一同名 → 直接采用
            shp = matching[0]
            ignored = [path for path in candidates if path != shp]  # 其余作为忽略项
            if ignored:
                logger.info("区域 %s 优先采用同名矢量 %s；忽略：%s",
                            area_name, shp, ", ".join(ignored))
        elif len(matching) > 1:                              # 多个同名 → 无法判定
            raise ValueError(
                "区域 %s 有多个同名矢量，无法自动选择：%s"
                % (area_name, ", ".join(matching)))
        elif len(candidates) == 1:                           # 仅一个候选 → 采用
            shp = candidates[0]
        else:
            raise ValueError(
                "区域 %s 有多个 .shp 且没有与区域同名的文件，请整理矢量目录：%s"
                % (area_name, ", ".join(candidates)))

        dataset = ogr.Open(shp, 0)                           # 以只读方式打开矢量
        if dataset is None:
            raise RuntimeError("无法打开研究区矢量：%s" % shp)
        layer = dataset.GetLayer(0)                          # 取第一个图层
        if layer is None or layer.GetSpatialRef() is None:   # 图层与坐标系必须齐全
            dataset = None
            raise ValueError("研究区矢量缺少图层或坐标系：%s" % shp)
        if layer.GetFeatureCount() <= 0:                     # 必须含有要素
            dataset = None
            raise ValueError("研究区矢量没有要素：%s" % shp)
        dataset = None                                       # 校验完即释放数据源
        areas.append({"name": area_name, "path": os.path.abspath(shp)})  # 记录有效区域

    if not areas:                                            # 筛选后一个都没有
        raise FileNotFoundError("没有符合筛选条件的有效研究区矢量")
    return areas                                             # 返回区域列表


def open_raster(path):
    """打开输入影像并校验波段、CRS、仿射变换与各波段 NoData 一致性。"""
    dataset = gdal.Open(path, gdal.GA_ReadOnly)              # 只读打开影像
    if dataset is None:
        raise RuntimeError("无法打开输入影像：%s" % path)
    if dataset.RasterCount <= 0:                             # 波段数必须大于 0
        dataset = None
        raise ValueError("输入影像没有波段：%s" % path)
    projection = dataset.GetProjection()                     # 取 CRS 的 WKT 文本
    if not projection:
        dataset = None
        raise ValueError("输入影像缺少坐标系：%s" % path)
    transform = dataset.GetGeoTransform()                    # 取仿射变换六参数
    if abs(transform[2]) > 1e-12 or abs(transform[4]) > 1e-12:  # 含旋转项则不支持
        dataset = None
        raise ValueError("暂不支持旋转/斜切的输入栅格：%s" % path)
    if transform[1] == 0 or transform[5] == 0:               # 像元大小不能为 0
        dataset = None
        raise ValueError("输入影像像元分辨率无效：%s" % path)
    if transform[1] < 0 or transform[5] > 0:                 # 仅接受常规北向上网格
        dataset = None
        raise ValueError("仅支持常规北向上栅格（像元宽为正、像元高为负）：%s"
                         % path)

    nodata_values = [
        dataset.GetRasterBand(index).GetNoDataValue()
        for index in range(1, dataset.RasterCount + 1)
    ]                                                        # 收集各波段 NoData
    present = [value for value in nodata_values if value is not None]  # 已设置的 NoData
    if present and (len(present) != len(nodata_values)
                    or any(value != present[0] for value in present)):  # 存在不一致
        dataset = None
        raise ValueError(
            "影像波段 NoData 不一致，无法安全统一写出：%s" % path)
    nodata = present[0] if present else OUTPUT_NODATA        # 统一 NoData 取值
    return dataset, transform, projection, nodata            # 返回校验后的句柄与参数


def make_spatial_ref(wkt):
    """由 WKT 构造 OGR 空间参考，并统一为传统 GIS 轴序（经度,纬度）。"""
    srs = osr.SpatialReference()                             # 新建空空间参考
    if srs.ImportFromWkt(wkt) != 0:                          # 解析失败返回非 0
        raise ValueError("无法解析坐标系")
    if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):           # 老版本无此常量
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)  # 固定经纬度轴序
    return srs                                               # 返回空间参考对象


def intersection_pixel_bounds(dataset, transform, projection, shp):
    """求研究区与影像交集的像元窗口，返回对齐到像元边界的 bounds。

    返回 (min_x, min_y, max_x, max_y, 宽, 高)；无重叠返回 None。
    行方向按“北向上”处理：pixel_y0 = floor((y_origin - ymax)/|dy|)。
    """
    vector_ds = ogr.Open(shp, 0)                             # 只读打开研究区矢量
    if vector_ds is None:
        raise RuntimeError("无法打开研究区矢量：%s" % shp)
    layer = vector_ds.GetLayer(0)                            # 取第一个图层
    vector_srs = layer.GetSpatialRef()                       # 矢量坐标系
    if vector_srs is None:
        vector_ds = None
        raise ValueError("研究区矢量缺少坐标系：%s" % shp)
    vector_srs = vector_srs.Clone()                          # 克隆后再改，避免污染原对象
    raster_srs = make_spatial_ref(projection)                # 由影像 WKT 构造坐标系
    if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):           # 老版本兼容判断
        vector_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)  # 统一轴序
    coordinate_transform = osr.CoordinateTransformation(vector_srs, raster_srs)  # 矢量→影像

    geometries = []
    layer.ResetReading()                                     # 从头开始遍历要素
    for feature in layer:                                    # 逐个要素取几何
        geometry = feature.GetGeometryRef()                  # 取要素几何
        if geometry is None or geometry.IsEmpty():           # 空几何直接跳过
            continue
        geometry = geometry.Clone()                          # 克隆后再做坐标转换
        if geometry.Transform(coordinate_transform) != 0:    # 转换失败返回非 0
            vector_ds = None
            raise RuntimeError("研究区坐标转换失败：%s" % shp)
        geometries.append(geometry)                          # 收进影像坐标系列表
    vector_ds = None                                         # 几何取完即释放
    if not geometries:
        raise ValueError("研究区矢量没有有效几何：%s" % shp)

    area_geometry = geometries[0]                            # 以首个几何为基准
    for geometry in geometries[1:]:                          # 其余逐个合并
        area_geometry = area_geometry.Union(geometry)        # 求并集得到完整范围
        if area_geometry is None:                            # 合并失败
            raise RuntimeError("无法合并研究区几何：%s" % shp)

    x_origin, pixel_width, _, y_origin, _, pixel_height = transform  # 解包仿射参数
    raster_x1 = x_origin + dataset.RasterXSize * pixel_width         # 影像右边界 X
    raster_y1 = y_origin + dataset.RasterYSize * pixel_height        # 影像下边界 Y
    raster_ring = ogr.Geometry(ogr.wkbLinearRing)                    # 构造线性环
    for x, y in (
            (x_origin, y_origin), (raster_x1, y_origin),
            (raster_x1, raster_y1), (x_origin, raster_y1),
            (x_origin, y_origin)):                               # 四角按序闭合
        raster_ring.AddPoint_2D(x, y)                            # 逐点加入环
    raster_polygon = ogr.Geometry(ogr.wkbPolygon)                # 构造多边形
    raster_polygon.AddGeometry(raster_ring)                      # 挂上外环
    overlap = area_geometry.Intersection(raster_polygon)         # 求研究区与影像交集
    if overlap is None or overlap.IsEmpty() or overlap.GetArea() == 0:  # 无重叠
        return None

    xmin, xmax, ymin, ymax = overlap.GetEnvelope()           # 交集外接矩形
    def snap_pixel_coordinate(value):
        nearest = round(value)                               # 取最近整数
        if abs(value - nearest) <= 1e-6:                     # 已落在像元边界上
            return float(nearest)
        return value                                         # 否则保留原值

    pixel_x0 = math.floor(snap_pixel_coordinate(
        (xmin - x_origin) / pixel_width))                    # 左列号向下取整
    pixel_x1 = math.ceil(snap_pixel_coordinate(
        (xmax - x_origin) / pixel_width))                    # 右列号向上取整
    pixel_y0 = math.floor(snap_pixel_coordinate(
        (y_origin - ymax) / abs(pixel_height)))              # 上行号（北向上）
    pixel_y1 = math.ceil(snap_pixel_coordinate(
        (y_origin - ymin) / abs(pixel_height)))              # 下行号

    pixel_x0 = max(0, min(dataset.RasterXSize, pixel_x0))    # 左列限制在影像内
    pixel_x1 = max(0, min(dataset.RasterXSize, pixel_x1))    # 右列限制在影像内
    pixel_y0 = max(0, min(dataset.RasterYSize, pixel_y0))    # 上行限制在影像内
    pixel_y1 = max(0, min(dataset.RasterYSize, pixel_y1))    # 下行限制在影像内
    if pixel_x1 <= pixel_x0 or pixel_y1 <= pixel_y0:         # 窗口退化即视为无交集
        return None

    left = x_origin + pixel_x0 * pixel_width                 # 左边界地理 X
    right = x_origin + pixel_x1 * pixel_width                # 右边界地理 X
    if pixel_height < 0:                                     # 北向上时行高为负
        top = y_origin + pixel_y0 * pixel_height             # 上边界地理 Y
        bottom = y_origin + pixel_y1 * pixel_height          # 下边界地理 Y
    else:
        top = y_origin + pixel_y1 * pixel_height             # 反向时上下互换
        bottom = y_origin + pixel_y0 * pixel_height
    return (min(left, right), min(bottom, top), max(left, right),
            max(bottom, top), pixel_x1 - pixel_x0, pixel_y1 - pixel_y0)  # 窗口+宽高


def make_progress_logger(logger, scene_name, area_name):
    """构造 gdal.Warp 进度回调，每 10% 打印一次裁剪进度。"""
    state = {"last": -10}                                    # 上次已打印的进度档位

    def report(complete, message, unused_data):
        del message, unused_data                             # 回调形参未使用，直接丢弃
        percent = int(complete * 100)                        # 完成度换算成百分比
        milestone = (percent // 10) * 10                     # 归到每 10% 一档
        if milestone >= state["last"] + 10 or percent >= 100:  # 进位到下一档或已完成
            state["last"] = milestone                        # 更新已打印档位
            logger.info("  裁剪进度 [%s / %s]：%d%%",
                        scene_name, area_name, percent)
        return 1                                             # 返回 1 表示继续处理
    return report                                            # 返回进度回调


def clip_raster(scene, area, out_root, overwrite, logger):
    """裁剪一个“场景 × 区域”组合，返回 "ok" 或 "skip"。"""
    # ① 打开输入影像并校验波段/CRS/仿射变换/NoData
    dataset, transform, projection, nodata = open_raster(scene["path"])  # 读入并校验
    temporary_path = None                                    # 临时文件路径，finally 清理
    try:
        # ② 计算研究区与影像交集的像元窗口
        bounds = intersection_pixel_bounds(
            dataset, transform, projection, area["path"])    # 求交集像元窗口
        # ③ 交集为空则跳过，不生成空影像
        if bounds is None:
            logger.info("  无空间重叠，跳过：%s / %s",
                        scene["name"], area["name"])
            return "skip"                                    # 无重叠直接跳过

        area_dir = os.path.join(out_root, area["name"])      # 按区域建输出目录
        os.makedirs(area_dir, exist_ok=True)                 # 目录不存在则创建
        output_path = os.path.join(
            area_dir, scene["name"] + OUTPUT_SUFFIX)         # 正式成果完整路径
        if os.path.isfile(output_path) and not overwrite:    # 已有成果且未要求覆盖
            logger.info("  已存在，跳过：%s", output_path)
            return "skip"

        # ④ 按窗口裁剪：cutline 保留区内像元，near 重采样保持原值，不重投影
        min_x, min_y, max_x, max_y, width, height = bounds   # 解包窗口与宽高
        logger.info("  开始裁剪：%s × %s；像元窗口 %d×%d",
                    scene["name"], area["name"], width, height)
        temporary_path = output_path + ".tmp.tif"            # 先写临时文件
        if os.path.exists(temporary_path):                   # 清理历史残留临时文件
            os.remove(temporary_path)

        options = gdal.WarpOptions(
            format="GTiff",                                  # 输出格式为 GeoTIFF
            outputBounds=(min_x, min_y, max_x, max_y),       # 输出地理范围
            xRes=abs(transform[1]),                          # 沿用源 X 分辨率
            yRes=abs(transform[5]),                          # 沿用源 Y 分辨率
            resampleAlg="near",                              # 最近邻，像元原值不变
            cutlineDSName=area["path"],                      # 以研究区矢量作 cutline
            cropToCutline=False,                             # 不按 cutline 外接框裁
            srcNodata=nodata,                                # 源 NoData 值
            dstNodata=nodata,                                # 目标 NoData 值
            multithread=True,                                # 启用多线程 warp
            warpOptions=["NUM_THREADS=ALL_CPUS"],            # 使用全部 CPU 核心
            creationOptions=GTIFF_OPTIONS,                   # 复用 GeoTIFF 写选项
            outputType=gdal.GDT_Float32,                     # 输出 Float32 类型
            callback=make_progress_logger(logger, scene["name"], area["name"]),  # 进度回调
        )
        result = gdal.Warp(temporary_path, scene["path"], options=options)  # 执行裁剪
        if result is None:                                   # 裁剪失败
            raise RuntimeError("GDAL 裁剪失败：%s × %s"
                               % (scene["name"], area["name"]))
        result = None                                        # 及时释放结果数据集

        # ⑤ 回写波段 NoData、描述、元数据以及裁剪来源信息
        output = gdal.Open(temporary_path, gdal.GA_Update)   # 以更新模式重开结果
        if output is None:
            raise RuntimeError("无法重新打开裁剪结果：%s" % temporary_path)
        if output.RasterCount != dataset.RasterCount:        # 波段数须与输入一致
            output = None
            raise RuntimeError("裁剪结果波段数与输入不一致：%s" % temporary_path)
        for index in range(1, output.RasterCount + 1):       # 逐波段回写属性
            source_band = dataset.GetRasterBand(index)       # 源波段
            output_band = output.GetRasterBand(index)        # 目标波段
            output_band.SetNoDataValue(nodata)               # 统一写入 NoData
            if source_band.GetDescription():                 # 源波段有描述才搬
                output_band.SetDescription(source_band.GetDescription())  # 保留波段描述
            for key, value in source_band.GetMetadata().items():  # 逐项复制元数据项
                output_band.SetMetadataItem(key, value)
        output.SetMetadataItem("SOURCE_SCENE", scene["name"])  # 记录来源场景名
        output.SetMetadataItem("SOURCE_RASTER", os.path.basename(scene["path"]))  # 来源影像
        output.SetMetadataItem("CLIP_AREA", area["name"])    # 记录裁剪区域名
        output.SetMetadataItem("PROCESSING_STEP",
                               "Step 3: vector clip of step-2 reflectance")  # 处理步骤
        output.FlushCache()                                  # 强制刷盘
        output = None                                        # 关闭写句柄

        # ⑥ 校验尺寸/波段/像元对齐，全部通过才原子替换为正式成果
        check = gdal.Open(temporary_path, gdal.GA_ReadOnly)  # 重新打开做校验
        if (check is None or check.RasterXSize != width
                or check.RasterYSize != height
                or check.RasterCount != dataset.RasterCount):  # 尺寸或波段不符
            check = None
            raise RuntimeError("裁剪结果尺寸/波段校验失败：%s" % temporary_path)
        expected_transform = (min_x, transform[1], 0.0, max_y,
                             0.0, transform[5])              # 期望的仿射变换
        actual_transform = check.GetGeoTransform()           # 实际仿射变换
        if any(abs(actual - expected) > 1e-7
               for actual, expected in zip(actual_transform,
                                           expected_transform)):  # 未与像元网格对齐
            check = None
            raise RuntimeError("裁剪结果未与输入像元网格对齐：%s"
                               % temporary_path)
        check = None                                         # 校验完释放句柄
        os.replace(temporary_path, output_path)              # 原子替换为正式成果
        logger.info("  完成：%s", output_path)
        return "ok"                                          # 返回成功状态
    finally:
        dataset = None                                       # 释放输入数据集
        if temporary_path and os.path.exists(temporary_path):  # 有残留临时文件
            os.remove(temporary_path)                        # 失败或异常时清理


def main(argv=None):
    """批处理入口：扫描影像与矢量 → （dry-run 仅打印计划）→ 逐组合裁剪 → 汇总。"""
    args = parse_args(argv)                                  # 解析命令行参数
    log_path = args.log or os.path.join(args.out_root, "clip_batch.log")  # 默认日志路径
    logger = setup_logging(log_path, args.verbose)           # 初始化日志
    started = time.time()                                    # 总耗时计时起点
    logger.info("=" * 68)
    logger.info("ZY1E AHSI 步骤3：矢量裁剪大气校正影像")
    logger.info("步骤2输入：%s", args.in_root)
    logger.info("矢量目录：%s", args.area_root)
    logger.info("输出目录：%s（按区域直接存放，不建额外子目录）",
                args.out_root)
    logger.info("输出命名：<原始场景名>%s", OUTPUT_SUFFIX)

    scenes = discover_scenes(args.in_root, args.scene)       # 扫描待处理影像
    areas = discover_areas(args.area_root, args.area, logger)  # 扫描研究区矢量
    logger.info("待处理影像：%d 景；研究区：%d 个",
                len(scenes), len(areas))
    for area in areas:                                       # 逐个打印区域明细
        logger.info("  区域 %-12s 矢量：%s", area["name"], area["path"])

    pairs = [(scene, area) for scene in scenes for area in areas]  # 场景×区域全组合
    if args.dry_run:                                         # 仅打印计划不写盘
        for scene, area in pairs:
            output_path = os.path.join(
                args.out_root, area["name"], scene["name"] + OUTPUT_SUFFIX)  # 计划产物路径
            logger.info("计划：%s + %s -> %s",
                        scene["name"], area["name"], output_path)
        logger.info("dry-run 完成：共 %d 个场景-区域组合，未处理影像。",
                    len(pairs))
        return 0                                             # dry-run 正常退出

    counts = {"ok": 0, "skip": 0, "failed": 0}               # 三类结果计数
    for index, (scene, area) in enumerate(pairs, start=1):   # 逐组合处理
        logger.info("处理组合 %d/%d：%s × %s",
                    index, len(pairs), scene["name"], area["name"])
        pair_started = time.time()                           # 单组合计时起点
        try:
            status = clip_raster(
                scene, area, args.out_root, args.overwrite, logger)  # 执行一次裁剪
            counts[status] += 1                              # 按状态累加计数
        except Exception:                                    # 单组合失败不中断批处理
            counts["failed"] += 1
            logger.exception("处理失败：%s × %s",
                             scene["name"], area["name"])
        logger.info("组合耗时：%.1f 秒", time.time() - pair_started)

    logger.info("=" * 68)
    logger.info("处理结束：成功 %d，跳过 %d，失败 %d；总耗时 %.1f 秒",
                counts["ok"], counts["skip"], counts["failed"],
                time.time() - started)
    logger.info("裁剪成果：%s", args.out_root)
    logger.info("日志文件：%s", log_path)
    return 1 if counts["failed"] else 0                      # 有失败返回 1，否则 0


if __name__ == "__main__":
    sys.exit(main())                                         # 以 main 返回值作退出码