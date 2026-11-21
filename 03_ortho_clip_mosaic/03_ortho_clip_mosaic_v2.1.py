#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Clip step-2 reflectance rasters by research-area shapefiles (v2.1).

按研究区矢量裁剪步骤2 大气校正影像：不重投影、不重采样，保留源 CRS、
像元大小、波段顺序与像元对齐。流程：扫描影像/矢量 → 求交集像元窗口
→ gdal.Warp 按 cutline 裁剪 → 回写波段元数据 → 校验网格后原子替换。

输出直接写到 <output-root>/<area-name>/<scene-name>_cal_badband_ortho_atmcorr_clip.tif，
不建立中间区域子目录。

Examples:
    python 03_clip_by_area_v2.1.py --dry-run
    python 03_clip_by_area_v2.1.py
    python 03_clip_by_area_v2.1.py --area JLHK --scene "*20250121*"
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

gdal.UseExceptions()
ogr.UseExceptions()

IN_ROOT_STEP2 = r"E:\ZhouYang\test\ZY1E_output_step2"
OUT_ROOT_STEP3 = r"E:\ZhouYang\test\ZY1E_output_step3"
AREA_ROOT = r"E:\ZhouYang\test\Hyperarea2"

IN_REFLECTANCE_SUFFIX = "_cal_badband_ortho_atmcorr.tif"
LEGACY_REFLECTANCE_SUFFIXES = (
    "_ortho_atmcorr.tif",
    "reflectance_atm.tif",
)
OUTPUT_SUFFIX = "_cal_badband_ortho_atmcorr_clip.tif"
OUTPUT_NODATA = -9999.0
LOG_NAME = "zy1e.step3.clip"
GTIFF_OPTIONS = [
    "BIGTIFF=IF_SAFER",
    "TILED=YES",
    "INTERLEAVE=BAND",
    "COMPRESS=DEFLATE",
    "PREDICTOR=3",
    "ZLEVEL=6",
]


def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(
        description="按研究区矢量裁剪步骤2大气校正影像（不重投影、不重采样）")
    parser.add_argument("--in", dest="in_root", default=IN_ROOT_STEP2,
                        help="步骤2成果根目录（默认：%(default)s）")
    parser.add_argument("--out", dest="out_root", default=OUT_ROOT_STEP3,
                        help="裁剪成果根目录（默认：%(default)s）")
    parser.add_argument("--areas", dest="area_root", default=AREA_ROOT,
                        help="研究区矢量根目录（默认：%(default)s）")
    parser.add_argument("--area", action="append", default=None,
                        help="只处理匹配的区域名，可多次指定，支持通配符")
    parser.add_argument("--scene", action="append", default=None,
                        help="只处理匹配的场景名，可多次指定，支持通配符")
    parser.add_argument("--overwrite", action="store_true",
                        help="覆盖已有裁剪成果（默认跳过）")
    parser.add_argument("--dry-run", action="store_true",
                        help="只显示处理计划，不写出影像")
    parser.add_argument("--log", default=None,
                        help="日志文件（默认：<输出根目录>\\clip_batch.log）")
    parser.add_argument("--verbose", action="store_true",
                        help="输出调试日志")
    return parser.parse_args(argv)


def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(logging.DEBUG)
    for existing_handler in list(logger.handlers):
        logger.removeHandler(existing_handler)
        existing_handler.close()
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S")
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(formatter)
    logger.addHandler(console)
    if log_path:
        parent = os.path.dirname(os.path.abspath(log_path))
        if parent:
            os.makedirs(parent, exist_ok=True)
        handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
        handler.setLevel(logging.DEBUG)
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    logger.propagate = False
    return logger


def normalize_area_name(name):
    """规范化区域名：去序号前缀、空白与非法字符，作为输出目录名。"""
    value = os.path.basename(os.path.normpath(str(name)))
    value = re.sub(r"^\s*\d+\s*[\.\-_]?\s*", "", value)
    value = re.sub(r"\s+", "", value)
    value = re.sub(r'[<>:"/\\|?*]', "_", value)
    return value.strip(" ._")


def matches_filters(name, filters):
    """名称是否匹配任一通配符筛选条件；filters 为空表示全部匹配。"""
    if not filters:
        return True
    return any(fnmatch.fnmatchcase(name.lower(), pattern.lower())
               for pattern in filters)


def discover_scenes(in_root, filters):
    """扫描步骤2 成果（新后缀优先、兼容旧后缀），返回场景列表。"""
    if not os.path.isdir(in_root):
        raise NotADirectoryError("步骤2成果目录不存在：%s" % in_root)
    paths_by_scene = {}
    for path in sorted(glob.glob(
            os.path.join(in_root, "**", "*" + IN_REFLECTANCE_SUFFIX),
            recursive=True)):
        scene_dir = os.path.dirname(path)
        paths_by_scene[scene_dir] = path
    for suffix in LEGACY_REFLECTANCE_SUFFIXES:
        legacy_name = suffix if suffix.endswith(".tif") else suffix + ".tif"
        pattern_name = ("*" + suffix if suffix.startswith("_") else legacy_name)
        for path in sorted(glob.glob(
                os.path.join(in_root, "**", pattern_name), recursive=True)):
            paths_by_scene.setdefault(os.path.dirname(path), path)

    scenes = []
    for scene_dir, path in sorted(paths_by_scene.items()):
        scene_name = os.path.basename(scene_dir)
        if matches_filters(scene_name, filters):
            scenes.append({"name": scene_name, "path": os.path.abspath(path)})
    if not scenes:
        raise FileNotFoundError(
            "未在步骤2目录中找到 *_cal_badband_ortho_atmcorr.tif "
            "或兼容旧格式成果：%s" % in_root)
    return scenes


def discover_areas(area_root, filters, logger):
    """扫描研究区矢量：每区确定唯一 .shp，校验图层、坐标系与要素。"""
    if not os.path.isdir(area_root):
        raise NotADirectoryError("研究区目录不存在：%s" % area_root)
    shapefiles = sorted(glob.glob(
        os.path.join(area_root, "**", "*.shp"), recursive=True))
    if not shapefiles:
        raise FileNotFoundError("研究区目录下未找到 .shp 文件：%s" % area_root)

    by_folder = {}
    for shp in shapefiles:
        folder_name = os.path.basename(os.path.dirname(shp))
        if os.path.abspath(os.path.dirname(shp)) == os.path.abspath(area_root):
            folder_name = os.path.splitext(os.path.basename(shp))[0]
        area_name = normalize_area_name(folder_name)
        if not area_name or not matches_filters(area_name, filters):
            continue
        by_folder.setdefault(area_name, []).append(shp)

    areas = []
    for area_name, candidates in sorted(by_folder.items()):
        matching = [path for path in candidates
                    if normalize_area_name(
                        os.path.splitext(os.path.basename(path))[0]).lower()
                    == area_name.lower()]
        if len(matching) == 1:
            shp = matching[0]
            ignored = [path for path in candidates if path != shp]
            if ignored:
                logger.info("区域 %s 优先采用同名矢量 %s；忽略：%s",
                            area_name, shp, ", ".join(ignored))
        elif len(matching) > 1:
            raise ValueError(
                "区域 %s 有多个同名矢量，无法自动选择：%s"
                % (area_name, ", ".join(matching)))
        elif len(candidates) == 1:
            shp = candidates[0]
        else:
            raise ValueError(
                "区域 %s 有多个 .shp 且没有与区域同名的文件，请整理矢量目录：%s"
                % (area_name, ", ".join(candidates)))

        dataset = ogr.Open(shp, 0)
        if dataset is None:
            raise RuntimeError("无法打开研究区矢量：%s" % shp)
        layer = dataset.GetLayer(0)
        if layer is None or layer.GetSpatialRef() is None:
            dataset = None
            raise ValueError("研究区矢量缺少图层或坐标系：%s" % shp)
        if layer.GetFeatureCount() <= 0:
            dataset = None
            raise ValueError("研究区矢量没有要素：%s" % shp)
        dataset = None
        areas.append({"name": area_name, "path": os.path.abspath(shp)})

    if not areas:
        raise FileNotFoundError("没有符合筛选条件的有效研究区矢量")
    return areas


def open_raster(path):
    """打开输入影像并校验波段、CRS、仿射变换与各波段 NoData 一致性。"""
    dataset = gdal.Open(path, gdal.GA_ReadOnly)
    if dataset is None:
        raise RuntimeError("无法打开输入影像：%s" % path)
    if dataset.RasterCount <= 0:
        dataset = None
        raise ValueError("输入影像没有波段：%s" % path)
    projection = dataset.GetProjection()
    if not projection:
        dataset = None
        raise ValueError("输入影像缺少坐标系：%s" % path)
    transform = dataset.GetGeoTransform()
    if abs(transform[2]) > 1e-12 or abs(transform[4]) > 1e-12:
        dataset = None
        raise ValueError("暂不支持旋转/斜切的输入栅格：%s" % path)
    if transform[1] == 0 or transform[5] == 0:
        dataset = None
        raise ValueError("输入影像像元分辨率无效：%s" % path)
    if transform[1] < 0 or transform[5] > 0:
        dataset = None
        raise ValueError("仅支持常规北向上栅格（像元宽为正、像元高为负）：%s"
                         % path)

    nodata_values = [
        dataset.GetRasterBand(index).GetNoDataValue()
        for index in range(1, dataset.RasterCount + 1)
    ]
    present = [value for value in nodata_values if value is not None]
    if present and (len(present) != len(nodata_values)
                    or any(value != present[0] for value in present)):
        dataset = None
        raise ValueError(
            "影像波段 NoData 不一致，无法安全统一写出：%s" % path)
    nodata = present[0] if present else OUTPUT_NODATA
    return dataset, transform, projection, nodata


def make_spatial_ref(wkt):
    """由 WKT 构造 OGR 空间参考，并统一为传统 GIS 轴序（经度,纬度）。"""
    srs = osr.SpatialReference()
    if srs.ImportFromWkt(wkt) != 0:
        raise ValueError("无法解析坐标系")
    if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):
        srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    return srs


def intersection_pixel_bounds(dataset, transform, projection, shp):
    """求研究区与影像交集的像元窗口，返回对齐到像元边界的 bounds。

    返回 (min_x, min_y, max_x, max_y, 宽, 高)；无重叠返回 None。
    行方向按“北向上”处理：pixel_y0 = floor((y_origin - ymax)/|dy|)。
    """
    vector_ds = ogr.Open(shp, 0)
    if vector_ds is None:
        raise RuntimeError("无法打开研究区矢量：%s" % shp)
    layer = vector_ds.GetLayer(0)
    vector_srs = layer.GetSpatialRef()
    if vector_srs is None:
        vector_ds = None
        raise ValueError("研究区矢量缺少坐标系：%s" % shp)
    vector_srs = vector_srs.Clone()
    raster_srs = make_spatial_ref(projection)
    if hasattr(osr, "OAMS_TRADITIONAL_GIS_ORDER"):
        vector_srs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    coordinate_transform = osr.CoordinateTransformation(vector_srs, raster_srs)

    geometries = []
    layer.ResetReading()
    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry is None or geometry.IsEmpty():
            continue
        geometry = geometry.Clone()
        if geometry.Transform(coordinate_transform) != 0:
            vector_ds = None
            raise RuntimeError("研究区坐标转换失败：%s" % shp)
        geometries.append(geometry)
    vector_ds = None
    if not geometries:
        raise ValueError("研究区矢量没有有效几何：%s" % shp)

    area_geometry = geometries[0]
    for geometry in geometries[1:]:
        area_geometry = area_geometry.Union(geometry)
        if area_geometry is None:
            raise RuntimeError("无法合并研究区几何：%s" % shp)

    x_origin, pixel_width, _, y_origin, _, pixel_height = transform
    raster_x1 = x_origin + dataset.RasterXSize * pixel_width
    raster_y1 = y_origin + dataset.RasterYSize * pixel_height
    raster_ring = ogr.Geometry(ogr.wkbLinearRing)
    for x, y in (
            (x_origin, y_origin), (raster_x1, y_origin),
            (raster_x1, raster_y1), (x_origin, raster_y1),
            (x_origin, y_origin)):
        raster_ring.AddPoint_2D(x, y)
    raster_polygon = ogr.Geometry(ogr.wkbPolygon)
    raster_polygon.AddGeometry(raster_ring)
    overlap = area_geometry.Intersection(raster_polygon)
    if overlap is None or overlap.IsEmpty() or overlap.GetArea() == 0:
        return None

    xmin, xmax, ymin, ymax = overlap.GetEnvelope()
    def snap_pixel_coordinate(value):
        nearest = round(value)
        if abs(value - nearest) <= 1e-6:
            return float(nearest)
        return value

    pixel_x0 = math.floor(snap_pixel_coordinate(
        (xmin - x_origin) / pixel_width))
    pixel_x1 = math.ceil(snap_pixel_coordinate(
        (xmax - x_origin) / pixel_width))
    pixel_y0 = math.floor(snap_pixel_coordinate(
        (y_origin - ymax) / abs(pixel_height)))
    pixel_y1 = math.ceil(snap_pixel_coordinate(
        (y_origin - ymin) / abs(pixel_height)))

    pixel_x0 = max(0, min(dataset.RasterXSize, pixel_x0))
    pixel_x1 = max(0, min(dataset.RasterXSize, pixel_x1))
    pixel_y0 = max(0, min(dataset.RasterYSize, pixel_y0))
    pixel_y1 = max(0, min(dataset.RasterYSize, pixel_y1))
    if pixel_x1 <= pixel_x0 or pixel_y1 <= pixel_y0:
        return None

    left = x_origin + pixel_x0 * pixel_width
    right = x_origin + pixel_x1 * pixel_width
    if pixel_height < 0:
        top = y_origin + pixel_y0 * pixel_height
        bottom = y_origin + pixel_y1 * pixel_height
    else:
        top = y_origin + pixel_y1 * pixel_height
        bottom = y_origin + pixel_y0 * pixel_height
    return (min(left, right), min(bottom, top), max(left, right),
            max(bottom, top), pixel_x1 - pixel_x0, pixel_y1 - pixel_y0)


def make_progress_logger(logger, scene_name, area_name):
    """构造 gdal.Warp 进度回调，每 10% 打印一次裁剪进度。"""
    state = {"last": -10}

    def report(complete, message, unused_data):
        del message, unused_data
        percent = int(complete * 100)
        milestone = (percent // 10) * 10
        if milestone >= state["last"] + 10 or percent >= 100:
            state["last"] = milestone
            logger.info("  裁剪进度 [%s / %s]：%d%%",
                        scene_name, area_name, percent)
        return 1
    return report


def clip_raster(scene, area, out_root, overwrite, logger):
    """裁剪一个“场景 × 区域”组合，返回 "ok" 或 "skip"。"""
    # ① 打开输入影像并校验波段/CRS/仿射变换/NoData
    dataset, transform, projection, nodata = open_raster(scene["path"])
    temporary_path = None
    try:
        # ② 计算研究区与影像交集的像元窗口
        bounds = intersection_pixel_bounds(
            dataset, transform, projection, area["path"])
        # ③ 交集为空则跳过，不生成空影像
        if bounds is None:
            logger.info("  无空间重叠，跳过：%s / %s",
                        scene["name"], area["name"])
            return "skip"

        area_dir = os.path.join(out_root, area["name"])
        os.makedirs(area_dir, exist_ok=True)
        output_path = os.path.join(
            area_dir, scene["name"] + OUTPUT_SUFFIX)
        if os.path.isfile(output_path) and not overwrite:
            logger.info("  已存在，跳过：%s", output_path)
            return "skip"

        # ④ 按窗口裁剪：cutline 保留区内像元，near 重采样保持原值，不重投影
        min_x, min_y, max_x, max_y, width, height = bounds
        logger.info("  开始裁剪：%s × %s；像元窗口 %d×%d",
                    scene["name"], area["name"], width, height)
        temporary_path = output_path + ".tmp.tif"
        if os.path.exists(temporary_path):
            os.remove(temporary_path)

        options = gdal.WarpOptions(
            format="GTiff",
            outputBounds=(min_x, min_y, max_x, max_y),
            xRes=abs(transform[1]),
            yRes=abs(transform[5]),
            resampleAlg="near",
            cutlineDSName=area["path"],
            cropToCutline=False,
            srcNodata=nodata,
            dstNodata=nodata,
            multithread=True,
            warpOptions=["NUM_THREADS=ALL_CPUS"],
            creationOptions=GTIFF_OPTIONS,
            outputType=gdal.GDT_Float32,
            callback=make_progress_logger(logger, scene["name"], area["name"]),
        )
        result = gdal.Warp(temporary_path, scene["path"], options=options)
        if result is None:
            raise RuntimeError("GDAL 裁剪失败：%s × %s"
                               % (scene["name"], area["name"]))
        result = None

        # ⑤ 回写波段 NoData、描述、元数据以及裁剪来源信息
        output = gdal.Open(temporary_path, gdal.GA_Update)
        if output is None:
            raise RuntimeError("无法重新打开裁剪结果：%s" % temporary_path)
        if output.RasterCount != dataset.RasterCount:
            output = None
            raise RuntimeError("裁剪结果波段数与输入不一致：%s" % temporary_path)
        for index in range(1, output.RasterCount + 1):
            source_band = dataset.GetRasterBand(index)
            output_band = output.GetRasterBand(index)
            output_band.SetNoDataValue(nodata)
            if source_band.GetDescription():
                output_band.SetDescription(source_band.GetDescription())
            for key, value in source_band.GetMetadata().items():
                output_band.SetMetadataItem(key, value)
        output.SetMetadataItem("SOURCE_SCENE", scene["name"])
        output.SetMetadataItem("SOURCE_RASTER", os.path.basename(scene["path"]))
        output.SetMetadataItem("CLIP_AREA", area["name"])
        output.SetMetadataItem("PROCESSING_STEP",
                               "Step 3: vector clip of step-2 reflectance")
        output.FlushCache()
        output = None

        # ⑥ 校验尺寸/波段/像元对齐，全部通过才原子替换为正式成果
        check = gdal.Open(temporary_path, gdal.GA_ReadOnly)
        if (check is None or check.RasterXSize != width
                or check.RasterYSize != height
                or check.RasterCount != dataset.RasterCount):
            check = None
            raise RuntimeError("裁剪结果尺寸/波段校验失败：%s" % temporary_path)
        expected_transform = (min_x, transform[1], 0.0, max_y,
                             0.0, transform[5])
        actual_transform = check.GetGeoTransform()
        if any(abs(actual - expected) > 1e-7
               for actual, expected in zip(actual_transform,
                                           expected_transform)):
            check = None
            raise RuntimeError("裁剪结果未与输入像元网格对齐：%s"
                               % temporary_path)
        check = None
        os.replace(temporary_path, output_path)
        logger.info("  完成：%s", output_path)
        return "ok"
    finally:
        dataset = None
        if temporary_path and os.path.exists(temporary_path):
            os.remove(temporary_path)


def main(argv=None):
    """批处理入口：扫描影像与矢量 → （dry-run 仅打印计划）→ 逐组合裁剪 → 汇总。"""
    args = parse_args(argv)
    log_path = args.log or os.path.join(args.out_root, "clip_batch.log")
    logger = setup_logging(log_path, args.verbose)
    started = time.time()
    logger.info("=" * 68)
    logger.info("ZY1E AHSI 步骤3：矢量裁剪大气校正影像")
    logger.info("步骤2输入：%s", args.in_root)
    logger.info("矢量目录：%s", args.area_root)
    logger.info("输出目录：%s（按区域直接存放，不建额外子目录）",
                args.out_root)
    logger.info("输出命名：<原始场景名>%s", OUTPUT_SUFFIX)

    scenes = discover_scenes(args.in_root, args.scene)
    areas = discover_areas(args.area_root, args.area, logger)
    logger.info("待处理影像：%d 景；研究区：%d 个",
                len(scenes), len(areas))
    for area in areas:
        logger.info("  区域 %-12s 矢量：%s", area["name"], area["path"])

    pairs = [(scene, area) for scene in scenes for area in areas]
    if args.dry_run:
        for scene, area in pairs:
            output_path = os.path.join(
                args.out_root, area["name"], scene["name"] + OUTPUT_SUFFIX)
            logger.info("计划：%s + %s -> %s",
                        scene["name"], area["name"], output_path)
        logger.info("dry-run 完成：共 %d 个场景-区域组合，未处理影像。",
                    len(pairs))
        return 0

    counts = {"ok": 0, "skip": 0, "failed": 0}
    for index, (scene, area) in enumerate(pairs, start=1):
        logger.info("处理组合 %d/%d：%s × %s",
                    index, len(pairs), scene["name"], area["name"])
        pair_started = time.time()
        try:
            status = clip_raster(
                scene, area, args.out_root, args.overwrite, logger)
            counts[status] += 1
        except Exception:
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
    return 1 if counts["failed"] else 0


if __name__ == "__main__":
    sys.exit(main())