# -*- coding: utf-8 -*-
"""
00_extract_georeferenced_jpg.py
提取快视图 JPG，将 SHP 边界绘制进 JPG，并为图片赋予地理参考（v1.0.1）

处理顺序：定位 JPG/SHP → 提取/复制 → 读 SHP 范围与投影 → 写 .jgw/.prj
→ 绘制边界 → 输出带坐标的 GeoTIFF。
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

gdal.UseExceptions()
gdal.PushErrorHandler("CPLQuietErrorHandler")
ogr.UseExceptions()

# =================================================================
# 配置区（按需修改）
# =================================================================
SRC_DIR = r"E:\ZhouYang\test\ZY1E_AHSI"
OUT_ROOT = r"E:\ZhouYang\test\ZY1E_Georef_JPG"

# 支持的压缩包格式
TAR_PATTERNS = ("*.tar", "*.tar.gz", "*.tgz")

# 输出设置
OVERWRITE = False
LOG_NAME = "zy1e.georef_jpg_shp"
BOUNDARY_COLOR = (255, 0, 0)
BOUNDARY_WIDTH = 2
JPEG_QUALITY = 95
GEOTIFF_JPEG_QUALITY = 90
BOUNDARY_RENDER_MARKER = ".boundary_rendered"

# =================================================================
# 通用工具
# =================================================================
def parse_args(argv=None):
    """解析命令行参数；argv 为 None 时取 sys.argv[1:]。"""
    parser = argparse.ArgumentParser(description="ZY1E 利用SHP提取并赋予JPG坐标")
    parser.add_argument("--src", default=SRC_DIR)
    parser.add_argument("--out", default=OUT_ROOT)
    parser.add_argument("--pattern", action="append", default=None)
    parser.add_argument("--scene", action="append", default=None)
    parser.add_argument("--overwrite", action="store_true", default=OVERWRITE)
    parser.add_argument("--limit", type=int, default=0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--log", default=None)
    parser.add_argument("--verbose", action="store_true")
    return parser.parse_args(argv)

def setup_logging(log_path, verbose=False):
    """配置日志：控制台（INFO/verbose=DEBUG）+ 文件（DEBUG，追加写）。"""
    logger = logging.getLogger(LOG_NAME)
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    
    console = logging.StreamHandler(stream=sys.stdout)
    console.setLevel(logging.DEBUG if verbose else logging.INFO)
    console.setFormatter(fmt)
    logger.addHandler(console)
    
    if log_path:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(log_path)), exist_ok=True)
            file_handler = logging.FileHandler(log_path, mode="a", encoding="utf-8")
            file_handler.setLevel(logging.DEBUG)
            file_handler.setFormatter(fmt)
            logger.addHandler(file_handler)
        except OSError:
            logger.exception("无法创建日志文件: %s", log_path)
            raise
    logger.propagate = False
    return logger

def discover_scenes(src_dir, patterns):
    """扫描源目录：返回压缩包路径，以及同时含 .jpg 和 .shp 的目录。"""
    if not os.path.isdir(src_dir):
        raise NotADirectoryError("源目录不存在: %s" % src_dir)
    found = set()
    for pattern in patterns:
        found.update(glob.glob(os.path.join(src_dir, pattern), recursive=True))
        found.update(glob.glob(os.path.join(src_dir, "**", pattern), recursive=True))
    for current_dir, directories, filenames in os.walk(src_dir):
        directories[:] = [name for name in directories if not name.startswith(".")]
        extensions = {os.path.splitext(name)[1].lower() for name in filenames}
        if ".jpg" in extensions and ".shp" in extensions:
            found.add(current_dir)
    return sorted(f for f in found if os.path.isfile(f) or os.path.isdir(f))

def get_scene_name(path):
    """由压缩包或目录路径推导景名（去掉 .tar/.tar.gz/.tgz 等后缀）。"""
    name = os.path.basename(path)
    if os.path.isdir(path):
        return name
    lower = name.lower()
    for ext in (".tar.gz", ".tgz", ".tar"):
        if lower.endswith(ext): return name[: -len(ext)]
    return os.path.splitext(name)[0]

def extract_tar_member(tar, member_name, out_path):
    """从已打开的 tar 包中流式提取单个成员到 out_path。"""
    source = tar.extractfile(member_name)
    if source is None:
        raise OSError("无法从 TAR 包读取文件: %s" % member_name)
    with source, open(out_path, "wb") as destination:
        shutil.copyfileobj(source, destination)


def matching_sidecars(files, shp_name):
    """返回与 .shp 同名、同目录的标准 sidecar（.dbf/.shx/.prj/.cpg）。"""
    shp_stem = os.path.splitext(os.path.basename(shp_name))[0].casefold()
    shp_dir = os.path.dirname(shp_name).casefold()
    sidecar_extensions = {".dbf", ".shx", ".prj", ".cpg"}
    return [
        path for path in files
        if os.path.dirname(path).casefold() == shp_dir
        and os.path.splitext(os.path.basename(path))[0].casefold() == shp_stem
        and os.path.splitext(path)[1].lower() in sidecar_extensions
    ]


def burn_shapefile_boundary(jpg_path, shp_path, extent, logger):
    """把 SHP 的线/面边界按 extent 线性映射到 JPG 像素坐标，用红线绘制。

    像素坐标换算（extent 为 SHP 外接矩形 minX/maxX/minY/maxY）：
        px = (x - min_x) / xres
        py = (max_y - y) / yres      # 图像 y 轴向下，故用 max_y 减
    结果先写临时 JPG 再 os.replace 覆盖原图，保证原子性。
    """
    min_x, max_x, min_y, max_y = extent
    image = Image.open(jpg_path).convert("RGB")
    width, height = image.size
    xres = (max_x - min_x) / width
    yres = (max_y - min_y) / height
    draw = ImageDraw.Draw(image)

    vector = ogr.Open(shp_path, 0)
    if vector is None:
        raise RuntimeError("绘制边界时无法重新打开 SHP: %s" % shp_path)
    layer = vector.GetLayer()
    if layer is None:
        vector = None
        raise RuntimeError("SHP 数据源中没有可绘制的图层")

    feature_count = 0
    try:
        for feature in layer:
            geometry = feature.GetGeometryRef()
            if geometry is None or geometry.IsEmpty():
                continue
            boundary = geometry.Boundary() if geometry.GetGeometryName().upper() in (
                "POLYGON", "MULTIPOLYGON"
            ) else geometry.Clone()
            if boundary is None or boundary.IsEmpty():
                continue

            parts = [boundary]
            while parts:
                part = parts.pop()
                if part.GetGeometryCount():
                    parts.extend(
                        part.GetGeometryRef(index)
                        for index in range(part.GetGeometryCount())
                    )
                    continue
                if part.GetPointCount() < 2:
                    continue
                pixels = [
                    (
                        (part.GetX(index) - min_x) / xres,
                        (max_y - part.GetY(index)) / yres,
                    )
                    for index in range(part.GetPointCount())
                ]
                draw.line(
                    pixels,
                    fill=BOUNDARY_COLOR,
                    width=BOUNDARY_WIDTH,
                    joint="curve",
                )
                feature_count += 1
    finally:
        vector = None

    if feature_count == 0:
        raise ValueError("SHP 中没有可绘制的线或面边界")

    temp_path = jpg_path + ".boundary.tmp.jpg"
    try:
        image.save(temp_path, format="JPEG", quality=JPEG_QUALITY)
        os.replace(temp_path, jpg_path)
    except Exception:
        if os.path.exists(temp_path):
            os.remove(temp_path)
        raise
    logger.info(
        "  - 已将 SHP 边界绘制到 JPG（红色，%d 像素宽）",
        BOUNDARY_WIDTH,
    )


def create_georeferenced_geotiff(jpg_path, tif_path, extent, projection_wkt,
                                 logger):
    """给已绘制边界的 JPG 写入仿射变换与投影，输出 JPEG 压缩的 GeoTIFF。

    geotransform = (min_x, 像元宽, 0, max_y, 0, -像元高)
    其中 像元宽=(max_x-min_x)/宽，像元高=(max_y-min_y)/高；写出后重新打开
    校验波段数、尺寸、投影与仿射变换，全部通过才原子替换为正式文件。
    """
    min_x, max_x, min_y, max_y = extent
    source = gdal.Open(jpg_path, gdal.GA_ReadOnly)
    if source is None:
        raise RuntimeError("无法打开已绘制边界的 RGB 预览图: %s" % jpg_path)
    if source.RasterCount != 3:
        source = None
        raise ValueError("RGB 预览图应为 3 个波段，实际波段数不符")

    width, height = source.RasterXSize, source.RasterYSize
    geotransform = (
        min_x,
        (max_x - min_x) / width,
        0.0,
        max_y,
        0.0,
        -(max_y - min_y) / height,
    )
    source.SetGeoTransform(geotransform)
    source.SetProjection(projection_wkt)
    temp_path = tif_path + ".tmp.tif"
    driver = gdal.GetDriverByName("GTiff")
    if driver is None:
        source = None
        raise RuntimeError("GDAL 未提供 GeoTIFF 驱动")

    try:
        output = driver.CreateCopy(
            temp_path,
            source,
            strict=1,
            options=[
                "COMPRESS=JPEG",
                "JPEG_QUALITY=%d" % GEOTIFF_JPEG_QUALITY,
                "PHOTOMETRIC=YCBCR",
                "TILED=YES",
            ],
        )
        if output is None:
            raise RuntimeError("创建压缩 GeoTIFF 失败: %s" % temp_path)
        output.FlushCache()
        output = None
        source = None

        check = gdal.Open(temp_path, gdal.GA_ReadOnly)
        if check is None:
            raise RuntimeError("生成的 GeoTIFF 无法重新打开")
        try:
            if (check.RasterCount != 3
                    or check.RasterXSize != width
                    or check.RasterYSize != height):
                raise RuntimeError("GeoTIFF 波段数或图像尺寸与预览 JPG 不一致")
            if not check.GetProjection():
                raise RuntimeError("GeoTIFF 未写入坐标参考系")
            if check.GetGeoTransform() != geotransform:
                raise RuntimeError("GeoTIFF 仿射变换与预期不一致")
        finally:
            check = None
        os.replace(temp_path, tif_path)
    except Exception:
        source = None
        if os.path.exists(temp_path):
            try:
                os.remove(temp_path)
            except OSError:
                logger.exception("清理临时 GeoTIFF 失败: %s", temp_path)
        raise

    logger.info(
        "  - 已生成带坐标系的 RGB GeoTIFF: %s (JPEG 压缩)",
        os.path.basename(tif_path),
    )

# =================================================================
# 单景处理逻辑
# =================================================================
def process_scene(scene_path, out_root, args, logger):
    """处理单景：①定位文件 ②提取/复制 ③读 SHP 范围 ④算分辨率

    ⑤写 .jgw ⑥写 .prj ⑦绘边界并输出 GeoTIFF。
    单景异常在此捕获并记录堆栈，返回状态字典，不中断批处理。
    """
    scene = get_scene_name(scene_path)
    scene_dir = os.path.join(out_root, scene)

    logger.info("开始处理场景: %s", scene)
    started = time.time()
    stage_started = started

    try:
        os.makedirs(scene_dir, exist_ok=True)
        jpg_member = None
        shp_member = None
        shp_related = []
        is_tar = not os.path.isdir(scene_path)
        
        out_jpg_path = None
        out_shp_path = None

        # ① 定位文件：压缩包内按后缀查找 JPG/SHP；目录则直接列文件
        if is_tar:
            with tarfile.open(scene_path, "r:*") as tar:
                members = [m for m in tar.getmembers() if m.isfile()]
                jpg_members = sorted(m.name for m in members
                                     if m.name.lower().endswith(".jpg"))
                shp_members = sorted(m.name for m in members
                                     if m.name.lower().endswith(".shp"))
                if not jpg_members:
                    raise FileNotFoundError("包内未找到 JPG")
                if not shp_members:
                    raise FileNotFoundError("包内未找到 SHP 矢量边界")
                if len(jpg_members) > 1 or len(shp_members) > 1:
                    logger.warning("包内 JPG/SHP 不止一个，按排序选择首个文件。")
                jpg_member = jpg_members[0]
                shp_member = shp_members[0]
                shp_dir = os.path.dirname(shp_member)
                shp_related = matching_sidecars(
                    [m.name for m in members
                     if os.path.dirname(m.name) == shp_dir],
                    shp_member,
                )

                out_jpg_path = os.path.join(scene_dir, os.path.basename(jpg_member))
                out_shp_path = os.path.join(scene_dir, os.path.basename(shp_member))
                out_tif_path = os.path.splitext(out_jpg_path)[0] + ".tif"
                jgw_path = os.path.splitext(out_jpg_path)[0] + ".jgw"
                jpg_prj_path = os.path.splitext(out_jpg_path)[0] + ".prj"
                expected_outputs = [
                    out_jpg_path, out_shp_path, out_tif_path, jgw_path, jpg_prj_path,
                    os.path.join(scene_dir, BOUNDARY_RENDER_MARKER),
                    *(os.path.join(scene_dir, os.path.basename(path))
                      for path in shp_related),
                ]
                if not args.overwrite and all(os.path.isfile(path)
                                              for path in expected_outputs):
                    logger.info("跳过（成果已存在）: %s", scene)
                    return {"scene": scene, "status": "跳过", "elapsed": 0.0}

                logger.info("  - 正在提取 JPG 和 SHP 矢量...")
                extract_tar_member(tar, jpg_member, out_jpg_path)
                extract_tar_member(tar, shp_member, out_shp_path)
                for rel in shp_related:
                    extract_tar_member(tar, rel, os.path.join(scene_dir, os.path.basename(rel)))
        else:
            source_files = [
                os.path.join(scene_path, name) for name in os.listdir(scene_path)
                if os.path.isfile(os.path.join(scene_path, name))
            ]
            jpg_files = sorted(path for path in source_files
                               if path.lower().endswith(".jpg"))
            shp_files = sorted(path for path in source_files
                               if path.lower().endswith(".shp"))
            if not jpg_files or not shp_files:
                raise FileNotFoundError("文件夹内缺少 JPG 或 SHP")
            if len(jpg_files) > 1 or len(shp_files) > 1:
                logger.warning("场景目录内 JPG/SHP 不止一个，按排序选择首个文件。")

            out_jpg_path = os.path.join(scene_dir, os.path.basename(jpg_files[0]))
            out_shp_path = os.path.join(scene_dir, os.path.basename(shp_files[0]))
            out_tif_path = os.path.splitext(out_jpg_path)[0] + ".tif"
            shp_related = matching_sidecars(source_files, shp_files[0])
            jgw_path = os.path.splitext(out_jpg_path)[0] + ".jgw"
            jpg_prj_path = os.path.splitext(out_jpg_path)[0] + ".prj"
            expected_outputs = [
                out_jpg_path, out_shp_path, out_tif_path, jgw_path, jpg_prj_path,
                os.path.join(scene_dir, BOUNDARY_RENDER_MARKER),
                *(os.path.join(scene_dir, os.path.basename(path))
                  for path in shp_related),
            ]
            if not args.overwrite and all(os.path.isfile(path)
                                          for path in expected_outputs):
                logger.info("跳过（成果已存在）: %s", scene)
                return {"scene": scene, "status": "跳过", "elapsed": 0.0}

            shutil.copy2(jpg_files[0], out_jpg_path)
            shutil.copy2(shp_files[0], out_shp_path)
            for sidecar in shp_related:
                shutil.copy2(sidecar, scene_dir)

        logger.info("  - 文件准备完成，用时 %.1f s", time.time() - stage_started)
        stage_started = time.time()

        # ② 读取 SHP 获取坐标范围和投影
        logger.info("  - 正在解析 SHP 获取空间坐标...")
        shp_ds = ogr.Open(out_shp_path, 0)
        if shp_ds is None: raise RuntimeError(f"无法打开矢量: {out_shp_path}")
        
        layer = shp_ds.GetLayer()
        if layer is None:
            raise RuntimeError("SHP 数据源中没有可读取的图层")
        extent = layer.GetExtent(force=1) # (minX, maxX, minY, maxY)
        spatial_ref = layer.GetSpatialRef()
        if extent is None or not all(math.isfinite(value) for value in extent):
            raise ValueError("SHP 的空间范围为空或包含无效坐标")
        
        if spatial_ref:
            # 转换为 WKT 以写入 .prj
            prj_wkt = spatial_ref.ExportToWkt()
            # 检查是不是地理坐标 (经纬度)
            is_geo = spatial_ref.IsGeographic()
            unit_name = "度" if is_geo else "米"
            logger.info("  - SHP 坐标系: %s (单位: %s)", spatial_ref.GetName(), unit_name)
        else:
            raise RuntimeError("SHP 文件缺少 .prj 投影信息，无法确定坐标系统")
            
        shp_ds = None # 关闭

        # ③ 获取 JPG 像素尺寸
        jpg_ds = gdal.Open(out_jpg_path)
        if jpg_ds is None: raise RuntimeError("无法打开 JPG")
        width_jpg, height_jpg = jpg_ds.RasterXSize, jpg_ds.RasterYSize
        jpg_ds = None

        # ④ 计算像元分辨率 (根据 SHP 范围计算)
        min_x, max_x, min_y, max_y = extent
        if (max_x <= min_x or max_y <= min_y
                or width_jpg <= 0 or height_jpg <= 0):
            raise ValueError(f"SHP 坐标范围异常: {extent}")

        xres_jpg = (max_x - min_x) / width_jpg
        yres_jpg = (max_y - min_y) / height_jpg
        logger.info(
            "  - 按 SHP 外接矩形计算 JPG 地理参考（假设 JPG 覆盖整个矩形且无旋转）"
        )

        # ⑤ 生成 .jgw (World File)：A=像元宽 D=旋转 B=旋转 E=-像元高 C/F=左上像元中心
        with open(jgw_path, "w") as f:
            f.write(f"{xres_jpg}\n")  # A: X 分辨率
            f.write("0.0\n")          # D: 旋转
            f.write("0.0\n")          # B: 旋转
            f.write(f"{-yres_jpg}\n") # E: Y 分辨率 (负值)
            f.write(f"{min_x + xres_jpg / 2.0:.17g}\n") # C: 左上像元中心 X
            f.write(f"{max_y - yres_jpg / 2.0:.17g}\n") # F: 左上像元中心 Y
        logger.info("  - 成功生成坐标: %s", os.path.basename(jgw_path))

        # ⑥ 生成 .prj (与 JPG 同名)，内容为 SHP 的 WKT 投影
        prj_path = jpg_prj_path
        with open(prj_path, "w", encoding="utf-8") as f:
            f.write(prj_wkt)
        logger.info("  - 成功生成投影: %s", os.path.basename(prj_path))
        logger.info("  - 地理参考写入完成，用时 %.1f s", time.time() - stage_started)

        # ⑦ 绘制 SHP 边界并输出带地理参考的 RGB GeoTIFF
        burn_shapefile_boundary(out_jpg_path, out_shp_path, extent, logger)
        create_georeferenced_geotiff(
            out_jpg_path, out_tif_path, extent, prj_wkt, logger
        )
        with open(os.path.join(scene_dir, BOUNDARY_RENDER_MARKER), "w",
                  encoding="utf-8") as marker:
            marker.write("SHP boundary burned into JPG\n")

        logger.info("场景处理完成: %s (耗时 %.1fs)", scene, time.time() - started)
        return {"scene": scene, "status": "成功", "elapsed": time.time() - started}

    except Exception as e:
        logger.exception("处理失败: %s", e)
        return {"scene": scene, "status": "失败", "elapsed": time.time() - started}

# =================================================================
# 批处理主流程
# =================================================================
def main(argv=None):
    """批处理入口：扫描 → （dry-run 仅打印计划）→ 逐景处理 → 汇总。"""
    args = parse_args(argv)
    out_root = os.path.abspath(args.out)
    os.makedirs(out_root, exist_ok=True)
    log_path = args.log or os.path.join(out_root, "step0_georef_jpg_shp.log")
    logger = setup_logging(log_path, args.verbose)

    logger.info("=" * 60)
    logger.info("ZY1E AHSI 步骤0：快视图 JPG 边界绘制与地理参考")
    logger.info("源目录  : %s", os.path.abspath(args.src))
    logger.info("输出目录: %s", out_root)
    logger.info("日志文件: %s", log_path)
    logger.info("=" * 60)

    patterns = tuple(args.pattern) if args.pattern else TAR_PATTERNS
    try:
        scenes = discover_scenes(os.path.abspath(args.src), patterns)
    except Exception as e:
        logger.error("扫描失败: %s", e)
        return 2

    if args.scene:
        import fnmatch
        scenes = [p for p in scenes if any(fnmatch.fnmatch(get_scene_name(p), pat) for pat in args.scene)]
    if args.limit > 0: scenes = scenes[:args.limit]
    if not scenes: 
        logger.error("未找到任何包含 JPG 和 SHP 的文件夹或压缩包")
        return 2

    if args.dry_run:
        for idx, p in enumerate(scenes, 1):
            logger.info("[%d/%d] %s -> %s", idx, len(scenes), get_scene_name(p), os.path.join(out_root, get_scene_name(p)))
        logger.info("dry-run 完成：共 %d 景，未处理文件", len(scenes))
        return 0

    results = []
    started = time.time()
    
    for idx, scene_path in enumerate(scenes, start=1):
        logger.info("="*50)
        logger.info("正在处理第 %d/%d 景数据", idx, len(scenes))
        results.append(process_scene(scene_path, out_root, args, logger))

    success = sum(1 for r in results if r["status"] == "成功")
    skipped = sum(1 for r in results if r["status"].startswith("跳过"))
    failed = sum(1 for r in results if r["status"] == "失败")
    
    logger.info("="*50)
    logger.info("批处理结束: 成功 %d, 跳过 %d, 失败 %d，总用时 %.1f s", success, skipped, failed, time.time() - started)
    for result in results:
        logger.info("  [%s] %s（%.1f s）", result["status"], result["scene"], result["elapsed"])
    if failed:
        logger.error("有 %d 景处理失败，详见日志: %s", failed, log_path)
        return 1
    logger.info("全部完成，成果位于: %s", out_root)
    return 0

if __name__ == "__main__":
    sys.exit(main())