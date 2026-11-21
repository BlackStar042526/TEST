# ZY1E AHSI 三步处理流水线 · v2.1.2 迭代副本

> 迭代自 `../pipeline/`（v2.0.0 原码，未改动）。本目录是**可运行的迭代副本**，所有改动只在本目录内。
> 重建日期：2026-11-20。测试解释器：`C:\ProgramData\anaconda3\python.exe`（GDAL 3.11.4 + Py6S + numpy 2.5.5）。

## 三步脚本与运行顺序

| 顺序 | 脚本 | 作用 | 关键输出 |
| --- | --- | --- | --- |
| 1 | `01_calibration_v2.1.0.py` | 辐射定标 L=Gain·DN+Bias；**全 166 波段保留**（坏波段剔除后移）；RPC 校正 | `*_native_radiance.tif`（native 网格，无投影仅 RPC）+ `georef.json`（含 `deferred_bad_bands`、`target_crs_epsg=32650`） |
| 2 | `02_atmospheric_correction_v2.1.0.py` | K-T 暗像元 AOT 反演 + 逐波段 6S 校正 + **反射率域波段检测后剔除** | 158 波段反射率 tif + `band_detection_report.csv` + `coef_repair_log.json` + `preview_quality.png` |
| 3 | `03_ortho_clip_mosaic_v2.1.0.py` | RPC 正射 + 研究区裁剪 + 云感知镶嵌 | 各区 `mosaic_reflectance.tif`（EPSG:32650/30m）+ 质检图 + `_catalog` 编目 |

标准运行（单景示例，路径已指向本机真实数据）：

```bat
set PY=C:\ProgramData\anaconda3\python.exe
set SCENE=ZY1E_AHSI_E116.45_N28.17_20250121_028090_L1A0000901923

%PY% 01_calibration_v2.1.0.py --scene %SCENE%
%PY% 02_atmospheric_correction_v2.1.0.py --scene %SCENE%
%PY% 03_ortho_clip_mosaic_v2.1.0.py --scene %SCENE%
```

先 `--dry-run` 演练、`--limit 1` 限量、`--verbose` 看细节；确认无误再去掉限制跑全量。

## v2.1.2 相对原码（v2.0.0）的改动要点

- **波段剔除后移**：Step 1 不再物理删波段，Step 2 大气校正后在反射率域检测再剔（先出报告、按数据判定）。
- **坐标单一真源**：Step 1 由足迹推导目标 CRS（本批 EPSG:32650）写入 `georef.json`；Step 2/3 一律读取，`--epsg`/`--res` 可覆盖；不再硬编码。
- **错误即停**：关键路径失败立即非零退出（默认），`--continue-on-error` 可显式关闭；质检图必需，失败即停。
- **K-T AOT 反演修正（v2.1.2）**：修复旧版"地表反射率=0 的 6S 正向 LUT"致 AOT 恒顶 0.6 的 bug（实测 rho=0→6S辐亮度=0）；改用 6S 系数反演 ρ_swir/ρ_red 并匹配 FLAASH 经验比值 0.45（比值法），暗像元加 red/NIR<0.6 判据排除阴影/水体。单景实测 AOT 从 0.60 → **0.38**（全幅），red/SWIR=0.436 命中期望 0.45，物理自洽。不做外部定量验证（用户决定），反演失败统一走初始能见度回退链路（默认40km，可 --visibility 覆盖，来源留痕可追溯）。
- **水汽吸收带剔除（v2.1.2）**：显式剔除 0.94/1.14/1.38/1.88μm 强吸收带（对齐 FLAASH bad bands），消除校正后反射率尖峰；`band_detection_report.csv` 分列 `water_vapor` / `dn_zero` / `low_snr` 来源。
- **FLAASH Initial Visibility 对齐**：K-T 反演失败时用初始能见度（默认 **40 km**）作为初始值参与校正，`--visibility` 可覆盖，来源记入 `used_parameters.json`（`aot_source=initial_visibility_fallback`），非静默。
- **6S 运行保障**：`set_sensor_satellite_level()`（ZY1E 778km 正确接口）；启动注入 `mingw-w64/bin` 解决 `libgfortran-3.dll` 缺失；NaN 坏点走"±1nm 微调 → 邻带一致性校验 → 插值 → raise"，全程记入 `coef_repair_log.json`。
- **native 网格识别**：Step 1 输出无投影网格（设计如此，正射在 Step 3），Step 2 正确识别、跳过 CRS/res 校验，不再对空投影调用 `ImportFromWkt` 崩溃。
- **路径改为本机真实数据**：`C:\Users\yangz\Desktop\TEST\ZY1E_AHSI` / `TEST\DEM` / `TEST\Hyperarea2`。

## 单景试跑状态（2026-11-20，E116.45/N28.17/2025-01-21，v2.1.2 定稿）

- Step1：166 波段全保留 + 目标 CRS=EPSG:32650。
- Step2：K-T 反演 AOT550=**0.38**（red/SWIR=0.436 命中 FLAASH 0.45）；显式剔除水汽带 19 个 → 输出 **147 波段**；质检图确认反射率光谱尖峰消失。
- Step3：mosaic EPSG:32650/30m、清晰覆盖 93.8%。
- 回退链路：`resolve_aot550` 失败注入测试通过（暗像元不足 / 6S 异常 / 用户通过 --visibility 指定能见度）——三种失败均正确回退到初始能见度（默认 40km），带 warning + source + reason + 参数文件留痕。
- 不做外部 AOT 定量验证（用户 21:22 决定），不实现方案 B（FLAASH-2007 迭代）。
- 详见交付文档：`C:\Users\yangz\Desktop\work\WorkBuddy\2026-11-20_1712_分析TRAE分享代码与网页流程\01_交付物\`（以 `..._v2.md` 为准）。

## 数据质量问题（均已修复，无需再决策）

1. ✅ AOT 顶格 0.6：根因=6S 正向 LUT 地表反射率=0（rho=0→辐亮度=0）；已改 FLAASH 系数比值法，反演物理自洽。
2. ✅ 水汽带尖峰：已显式剔除 0.94/1.14/1.38/1.88μm 四个强吸收带，报告分列 water_vapor 来源。

## 目录说明

- 三脚本即本目录全部内容；`__pycache__` 为编译缓存，可安全删除。
- 原码与备份在 `../pipeline/`（未改动）；其它 `../00_preview_jpg/`、`../01_calibration_badbands/` 等为早期流程目录，未触碰。
