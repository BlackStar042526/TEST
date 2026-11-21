# ZY1E AHSI 三步处理流水线 · v2.2.0

> 迭代自 `pipeline_v2.1.0`。本目录是**可运行的迭代副本**，改动只在本目录内，`v2.1.0` 保持只读。
> 运行环境：`E:\ProgramData\anaconda3\envs\hsi\python.exe`（Python 3.9.25，GDAL 3.10.3，Py6S，numpy，matplotlib；6S 位于 `E:\ProgramData\anaconda3\envs\hsi\Library\bin\sixs.exe`）。

## 三步脚本与运行顺序

| 顺序 | 脚本 | 作用 | 关键输出 |
| --- | --- | --- | --- |
| 1 | `01_calibration_v2.2.0.py` | 辐射定标 `L=Gain·DN+Bias`；全 166 波段保留；**四角足迹反算的 RPC 校正**；推导目标 CRS | `radiance_native.tif`（native 网格，无投影仅 RPC）+ `georef.json`（含 `rpc_corrected`、`deferred_bad_bands`、`target_crs_epsg`） |
| 2 | `02_atmospheric_correction_v2.2.0.py` | K-T 暗像元 AOT 反演 + 逐波段 6S 校正 + 反射率域波段检测后剔除 | 147 波段反射率 tif + `band_detection_report.csv` + `coef_repair_log.json` + `used_parameters.json` |
| 3 | `03_ortho_clip_mosaic_v2.2.0.py` | RPC 正射 + 研究区裁剪 + 云感知镶嵌 | 各区 `mosaic/mosaic_reflectance.tif`（EPSG:32650/30m）+ 质检图 + `_catalog` |

标准运行（先单景验收，再全量）：

```bat
set PY=E:\ProgramData\anaconda3\envs\hsi\python.exe
set SCENE=ZY1E_AHSI_E116.45_N28.17_20250121_028090_L1A0000901923

%PY% 01_calibration_v2.2.0.py --scene %SCENE% --overwrite
%PY% 02_atmospheric_correction_v2.2.0.py --in E:\ZhouYang\TEST\ZY1E_output_step1_v21 --out E:\ZhouYang\TEST\ZY1E_output_step2_v21 --scene %SCENE% --overwrite
%PY% 03_ortho_clip_mosaic_v2.2.0.py --in E:\ZhouYang\TEST\ZY1E_output_step2_v21 --out E:\ZhouYang\TEST\ZY1E_output_step3_v21 --scene %SCENE% --overwrite
```

## v2.2.0 相对 v2.1.0 的改动

- **RPC 校正修复（几何精度，最重要）**：v2.0.0 起 `compute_rpc_grid_correction` 改成了「RPC 正变换→再反变换」往返拟合，数学上恒等（sx≈1），导致原始 RPC 与产品足迹之间约 2.4% 的格网偏差未被修正。本版恢复四角足迹反算：`sx/sy/tx/ty` 由足迹角点经 RPC 反变换后与影像四角拟合得到，实测 `sx≈1.0235, sy≈1.0215, tx≈-27.5, ty≈-39.5, 残差 0.0009 px`。
- **路径单一真源**：新增 `pipeline_config.py`，三步脚本的源目录/输出目录/环境/网格常量统一从它读取；默认路径指向 `E:\ZhouYang\TEST\...`，输出隔离到 `ZY1E_output_step{1,2,3}_v21`。
- **Step2 旧成果复用校验**：不再“旧文件名→新文件名”改名跳过；仅当已有成果波段数与当前输入波长一致才跳过，否则强制重算（避免 142/158 波段旧成果冒充 147）。
- **Step3 新鲜度/波段校验**：瓦片波段数与场景波长不一致、或上游反射率比瓦片新时自动重建瓦片；镶嵌前校验所有瓦片波段数一致；`mosaic_is_current` 增加波段一致性检查。
- **Step2 元数据契约**：`georef.json` 波段字段重写为剔除后口径；`coef_repair_log.json` 用剔除前波长（修复波段索引错位）；AOT/波段号口径留痕；新增 `--tar-dir`/`--dem-dir`；`--visibility 0` 显式校验。
- **Step3 研究区多要素**：按目录分组、优先主 shp、union 全部要素，避免只取第一个要素漏算多面研究区；修复 mosaic 组装时的 GDAL 句柄泄漏。
- **质检图失败即停**：Step2/Step3 补齐返回值检查，生成失败即抛异常（对齐 Step1 与 README 约定）。

## 数据与精度注意

- **DEM 缺口**：5 个研究区经度范围 116.5~119.5°E。现有 `DEM` 已覆盖 E116~E118；仅 `5.SC`（118.75~119.5°E）落入 E119，需补 **`ASTGTM_N28E119.img.zip`**（放入 `E:\ZhouYang\TEST\DEM\`）。E120 在研究区外，不需要；`N27E119 / N29E119` 仅为个别景足迹南北边缘，可选。
- **AOT**：K-T 暗像元为主；反演失败按 FLAASH「Initial Visibility」回退（默认 40 km，`--visibility` 覆盖），来源记入 `used_parameters.json`。
- **波段口径**：Step1 输出 166；Step2 反射率域剔除 19 个水汽带 → 147 波段；Step3 镶嵌时再剔除 6 个"完全无用"水汽带尾段（1.4 μm: 1.4253/1.4422/1.4590/1.4759；1.8 μm: 1.7954/1.8120）→ **141 波段**。有用波段（2.0~2.5 μm）保留，其残留 NoData 用光谱补值补救（连续缺口 ≤3 线性插值、>3 拷贝兜底），补值像元数记入 `mosaic_info.json`。
- **目标网格**：EPSG:32650 / 30 m（Step1 由足迹推导写入 georef，Step2/3 读取）。
- **重采样**：反射率 bilinear、掩膜 nearest。

## 验收要点（单景 E116.45 已跑通，2026-11-20）

- Step1：`rpc.grid_correction.sx=1.0235, sy=1.0215, tx=-27.5, ty=-39.5, 残差 0.0009 px`（2.4% 格网偏差已修正）；166 波段、目标 CRS EPSG:32650。
- Step2：DEM 平均海拔 55.4 m（修复多 UTM 带读取后，100% 覆盖）；AOT=0.38；166→**147 波段**；`coef_repair_log.json` 波段号与波长自洽；`georef.json` 波段口径重写为 147。
- Step3：`mosaic/mosaic_info.json.band_count=147`、EPSG:32650/30 m、覆盖 95.6%、清晰 95.6%。

全量前，仍需补齐 1 块 DEM：`ASTGTM_N28E119.img.zip`（覆盖 5.SC 研究区东半，放入 `E:\ZhouYang\TEST\DEM\`）。
