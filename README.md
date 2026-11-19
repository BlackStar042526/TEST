# ZY1E AHSI 高光谱蚀变矿物提取流水线

针对资源一号 02E（ZY1-02E / ZY1E）AHSI 高光谱数据的一套批处理脚本：从 L1A 压缩包出发，
完成辐射定标、坏波段剔除、正射、大气校正、研究区裁剪，以及 Al-OH 蚀变矿物提取与
MNF/PPI/MTMF 填图。

## 环境

本项目统一使用 conda 环境 `hsi`：

| 项 | 值 |
|---|---|
| 解释器 | `E:\ProgramData\anaconda3\envs\hsi\python.exe` |
| Python | 3.9.25 |
| GDAL | 3.10.3 |
| numpy | 2.0.2 |
| scipy | 1.13.1 |
| matplotlib | 3.9.2 |
| Pillow | 11.3.0 |
| Py6S | 已安装（步骤 2 还需要 6S 可执行程序） |

运行前激活环境：

```powershell
conda activate hsi
```

或在脚本里直接用解释器全路径：

```powershell
& E:\ProgramData\anaconda3\envs\hsi\python.exe .\01_calibration_badbands\01_calibration_badbands_v1.0.3.py
```


## 目录结构

每个处理步骤一个文件夹，文件夹内并列保存该步骤的全部版本，文件名带版本号，
因此可以同屏对比、按版本回退。

```
testcode/
├── 00_preview_jpg/
│   └── 00_extract_georeferenced_jpg_v1.0.0.py                       快视图 JPG + SHP 边界 + 地理参考
├── 01_calibration_badbands/
│   ├── 01_calibration_badbands_v1.0.0.py
│   ├── 01_calibration_badbands_v1.0.1.py
│   ├── 01_calibration_badbands_v1.0.2.py
│   └── 01_calibration_badbands_v1.0.3.py                            ← 当前版本
├── 02_atmospheric_correction/
│   ├── 02_atmospheric_correction_v1.0.0.py
│   ├── 02_atmcorrection_6s_v1.1.0.py
│   ├── 02_atmospheric_correction_v1.1.1.py
│   └── 02_atmospheric_correction_v1.1.2.py                          ← 当前版本
├── 03_clip_by_area/
│   ├── 03_clip_by_area_v2.0.py                                      ← 当前版本（仅裁剪）
│   └── 03_orthorectify_clip_mosaic_v1.0.py                          区域镶嵌实现（未参与主链）
├── 04_aloh_mineral_extraction/
│   └── 04_aloh_mineral_extraction.py                                Al-OH 蚀变提取
├── 05_mnf_ppi_mtmf/
│   └── 05_mnf_ppi_mtmf.py                                           MNF + PPI + MTMF
├── requirements.txt
└── README.md
```

04 / 05 目前各只有一个版本，原文件名未带版本号，故保持原样。

## 版本一览

| 步骤 | 当前使用 | 同目录其他版本 | 说明 |
|---|---|---|---|
| 00 | `00_extract_georeferenced_jpg_v1.0.0.py` | — | 唯一版本 |
| 01 | `01_calibration_badbands_v1.0.3.py` | v1.0.0 / v1.0.1 / v1.0.2 | 含 RPC 偏差修正与质检 |
| 02 | `02_atmospheric_correction_v1.1.2.py` | v1.1.1、v1.0.0、`02_atmcorrection_6s_v1.1.0.py` | 6S 大气校正 |
| 03 | `03_clip_by_area_v2.0.py` | `03_orthorectify_clip_mosaic_v1.0.py` | 主链只裁剪；镶嵌版未参与主链 |
| 04 | `04_aloh_mineral_extraction.py` | — | 唯一版本 |
| 05 | `05_mnf_ppi_mtmf.py` | — | 唯一版本 |

## 处理流程

| 步骤 | 主要输入 | 主要输出 |
|---|---|---|
| 00 | L1A `.tar` 快视图、研究区 `.shp` | `<景名>/` 带地理参考的 JPG |
| 01 | L1A `.tar` | `*_cal_badband_ortho.tif`（142 波段辐亮度，EPSG:32650 / 30 m）、`wavelengths_clean.npy`、`georef.json` |
| 02 | 步骤 1 输出 + MODIS AOD + DEM | `*_cal_badband_ortho_atmcorr.tif`（地表反射率，与步骤 1 同网格） |
| 03 | 步骤 2 输出 + 研究区 `.shp` | `<区名>/<景名>_cal_badband_ortho_atmcorr_clip.tif` |
| 04 | 步骤 3 区域镶嵌 + USGS splib07a 光谱库 | Al-OH 分类、置信度、吸收深度/位置/宽度、SAM/SFF/MF 成果 |
| 05 | 步骤 3 区域镶嵌 + 步骤 4 光谱库 | MNF 分量、PPI 端元、MTMF 丰度/不可行性/分类、与 SAM 对比 |

运行顺序固定为 00 → 01 → 02 → 03，再执行 04 / 05。每步都支持 `--dry-run` 和
`--in` / `--out` 覆盖默认目录。

```bash
conda activate hsi
cd E:\ZhouYang\test\code\testcode

python 00_preview_jpg/00_extract_georeferenced_jpg_v1.0.0.py
python 01_calibration_badbands/01_calibration_badbands_v1.0.3.py --dry-run
python 01_calibration_badbands/01_calibration_badbands_v1.0.3.py
python 02_atmospheric_correction/02_atmospheric_correction_v1.1.2.py
python 03_clip_by_area/03_clip_by_area_v2.0.py --area JLHK
python 04_aloh_mineral_extraction/04_aloh_mineral_extraction.py
python 05_mnf_ppi_mtmf/05_mnf_ppi_mtmf.py
```

## 数据目录约定

脚本中的默认路径统一为 `E:\ZhouYang\test\` 下的数据目录，**数据本身不纳入仓库**：

| 目录 | 用途 |
|---|---|
| `ZY1E_AHSI/` | L1A 原始压缩包（`.tar` / `.tar.gz`） |
| `ZY1E_output_step1/` | 步骤 1 成果 |
| `ZY1E_output_step2/` | 步骤 2 成果 |
| `ZY1E_output_step3/` | 步骤 3 裁剪成果 |
| `ZY1E_output_step4/` | 步骤 4 成果（含 `_library/`） |
| `ZY1E_output_step4_mnf/` | 步骤 5 成果 |
| `ZY1E_Georef_JPG/` | 步骤 0 成果 |
| `Hyperarea2/` | 研究区矢量（按区子目录存放 `.shp`） |
| `MOD04_L2/` | MODIS AOD 产品 |
| `DEM/` | ASTGTM DEM |
| `ASCIIdata_splib07a/` | USGS splib07a ASCII 光谱库 |

## 命名与版本约定

- 一步一个文件夹，命名 `<步骤号>_<功能>/`。
- 同一文件夹内，所有历史版本并列保留，文件名统一为 `<步骤号>_<功能>_v<主>.<次>.py`。
- 当前使用版本在本 README 的「版本一览」中标注，不额外加 `current` 前缀。
- 新增版本时复制上一个版本改名即可，不覆盖旧文件。

## 已知限制

- **步骤 3 / 4 / 5 的衔接**：主链的 `03_clip_by_area` 只做逐景裁剪，输出
  `<区名>/<景名>_..._clip.tif`；而 04 / 05 读取的是 `<区名>/mosaic/mosaic_reflectance.tif`
  （区域镶嵌）。若坚持"不拼接"的裁剪路线，需要调整 04 / 05 的输入读取逻辑。
  区域镶嵌实现见 `03_clip_by_area/03_orthorectify_clip_mosaic_v1.0.py`。
- 步骤 2 的 NASA LAADS 令牌请通过环境变量 `LAADS_TOKEN` 提供，不要写入源码。
