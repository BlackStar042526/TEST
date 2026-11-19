# ZY1E AHSI 高光谱蚀变矿物提取流水线

针对资源一号 02E（ZY1-02E / ZY1E）AHSI 高光谱数据的一套批处理脚本：从 L1A 压缩包出发，
完成辐射定标、坏波段剔除、正射、大气校正、研究区裁剪，以及 Al-OH 蚀变矿物提取与
MNF/PPI/MTMF 填图。

## 环境

- Python 3.9（conda 环境）
- GDAL（`osgeo`）、numpy、scipy、matplotlib、Pillow
- Py6S + 6S 可执行程序（仅步骤 2 需要）

建议用 conda 安装 GDAL，而非 pip：

```bash
conda create -n zy1e python=3.9
conda activate zy1e
conda install -c conda-forge gdal numpy scipy matplotlib pillow
pip install Py6S
```

## 目录结构

```
testcode/
├── 00_extract_georeferenced_jpg.py   快视图 JPG + SHP 边界 + 地理参考
├── 01_calibration_badbands.py         辐射定标 + 坏波段剔除 + RPC 正射
├── 02_atmospheric_correction.py       6S 大气校正（辐亮度 → 地表反射率）
├── 03_clip_by_area.py                 按研究区矢量裁剪（不重投影、不重采样）
├── 04_aloh_mineral_extraction.py      Al-OH 蚀变提取（SAM / SFF / MF）
├── 05_mnf_ppi_mtmf.py                 MNF + PPI + MTMF 填图
├── archive/                           历史版本（仅供追溯，不参与主链）
├── requirements.txt
└── README.md
```

## 处理流程

| 步骤 | 脚本 | 主要输入 | 主要输出 |
|---|---|---|---|
| 00 | `00_extract_georeferenced_jpg.py` | L1A `.tar` 快视图、研究区 `.shp` | `<景名>/` 带地理参考的 JPG |
| 01 | `01_calibration_badbands.py` | L1A `.tar` | `*_cal_badband_ortho.tif`（142 波段辐亮度，EPSG:32650 / 30 m）、`wavelengths_clean.npy`、`georef.json` |
| 02 | `02_atmospheric_correction.py` | 步骤 1 输出 + MODIS AOD + DEM | `*_cal_badband_ortho_atmcorr.tif`（地表反射率，与步骤 1 同网格） |
| 03 | `03_clip_by_area.py` | 步骤 2 输出 + 研究区 `.shp` | `<区名>/<景名>_cal_badband_ortho_atmcorr_clip.tif` |
| 04 | `04_aloh_mineral_extraction.py` | 步骤 3 区域镶嵌 + USGS splib07a 光谱库 | Al-OH 分类、置信度、吸收深度/位置/宽度、SAM/SFF/MF 成果 |
| 05 | `05_mnf_ppi_mtmf.py` | 步骤 3 区域镶嵌 + 步骤 4 光谱库 | MNF 分量、PPI 端元、MTMF 丰度/不可行性/分类、与 SAM 对比 |

运行顺序固定为 00 → 01 → 02 → 03，再执行 04 / 05。每步都支持 `--dry-run` 和
`--in` / `--out` 覆盖默认目录。

```bash
python 01_calibration_badbands.py --dry-run
python 01_calibration_badbands.py
python 02_atmospheric_correction.py
python 03_clip_by_area.py --area JLHK
python 04_aloh_mineral_extraction.py
python 05_mnf_ppi_mtmf.py
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

- 主链脚本统一为 `步骤号_功能.py`，不再在文件名中带版本号。
- 版本信息写在脚本头部 docstring 中，例如 `（v1.0.3）`。
- 历史版本放入 `archive/`，命名为 `<原功能名>_v<版本>.py`。

## 已知限制

- **步骤 3 / 4 / 5 的衔接**：当前主链的 `03_clip_by_area.py` 只做逐景裁剪，输出
  `<区名>/<景名>_..._clip.tif`；而 04 / 05 读取的是 `<区名>/mosaic/mosaic_reflectance.tif`
  （区域镶嵌）。若走"不拼接"的裁剪路线，需要调整 04 / 05 的输入读取逻辑。
  区域镶嵌实现见 `archive/03_orthorectify_clip_mosaic_v1.0.py`。
- 步骤 2 的 NASA LAADS 令牌请通过环境变量 `LAADS_TOKEN` 提供，不要写入源码。
