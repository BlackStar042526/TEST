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
& E:\ProgramData\anaconda3\envs\hsi\python.exe .\01_calibration_badbands\01_calibration_badbands_v1.0.5.py
```


## 主流水线（pipeline_v2.2.0，推荐）

`pipeline_v2.2.0/` 是当前主用的 3 步流水线，用集中配置文件 `pipeline_config.py`
统一管理路径、环境与常量，前一步产物直接传给下一步：

| 步骤 | 脚本 | 作用 | 输出 |
|---|---|---|---|
| 1 | `01_calibration_v2.2.0.py` | L1A 解包 → 辐射定标（`L=Gain·DN+Bias`）→ DN=0 掩膜 → **RPC 足迹角点反算网格修正**（修正 v2.0 仿射往返恒等造成的约 2.4% 网格偏差） | `ZY1E_output_step1_v21/<景>/radiance_native.tif`（166 波段，原生网格 + RPC） |
| 2 | `02_atmospheric_correction_v2.2.0.py` | **6S/Py6S 大气校正**：K-T 暗像元 AOT（Visual 兜底）→ **高程分块 6S**（按 DEM 逐像元海拔分 3 档分别算系数）→ 反射率反演 → 水汽带剔除 | `ZY1E_output_step2_v21/<景>/…_atmcorr.tif`（147 波段） |
| 3 | `03_ortho_clip_mosaic_v2.2.0.py` | RPC 正射 + 研究区裁剪 + 云感知镶嵌：**主景优先 + 逐波段兜底**、剔除 6 个无用水汽带（147→**141**）、残留 NoData **光谱补值**（缺口 ≤3 线性插值、>3 拷贝兜底）、**S-G 3 点平滑 + 诊断峰核验** | `ZY1E_output_step3_v21/<区>/mosaic/mosaic_reflectance.tif`（141 波段） |

v2.2.0 关键精度改进：

- **几何**：RPC 网格修正改为用足迹 4 角点反算 SAMP/LINE 的 scale+offset，消除 v2.0 仿射往返导致的约 2.4% 网格偏差（残差 ≈0.001 像元）。
- **大气**：高程分块 6S——山区像元用对应海拔的 6S 系数，而不是整景一个平均海拔。
- **光谱连贯**：镶嵌时每个像元优先使用同一"主景"，仅缺失波段才逐波段补，减少跨景拼接台阶。
- **光谱保真**：只剔除 6 个无信息水汽带（1.4253/1.4422/1.4590/1.4759/1.7954/1.8120 μm），2.0–2.5 μm 等有用波段保留并补值；S-G 平滑后 2200/2340/1100 nm 诊断吸收峰深度变化 0.0%。

运行（先在 `pipeline_config.py` 把 `ROOT` 改成本机数据目录）：

```powershell
& E:\ProgramData\anaconda3\envs\hsi\python.exe pipeline_v2.2.0\01_calibration_v2.2.0.py
& E:\ProgramData\anaconda3\envs\hsi\python.exe pipeline_v2.2.0\02_atmospheric_correction_v2.2.0.py
& E:\ProgramData\anaconda3\envs\hsi\python.exe pipeline_v2.2.0\03_ortho_clip_mosaic_v2.2.0.py
```

> `pipeline/`（v2.0.0）与 `pipeline_v2.1.0/` 为历史版本，保留用于对照。

---

## 早期版本（00~05，单步脚本）

### 目录结构

每个处理步骤一个文件夹，文件夹内并列保存该步骤的全部版本，文件名带版本号，
因此可以同屏对比、按版本回退。

```
testcode/
├── 00_preview_jpg/
│   ├── 00_extract_georeferenced_jpg_v1.0.0.py                       原始版
│   ├── 00_extract_georeferenced_jpg_v1.0.1.py                       日志补阶段耗时与汇总
│   └── 00_extract_georeferenced_jpg_v1.0.2.py                       ← 当前版本（逐行行尾备注）
├── 01_calibration_badbands/
│   ├── 01_calibration_badbands_v1.0.0.py
│   ├── 01_calibration_badbands_v1.0.1.py
│   ├── 01_calibration_badbands_v1.0.2.py
│   ├── 01_calibration_badbands_v1.0.3.py
│   ├── 01_calibration_badbands_v1.0.4.py                            日志/标注
│   └── 01_calibration_badbands_v1.0.5.py                            ← 当前版本（逐行行尾备注）
├── 02_atmospheric_correction/
│   ├── 02_atmospheric_correction_v1.0.0.py
│   ├── 02_atmcorrection_6s_v1.1.0.py
│   ├── 02_atmospheric_correction_v1.1.1.py
│   ├── 02_atmospheric_correction_v1.1.2.py
│   ├── 02_atmospheric_correction_v1.1.3.py                          日志/标注
│   └── 02_atmospheric_correction_v1.1.4.py                          ← 当前版本（逐行行尾备注）
├── 03_clip_by_area/
│   ├── 03_clip_by_area_v2.0.py
│   ├── 03_clip_by_area_v2.1.py                                      日志/标注
│   ├── 03_clip_by_area_v2.2.py                                      ← 当前版本（逐行行尾备注）
│   └── 03_orthorectify_clip_mosaic_v1.0.py                          区域镶嵌实现（未参与主链）
├── 04_aloh_mineral_extraction/
│   ├── 04_aloh_mineral_extraction_v1.0.0.py                         原始版（仅支持 mosaic 输入）
│   ├── 04_aloh_mineral_extraction_v1.1.0.py                         支持裁剪成果输入
│   ├── 04_aloh_mineral_extraction_v1.1.1.py                         增强运行日志
│   ├── 04_aloh_mineral_extraction_v1.1.2.py                         代码标注
│   ├── 04_aloh_mineral_extraction_v1.1.3.py                         精简开头注释 + 阶段耗时
│   ├── 04_aloh_mineral_extraction_v1.1.4.py                         逐行行尾备注
│   └── 04_aloh_mineral_extraction_v1.1.5.py                         ← 当前版本（精简矿物清单：云母族+绿泥石族+蒙脱石）
├── 05_mnf_ppi_mtmf/
│   ├── 05_mnf_ppi_mtmf_v1.0.0.py                                    原始版（仅支持 mosaic 输入）
│   ├── 05_mnf_ppi_mtmf_v1.1.0.py                                    支持裁剪成果输入
│   ├── 05_mnf_ppi_mtmf_v1.1.1.py                                    增强运行日志
│   ├── 05_mnf_ppi_mtmf_v1.1.2.py                                    代码标注
│   ├── 05_mnf_ppi_mtmf_v1.1.3.py                                    精简开头注释
│   └── 05_mnf_ppi_mtmf_v1.1.4.py                                    ← 当前版本（逐行行尾备注）
├── requirements.txt
└── README.md
```

04 / 05 的 v1.0.0 为接入裁剪成果输入前的原始实现，保留用于对照与回退；
v1.1.0 增加输入适配并统一代码风格；v1.1.1 增强运行日志（阶段耗时、进度、
错误堆栈、失败区域汇总）；v1.1.2 补充代码结构与算法注释（函数说明、步骤编号、
关键公式解释），处理算法未变。

00～03 的最新版本同样只做“整理”：补充函数中文说明与步骤编号（①…⑦）、
阶段耗时/进度日志、失败场景汇总与日志文件路径；算法、输出文件名和格式
与上一版一致，原版本全部保留。

04 / 05 的 v1.1.3 按同一风格把开头的大段说明压缩为简短 docstring（详细内容
保留在配置区分节注释与函数说明中），并在 04 的区域处理流程中补上各阶段耗时
日志；算法与输出未变。

当前 6 个脚本（00 v1.0.2、01 v1.0.5、02 v1.1.4、03 v2.2、04 v1.1.5、
05 v1.1.4）统一补充了**逐行行尾中文备注**：在赋值、判断、循环、GDAL/NumPy
调用、写盘、异常分支等实质代码行末标注用途，便于逐行阅读；原有的分节注释、
①…⑦ 步骤标记与函数说明全部保留。其中 **04 v1.1.5 进一步将目标矿物由 11 种
Al-OH 精简为 6 种（云母族 + 绿泥石族 + 蒙脱石）**，剔除黄玉等非目标矿物，消除
黄玉大面积误匹配；其余算法与输出文件名不变。

## 版本一览

| 步骤 | 当前使用 | 同目录其他版本 | 说明 |
|---|---|---|---|
| 00 | `00_extract_georeferenced_jpg_v1.0.2.py` | v1.0.0 / v1.0.1 | v1.0.2 逐行行尾备注 |
| 01 | `01_calibration_badbands_v1.0.5.py` | v1.0.0 / v1.0.1 / v1.0.2 / v1.0.3 / v1.0.4 | v1.0.5 逐行行尾备注 |
| 02 | `02_atmospheric_correction_v1.1.4.py` | v1.1.3、v1.1.2、v1.1.1、v1.0.0、`02_atmcorrection_6s_v1.1.0.py` | v1.1.4 逐行行尾备注 |
| 03 | `03_clip_by_area_v2.2.py` | v2.1、v2.0、`03_orthorectify_clip_mosaic_v1.0.py` | v2.2 逐行行尾备注；镶嵌版未参与主链 |
| 04 | `04_aloh_mineral_extraction_v1.1.5.py` | v1.0.0 / v1.1.0 / v1.1.1 / v1.1.2 / v1.1.3 / v1.1.4 | v1.1.5 精简矿物清单（云母族+绿泥石族+蒙脱石） |
| 05 | `05_mnf_ppi_mtmf_v1.1.4.py` | v1.0.0 / v1.1.0 / v1.1.1 / v1.1.2 / v1.1.3 | v1.1.4 逐行行尾备注 |

## 处理流程

| 步骤 | 主要输入 | 主要输出 |
|---|---|---|
| 00 | L1A `.tar` 快视图、研究区 `.shp` | `<景名>/` 带地理参考的 JPG |
| 01 | L1A `.tar` | `*_cal_badband_ortho.tif`（142 波段辐亮度，EPSG:32650 / 30 m）、`wavelengths_clean.npy`、`georef.json` |
| 02 | 步骤 1 输出 + MODIS AOD + DEM | `*_cal_badband_ortho_atmcorr.tif`（地表反射率，与步骤 1 同网格） |
| 03 | 步骤 2 输出 + 研究区 `.shp` | `<区名>/<景名>_cal_badband_ortho_atmcorr_clip.tif` |
| 04 | 步骤 3 裁剪成果（或区域镶嵌）+ USGS splib07a 光谱库 | Al-OH 分类、置信度、吸收深度/位置/宽度、SAM/SFF/MF 成果 |
| 05 | 步骤 3 裁剪成果（或区域镶嵌）+ 步骤 4 光谱库 | MNF 分量、PPI 端元、MTMF 丰度/不可行性/分类、与 SAM 对比 |

运行顺序固定为 00 → 01 → 02 → 03，再执行 04 / 05。每步都支持 `--dry-run` 和
`--in` / `--out` 覆盖默认目录。

04 / 05 的输入适配：若 `<区名>/mosaic/` 三件套存在则直接使用；否则自动扫描
`<区名>/*_clip.tif`（含旧后缀 `_ortho_atmcorr_clip.tif`），用 `gdal.BuildVRT`
合并成虚拟镶嵌并按有效值生成掩膜，缓存到 `<步骤4输出>/_cache/<区名>/`，
波长取自步骤 1/2 同景目录的 `wavelengths_clean.npy`。加 `--overwrite` 可重建缓存。

```bash
conda activate hsi
cd E:\ZhouYang\test\code\testcode

python 00_preview_jpg/00_extract_georeferenced_jpg_v1.0.2.py
python 01_calibration_badbands/01_calibration_badbands_v1.0.5.py --dry-run
python 01_calibration_badbands/01_calibration_badbands_v1.0.5.py
python 02_atmospheric_correction/02_atmospheric_correction_v1.1.4.py
python 03_clip_by_area/03_clip_by_area_v2.2.py --area JLHK
python 04_aloh_mineral_extraction/04_aloh_mineral_extraction_v1.1.5.py
python 05_mnf_ppi_mtmf/05_mnf_ppi_mtmf_v1.1.4.py
```

## 数据目录约定

脚本中的默认路径统一为 `E:\ZhouYang\test\` 下的数据目录，**数据本身不纳入仓库**：

| 目录 | 用途 |
|---|---|
| `ZY1E_AHSI/` | L1A 原始压缩包（`.tar` / `.tar.gz`） |
| `ZY1E_output_step1/` | 步骤 1 成果 |
| `ZY1E_output_step2/` | 步骤 2 成果 |
| `ZY1E_output_step3/` | 步骤 3 裁剪成果 |
| `ZY1E_output_step4/` | 步骤 4 成果（含 `_library/`、`_cache/`） |
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

- 步骤 4 / 5 已内置输入适配，可直接消费步骤 3 的逐景裁剪成果（首次运行会建立
  `_cache/<区名>/` 下的 VRT 与掩膜，不落地大栅格），无需先做区域镶嵌。
- 步骤 2 的 NASA LAADS 令牌请通过环境变量 `LAADS_TOKEN` 提供，不要写入源码。
