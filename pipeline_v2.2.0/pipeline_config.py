# -*- coding: utf-8 -*-
"""
pipeline_config.py —— 三步流水线统一配置（v2.2.0 单一真源）

把散落在 01/02/03 三个脚本里的默认路径、运行环境、输出网格常量收口到这里。
三个脚本 `import pipeline_config` 后取用；命令行参数仍可覆盖。

注意：
- 默认输出到 *_v21 新目录，与历史多代成果（142/158/166 波段混存）隔离；
- 观测几何目前只存在于 L1A 压缩包的 XML，因此 TAR_DIR 必须指向 ZY1E_AHSI；
- DEM 已覆盖 N27~N29 / E116~E118，还缺 N27E119、N28E119、N28E120、N29E119
  四块（补齐后 6S 海拔才不回落默认值）。
"""
import os

# ------------------------------------------------------------------
# 根目录与数据目录（按本机实际布局）
# ------------------------------------------------------------------
ROOT = r"E:\ZhouYang\TEST"

SRC_DIR = os.path.join(ROOT, "ZY1E_AHSI")            # L1A 压缩包目录（Step1 输入 / Step2 读观测几何）
DEM_DIR = os.path.join(ROOT, "DEM")                  # ASTGTM DEM（含 .img.zip）
AREA_ROOT = os.path.join(ROOT, "Hyperarea2")         # 研究区矢量（每子目录一个 .shp）
SPLIB_DIR = os.path.join(ROOT, "ASCIIdata_splib07a") # 下游 SAM/MTMF 光谱库（本三步暂不直接使用）
MOD04_DIR = os.path.join(ROOT, "MOD04_L2")           # MODIS AOD（仅留档，本版不启用）

# ------------------------------------------------------------------
# 输出根目录（与旧目录隔离）
# ------------------------------------------------------------------
OUT_ROOT_STEP1 = os.path.join(ROOT, "ZY1E_output_step1_v21")
OUT_ROOT_STEP2 = os.path.join(ROOT, "ZY1E_output_step2_v21")
OUT_ROOT_STEP3 = os.path.join(ROOT, "ZY1E_output_step3_v21")

# ------------------------------------------------------------------
# 运行环境（hsi 环境）
# ------------------------------------------------------------------
PYTHON = r"E:\ProgramData\anaconda3\envs\hsi\python.exe"
SIXS_EXE = r"E:\ProgramData\anaconda3\envs\hsi\Library\bin\sixs.exe"
# 6S 的 gfortran 运行库目录：仅在目录存在时注入 PATH；不存在则跳过（sixs.exe 若已
# 能自解 DLL 则无需注入）。留空表示完全跳过注入。
SIXS_DLL_DIR = r"E:\ProgramData\anaconda3\envs\hsi\Library\mingw-w64\bin"

# ------------------------------------------------------------------
# 输出网格与重采样
# ------------------------------------------------------------------
TARGET_EPSG = 32650        # WGS84 / UTM 50N（Step1 由足迹推导，一致）
TARGET_RES = 30.0          # 输出像元大小（米）
RESAMPLING = "bilinear"    # 反射率重采样（掩膜用 nearest）

# ------------------------------------------------------------------
# 波段参数
# ------------------------------------------------------------------
N_VNIR = 76
N_SWIR = 90
