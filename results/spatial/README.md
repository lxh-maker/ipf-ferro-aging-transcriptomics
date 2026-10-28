# E-MTAB-14121 IPF 空间转录组验证

## 数据与设计

- 全转录组 Visium：E-MTAB-14121，共 25 张切片、101,907 个组织内空间点。
- IPF：4 位患者，B1/B2/B3 代表轻、中、重纤维化区域；技术重复先合并，再以患者为生物学重复。
- 对照：4 位健康供体。
- 细胞组成：作者发布的 Habermann 2020 参考 cell2location 5% 分位丰度（S-BSST1410）。
- 空间邻域：固定 Visium 阵列坐标的 6、12、24、48、96 个最近空间点。

## 核心结果

1. 上皮分泌重塑、LCN2 和综合 LCN2-high/ACSL4-low 状态均在 4/4 位患者中随 B1→B3 上升；ACSL4 在 4/4 位患者中下降。由于仅 4 位患者，双侧精确符号翻转检验最小 P=0.125，方向一致性比 P 值更有解释价值。
2. 重度 IPF 相对健康对照的髓系加权 TGFB1 中位数为 0.515 对 0.220，Cliff's delta=1.0，精确 P=0.0286、多指标 FDR=0.114；但多尺度空间相关均未显示 TGFB1 邻域与异常上皮状态正相关。因此 TGFB1 可保留为 NicheNet 候选和严重度伴随信号，不能写成已获空间支持的局部上游驱动。
3. SPP1 髓系邻域与 cell2location 定义的异常上皮丰度随邻域扩大而增强。在 96 个最近空间点尺度，4/4 位重度 IPF 患者相关系数均为正，中位 rho=0.193；对照中位 rho=0.023，Cliff's delta=0.875，未校正精确 P=0.0571。该多尺度探索结果支持“SPP1-rich 髓系—异常上皮生态位”，但没有达到传统显著性阈值，也不是配体因果证明。
4. S100A8/铁衰老髓系模块在任何预设空间尺度均未与异常上皮状态呈稳定正相关。S100A8 最稳妥的定位仍是髓系细胞内部的溶酶体—脂质—铁稳态模块。

## 主要文件

- `Figure_1_spatial_severity_trajectories.*`：患者内 B1→B3 程序轨迹。
- `Figure_3_representative_B3_spatial_maps.*`：代表性重度 IPF 切片空间图。
- `Figure_4_spatial_scale_sensitivity.*`：TGFB1、S100A8、SPP1 多尺度敏感性。
- `Figure_5_SPP1_abnormal_epithelial_niche.*`：SPP1-rich 髓系—异常上皮生态位重点图。
- `severity_trend_exact_tests.csv`：患者内严重度趋势精确检验。
- `IPF_B3_vs_control_donor_tests.csv`：重度 IPF 与对照的供体层面比较。
- `spatial_scale_sensitivity_B3_vs_control.csv`：多尺度供体层面效应量与精确检验。
- `spot_scores_and_neighborhoods.csv.gz`：全部空间点模块、去卷积和邻域分数。

## 可重复性与限制

- 代码：`01_scripts/17_spatial_ETAB14121_analysis.py` 和 `01_scripts/18_spatial_scale_sensitivity.py`。
- 固定随机种子：20261027；每张切片核心邻域检验 200 次置换。
- cell2location 压缩包 MD5：`03a0baaf460ba7ff2c1940b42536463e`。
- EBI 中 `V10T03-280-D1_spatial_tissue_positions_list.csv` 错误指向 HDF5 内容。邻域计算使用 Visium 固定条形码—阵列坐标从有效切片恢复；该切片不用于组织图像叠加。处理方式记录在 `section_level_spatial_scores.csv` 的 `position_source` 字段。
- Visium 每个空间点混合多个细胞；去卷积和空间相关只能支持组织生态位关联，不能代替受体阻断、共培养或空间蛋白验证。

## 数据来源

- E-MTAB-14121：https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-14121
- S-BSST1410：https://www.ebi.ac.uk/biostudies/studies/S-BSST1410
- 论文：https://www.nature.com/articles/s41588-024-01819-2
