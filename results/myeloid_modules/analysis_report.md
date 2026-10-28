# IPF 髓系状态连续体补充分析报告

## 结论先行

预设的四阶段线性模型：

`Resident/alveolar → S100A8 transitional → Lysosome/lipid/iron → SPP1 effector`

**未在独立队列中复现，因此不能作为论文的确定性主轨迹。** 数据更支持两个髓系分支：

1. `Resident/alveolar → lysosome–lipid–iron processing → SPP1-rich macrophage` 巨噬细胞重塑分支；
2. 以经典单核细胞为主的 `S100A8-high inflammatory` 分支。

两者可共同出现在 IPF 病理生态位中，但目前没有证据证明 S100A8-high 细胞必然沿单线转化为 SPP1-high 巨噬细胞。

## 方法

- 所有队列使用同一套事先固定的四个基因模块；独立验证阶段未重新筛基因或拟合权重。
- GSE135893：12 名 IPF 供体，每名 100 个巨噬细胞，共 1,200 个细胞；在 982 基因表达空间进行 Harmony 校正、扩散图和拟时序分析。
- GSE136831：94,679 个 IPF 髓系细胞、32 名患者；按患者和作者细胞类型平衡抽样 10,302 个细胞，进行固定模块状态投射。该投射用于描述状态邻近关系，不作为谱系追踪证据。
- E-MTAB-14121：25 张空间切片，4 名 IPF 患者和 4 名对照；技术重复先合并，再以供体为统计单位。分析模块严重度趋势、6/18/36 邻域的上皮—髓系相关及髓系富集区域内的模块共定位。

## 主要结果

### 1. GSE135893 发现队列仅部分支持原假说

| 相邻状态比较 | 拟时序位置差中位数 | 方向一致供体 | P 值 | FDR |
|---|---:|---:|---:|---:|
| Resident/alveolar → S100A8 transitional | +0.343 | 12/12 | 0.000488 | 0.00147 |
| S100A8 transitional → Lysosome/lipid/iron | -0.00349 | 5/12 | 0.970 | 0.970 |
| Lysosome/lipid/iron → SPP1 effector | +0.0101 | 10/12 | 0.0425 | 0.0637 |

- 四状态总体位置不同：Friedman P=9.11×10⁻⁵。
- 仅 4/12 名供体满足完整四阶段顺序。
- 因此可认为驻留型位于较早端，处理模块与 SPP1 端存在连续趋势；但 S100A8 与处理模块没有稳定的先后关系。

### 2. GSE136831 独立队列否定严格四阶段线性顺序

- 预设模块在 32 名患者的作者细胞类型中均呈显著差异，FDR 均小于 3×10⁻¹⁷。
- S100A8 模块在经典单核细胞中最高（患者伪总体中位数 1.29）。
- Resident/alveolar 模块在肺泡巨噬细胞中最高（0.67）。
- Lysosome/lipid/iron 与 SPP1 模块相对富集于巨噬细胞区室，而非经典单核细胞。
- 固定模块投射中：
  - S100A8 晚于 Resident/alveolar：32/32 名患者；
  - Lysosome/lipid/iron 晚于 S100A8：0/32 名患者，实际方向相反；
  - SPP1 晚于 Lysosome/lipid/iron：32/32 名患者。
- 完整四阶段顺序：0/32 名患者。

这说明 S100A8-high 状态主要代表炎症单核细胞分支，不能在本数据中被定义为处理型/SPP1 型巨噬细胞的统一前体。

### 3. 空间数据支持处理模块与 SPP1 模块共同形成髓系生态位，但不支持直接上皮驱动

- 四个模块的 B1→B3 供体内严重度趋势经精确符号翻转检验均不显著。
- 6、18、36 邻域中，各髓系模块与 LCN2-high/ACSL4-low 上皮转化分数的相关均不显著。
- 在髓系富集 spot 中，Lysosome/lipid/iron 与 SPP1 模块共定位：
  - 重度 IPF B3 中位偏相关 rho=0.309；
  - 对照中位 rho=0.147；
  - 4/4 名 B3 患者为正相关；
  - B3 对照精确 P=0.0286，四项共定位检验校正后 FDR=0.114。

该结果应表述为“探索性空间共定位证据”，不能表述为已证实的因果机制。

## 对论文主线的影响

不建议使用：

> S100A8-high transitional macrophages progress through a lysosomal–lipid–iron program into SPP1-high macrophages.

建议改为：

> IPF myeloid remodeling comprised a recruited S100A8-high inflammatory branch and a macrophage-centered lysosomal–lipid–iron/SPP1 branch. The latter showed enhanced spatial colocalization in severe fibrosis, whereas a uniform S100A8-to-SPP1 lineage and direct coupling to the complete LCN2-high/ACSL4-low epithelial program were not supported.

中文：

> IPF 髓系重塑包含 S100A8 高表达炎症单核细胞分支，以及以巨噬细胞为中心的溶酶体—脂质—铁处理/SPP1 分支。后者在重度纤维化中呈增强的空间共定位；但现有数据不支持统一的 S100A8 到 SPP1 线性谱系，也不支持该髓系模块直接驱动完整的 LCN2-high/ACSL4-low 上皮程序。

因此，S100A8 适合保留为炎症分支和虚拟扰动结果，放在次要结果或补充材料；正文的髓系核心更适合集中在 `GPNMB/LIPA/FTL/CTSD/ASAH1/HEXB/LAMP1/CYP27A1 → SPP1-rich macrophage niche`。

## 建议标题

**Ferro-aging-informed epithelial and macrophage remodeling reveals a lysosomal–lipid–iron/SPP1 niche in idiopathic pulmonary fibrosis**

中文：

**铁衰老框架下的上皮与巨噬细胞重塑揭示特发性肺纤维化中的溶酶体—脂质—铁处理/SPP1生态位**

## 解释边界

- 拟时序为计算状态排序，不是谱系追踪。
- GSE136831 的投射基于固定的 33 基因模块，适合做外部复现与反证，不应被描述为全转录组轨迹。
- 空间转录组 spot 为多细胞混合；cell2location 加权和近邻相关不能证明配体—受体因果关系。
- 空间供体数仅为 4+4，FDR 后证据未达到显著，需用多重免疫荧光、原位杂交或独立空间队列加强。

