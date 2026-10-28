# Ferro-aging-informed transcriptomic analysis of idiopathic pulmonary fibrosis

This repository contains the analysis code and processed result tables supporting the manuscript **Ferro-aging-informed multi-cohort transcriptomic analysis of idiopathic pulmonary fibrosis identifies a severity-linked epithelial remodeling state and a macrophage iron-lipid handling program**.

## Study overview

The analysis integrates bulk lung transcriptomes, single-cell RNA sequencing, virtual perturbation, independent single-cell validation, and spatial transcriptomics. The main result is a reproducible LCN2/SCGB-high and ACSL4/alveolar-low epithelial remodeling state associated with histological severity. The S100A8-associated macrophage lysosomal, lipid-processing, and ferritin program is treated as a secondary cell-intrinsic hypothesis.

## Public datasets

- GSE53845: bulk lung discovery cohort
- GSE32537: bulk lung external replication cohort
- GSE135893: single-cell discovery and cell-compartment analysis
- GSE136831: independent single-cell validation
- E-MTAB-14121: spatial transcriptomic severity analysis

Raw public data are not redistributed in this repository. They can be obtained from GEO and ArrayExpress/BioStudies using the accession identifiers above.

## Repository structure

- `scripts/bulk_ml/`: LASSO, SVM-RFE, and fixed-signature bulk analyses.
- `scripts/cell_and_spatial/`: single-cell, virtual perturbation, validation, and spatial analysis scripts.
- `results/bulk_ml/`: machine-learning summaries and out-of-fold results.
- `results/single_cell/`: donor-level pseudobulk and compartment-specific statistics.
- `results/virtual_perturbation/`: virtual knockout and overexpression result tables.
- `results/epithelial_state/`: epithelial pseudotime and independent validation tables.
- `results/myeloid_modules/`: fixed-module discovery and validation summaries.
- `results/spatial/`: donor- and section-level spatial severity summaries.

Large raw matrices, intermediate binary files, RDS objects, H5AD objects, and cell- or spot-level working matrices are excluded. The included tables are processed outputs used to reproduce the reported statistical summaries.

## Software

The principal analyses used R 4.6.1 with glmnet 5.1, e1071 1.7-17, and scTenifoldKnk 2.0.0, and Python 3.13.7 with Scanpy 1.12.2, anndata 0.13.2, harmonypy 0.0.9, scikit-learn 1.9.0, statsmodels 0.14.6, NumPy 2.4.6, pandas 2.3.3, and SciPy 1.18.0. Individual scripts may require additional packages documented in their import statements or session-information files.

## Reproducibility notes

1. Download the relevant public dataset into the `00_data/` directory expected by the scripts, or adjust the input path for your environment.
2. Run scripts in numerical order where numbered dependencies apply.
3. Use donors, rather than individual cells or spatial spots, as biological replicates for inferential comparisons.
4. Randomized procedures use the seeds recorded in the scripts and result metadata.

The scripts retain the analytical thresholds used for the manuscript. Because public data providers may update download endpoints or file packaging, input paths may require adjustment after download.

## Authors

Xiaohan Liu and Dong Zhang, Shandong University of Traditional Chinese Medicine, Jinan, Shandong, China.

## Contact

Dong Zhang: zhangdong@sdutcm.edu.cn

