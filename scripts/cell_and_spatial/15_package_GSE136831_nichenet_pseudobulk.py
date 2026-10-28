#!/usr/bin/env python3
"""Package C++ full-transcriptome aggregates as gene-symbol pseudobulk tables."""

from pathlib import Path
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data" / "GSE136831"
OUT = ROOT / "06_mechanism_axis" / "nichenet_GSE136831"
N_GENES = 45947
N_DONORS = 32
CATEGORIES = ["ATII", "Secretory_Club_Goblet", "Myeloid"]

summary = pd.read_csv(OUT / "donor_category_summary.csv")
donor_order = (
    summary[["donor_index", "Subject_Identity"]]
    .drop_duplicates().sort_values("donor_index")["Subject_Identity"].tolist()
)
if len(donor_order) != N_DONORS:
    raise ValueError(f"Expected {N_DONORS} donors, found {len(donor_order)}")

counts = np.fromfile(OUT / "pseudobulk_counts.uint64.bin", dtype=np.uint64)
detected = np.fromfile(OUT / "pseudobulk_detected.uint32.bin", dtype=np.uint32)
expected = len(CATEGORIES) * N_GENES * N_DONORS
if len(counts) != expected or len(detected) != expected:
    raise ValueError(f"Binary size mismatch: counts={len(counts)}, detected={len(detected)}, expected={expected}")
counts = counts.reshape(len(CATEGORIES), N_GENES, N_DONORS)
detected = detected.reshape(len(CATEGORIES), N_GENES, N_DONORS)

gene_table = pd.read_csv(DATA / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t", compression="gzip")
gene_symbols = gene_table["HGNC_EnsemblAlt_GeneID"].fillna("").astype(str)
valid = gene_symbols.str.len().gt(0) & ~gene_symbols.eq("nan")

qc_rows = []
for category_i, category in enumerate(CATEGORIES):
    count_frame = pd.DataFrame(counts[category_i, valid.to_numpy(), :], columns=donor_order)
    count_frame.insert(0, "gene", gene_symbols.loc[valid].to_numpy())
    count_frame = count_frame.groupby("gene", sort=True, as_index=False).sum(numeric_only=True)

    detect_frame = pd.DataFrame(detected[category_i, valid.to_numpy(), :], columns=donor_order)
    detect_frame.insert(0, "gene", gene_symbols.loc[valid].to_numpy())
    detect_frame = detect_frame.groupby("gene", sort=True, as_index=False).sum(numeric_only=True)

    count_frame.to_csv(OUT / f"pseudobulk_counts_{category}.tsv.gz", sep="\t", index=False, compression="gzip")
    detect_frame.to_csv(OUT / f"pseudobulk_detected_{category}.tsv.gz", sep="\t", index=False, compression="gzip")

    expected_depth = (
        summary.loc[summary["category"].eq(category)]
        .set_index("Subject_Identity").reindex(donor_order)["total_nUMI"].to_numpy(dtype=np.uint64)
    )
    matrix_depth = count_frame[donor_order].sum(axis=0).to_numpy(dtype=np.uint64)
    n_cells = (
        summary.loc[summary["category"].eq(category)]
        .set_index("Subject_Identity").reindex(donor_order)["n_cells"].to_numpy(dtype=int)
    )
    for donor, cells, expected_umis, observed_umis in zip(donor_order, n_cells, expected_depth, matrix_depth):
        qc_rows.append({
            "category": category, "Subject_Identity": donor, "n_cells": cells,
            "metadata_total_nUMI": int(expected_umis), "matrix_total_counts": int(observed_umis),
            "matrix_minus_metadata": int(observed_umis) - int(expected_umis),
        })

qc = pd.DataFrame(qc_rows)
qc.to_csv(OUT / "pseudobulk_depth_QC.csv", index=False, encoding="utf-8-sig")
report = {
    "categories": CATEGORIES,
    "n_gene_rows_input": N_GENES,
    "n_unique_gene_symbols": int(count_frame.shape[0]),
    "n_IPF_donors": N_DONORS,
    "max_absolute_depth_difference": int(qc["matrix_minus_metadata"].abs().max()),
    "notes": "Duplicate HGNC symbols were collapsed by summing counts and detection events.",
}
(OUT / "pseudobulk_packaging_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
print(json.dumps(report, indent=2))
print(qc.groupby("category")["matrix_minus_metadata"].agg(["min", "max"]).to_string())
