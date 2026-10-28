#!/usr/bin/env python3
"""Prepare compact cell-to-donor/category maps for full-matrix NicheNet aggregation."""

from pathlib import Path
import json

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data" / "GSE136831"
OUT = ROOT / "06_mechanism_axis" / "nichenet_GSE136831"
OUT.mkdir(parents=True, exist_ok=True)

meta = pd.read_csv(
    DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
    sep="\t", compression="gzip",
)
barcodes = pd.read_csv(
    DATA / "GSE136831_AllCells.cellBarcodes.txt.gz",
    header=None, compression="gzip", dtype=str,
)[0].to_numpy()
if not np.array_equal(barcodes, meta["CellBarcode_Identity"].astype(str).to_numpy()):
    raise ValueError("Barcode and metadata order do not match")

donors = sorted(meta.loc[meta["Disease_Identity"].eq("IPF"), "Subject_Identity"].unique())
donor_to_index = {donor: i for i, donor in enumerate(donors)}

donor_index = np.full(len(meta), -1, dtype=np.int16)
category_bits = np.zeros(len(meta), dtype=np.uint8)

for i, row in meta.iterrows():
    if row["Disease_Identity"] != "IPF":
        continue
    donor_index[i] = donor_to_index[row["Subject_Identity"]]
    if row["Manuscript_Identity"] == "ATII":
        category_bits[i] |= 1
    if row["Manuscript_Identity"] in {"Club", "Goblet"}:
        category_bits[i] |= 2
    if row["CellType_Category"] == "Myeloid":
        category_bits[i] |= 4

donor_index.tofile(OUT / "cell_donor_index.int16.bin")
category_bits.tofile(OUT / "cell_category_bits.uint8.bin")

category_names = {1: "ATII", 2: "Secretory_Club_Goblet", 4: "Myeloid"}
rows = []
for bit, category in category_names.items():
    for donor, donor_i in donor_to_index.items():
        idx = np.flatnonzero((donor_index == donor_i) & ((category_bits & bit) != 0))
        rows.append({
            "donor_index": donor_i,
            "Subject_Identity": donor,
            "category_index": {1: 0, 2: 1, 4: 2}[bit],
            "category": category,
            "n_cells": len(idx),
            "total_nUMI": int(meta.loc[idx, "nUMI"].sum()),
        })
summary = pd.DataFrame(rows)
summary.to_csv(OUT / "donor_category_summary.csv", index=False, encoding="utf-8-sig")

run_info = {
    "dataset": "GSE136831",
    "n_cells": len(meta),
    "n_IPF_donors": len(donors),
    "category_bits": {"ATII": 1, "Secretory_Club_Goblet": 2, "Myeloid": 4},
    "minimum_cells_for_paired_DE": 5,
    "donors": donors,
}
(OUT / "map_run_info.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")

paired = summary.pivot(index="Subject_Identity", columns="category", values="n_cells").fillna(0)
paired = paired[(paired["ATII"] >= 5) & (paired["Secretory_Club_Goblet"] >= 5)]
print(summary.groupby("category").agg(n_cells=("n_cells", "sum"), n_donors=("n_cells", lambda x: (x > 0).sum())).to_string())
print(f"Paired IPF donors with >=5 ATII and >=5 secretory cells: {len(paired)}")
print(paired[["ATII", "Secretory_Club_Goblet"]].to_string())
