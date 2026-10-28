from __future__ import annotations

import gzip
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmwrite


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data"
GEO = DATA / "GSE135893"
UCSC = DATA / "GSE135893_UCSC"
OUT = ROOT / "03_virtual_knockout" / "inputs"
OUT.mkdir(parents=True, exist_ok=True)

MATRIX = UCSC / "exprMatrix.tsv.gz"
META = GEO / "GSE135893_IPF_metadata.csv.gz"

COMPARTMENTS = {
    "IPF_Macrophages": {
        "celltypes": ["Macrophages"],
        "targets": ["S100A8", "S100A9"],
    },
    "IPF_Epithelial": {
        "celltypes": ["AT2", "Transitional AT2", "MUC5B+", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "KRT5-/KRT17+"],
        "targets": ["LCN2", "ACSL4", "PTGS2", "HIF1A"],
    },
    "IPF_Fibroblast": {
        "celltypes": ["Fibroblasts", "PLIN2+ Fibroblasts", "HAS1 High Fibroblasts", "Myofibroblasts"],
        "targets": ["LIFR", "ACSL4"],
    },
}


def read_genes(path: Path) -> list[str]:
    return [line.strip() for line in path.read_text(encoding="utf-8-sig").splitlines() if line.strip()]


def balanced_cells(meta: pd.DataFrame, celltypes: list[str], max_per_donor: int = 100) -> list[str]:
    frame = meta[(meta.Diagnosis == "IPF") & meta.celltype.isin(celltypes)].copy()
    selected: list[str] = []
    for sample, sf in frame.groupby("Sample_Name", sort=True):
        n = min(max_per_donor, len(sf))
        selected.extend(sf.sample(n=n, random_state=20261027)["Cell"].astype(str).tolist())
    return selected


def main() -> None:
    meta = pd.read_csv(META, compression="gzip").rename(columns=lambda x: "Cell" if x == "Unnamed: 0" else x)
    compartment_cells = {
        name: balanced_cells(meta, cfg["celltypes"])
        for name, cfg in COMPARTMENTS.items()
    }
    union_cells = sorted(set(cell for cells in compartment_cells.values() for cell in cells))

    with gzip.open(MATRIX, "rt", encoding="utf-8") as handle:
        header = handle.readline().rstrip("\n\r").split("\t")
    if len(header) != 114_397:
        raise ValueError(f"Unexpected matrix column count: {len(header)}")
    col_index = {cell: i for i, cell in enumerate(header)}
    missing = [cell for cell in union_cells if cell not in col_index]
    if missing:
        raise ValueError(f"{len(missing)} selected cells are absent from matrix header")

    use_positions = [0] + sorted(col_index[cell] for cell in union_cells)
    print(f"Reading {len(union_cells)} balanced IPF cells from {len(header)-1} total cells", flush=True)
    frame = pd.read_csv(
        MATRIX,
        sep="\t",
        compression="gzip",
        usecols=use_positions,
        index_col=0,
        dtype={cell: np.uint32 for cell in union_cells},
        engine="c",
    )
    frame.index = [str(x).split("|")[-1] for x in frame.index]
    if frame.index.duplicated().any():
        frame = frame.groupby(level=0, sort=False).sum()

    ferro = read_genes(DATA / "ferro_aging_genes.txt")
    ferro = ["PTRF" if gene == "CAVIN1" else gene for gene in ferro]
    extension = read_genes(DATA / "analysis_extension_genes.txt")
    always_keep = set(ferro + extension + [g for cfg in COMPARTMENTS.values() for g in cfg["targets"]])

    report = {}
    for name, cfg in COMPARTMENTS.items():
        cells = [cell for cell in compartment_cells[name] if cell in frame.columns]
        counts = frame[cells]
        detected = (counts > 0).sum(axis=1)
        min_cells = max(5, int(np.ceil(0.03 * len(cells))))
        eligible = detected >= min_cells
        eligible_counts = counts.loc[eligible]

        totals = np.maximum(eligible_counts.sum(axis=0).to_numpy(dtype=float), 1.0)
        lognorm = np.log1p(eligible_counts.to_numpy(dtype=float) * (10_000.0 / totals[None, :]))
        variances = pd.Series(lognorm.var(axis=1, ddof=1), index=eligible_counts.index)
        top = variances.sort_values(ascending=False).head(900).index.tolist()
        extras = [g for g in always_keep if g in counts.index and detected.loc[g] >= 3]
        genes = list(dict.fromkeys(top + sorted(extras)))
        if len(genes) > 1_100:
            protected = set(extras)
            genes = extras + [g for g in top if g not in protected][: 1_100 - len(extras)]

        out_dir = OUT / name
        out_dir.mkdir(parents=True, exist_ok=True)
        selected = counts.loc[genes, cells].astype(np.uint32)
        # R's Matrix::readMM does not accept SciPy's non-standard
        # "unsigned-integer" Matrix Market field. Counts are small enough for
        # signed 32-bit integers, so write the portable standard integer type.
        matrix_for_r = selected.to_numpy().astype(np.int32, copy=False)
        mmwrite(str(out_dir / "counts.mtx"), sparse.coo_matrix(matrix_for_r))
        (out_dir / "genes.tsv").write_text("\n".join(selected.index.astype(str)) + "\n", encoding="utf-8")
        (out_dir / "cells.tsv").write_text("\n".join(selected.columns.astype(str)) + "\n", encoding="utf-8")
        cell_meta = meta.set_index("Cell").loc[cells].reset_index()
        cell_meta.to_csv(out_dir / "cell_metadata.csv", index=False, encoding="utf-8-sig")

        target_info = {}
        for target in cfg["targets"]:
            target_info[target] = {
                "present": target in selected.index,
                "detected_cells": int((selected.loc[target] > 0).sum()) if target in selected.index else 0,
                "detection_pct": float((selected.loc[target] > 0).mean() * 100) if target in selected.index else 0.0,
            }
        report[name] = {
            "celltypes": cfg["celltypes"],
            "n_cells": len(cells),
            "n_donors": int(cell_meta.Sample_Name.nunique()),
            "n_genes": len(genes),
            "minimum_detection_cells": min_cells,
            "targets": target_info,
        }
        print(json.dumps({name: report[name]}, ensure_ascii=False), flush=True)

    (OUT / "input_preparation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
