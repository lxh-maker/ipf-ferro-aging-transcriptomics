#!/usr/bin/env python3
"""Stream selected genes from the 692M-entry GSE136831 Matrix Market file.

The full sparse matrix is too large to materialize safely. This reader inflates
fixed-size blocks, parses them with NumPy, and retains only target rows. Output is
a compact dense target-gene x cell matrix (~tens of MB).
"""

import gzip
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
IN_DIR = ROOT / "00_data" / "GSE136831"
OUT = IN_DIR / "GSE136831_targeted_expression.npz"
BLOCK_BYTES = 64 * 1024 * 1024

TARGETS = [
    "LCN2", "ACSL4", "SCGB1A1", "SCGB3A1", "SFTPC", "SFTPA1",
    "KRT8", "KRT17", "KRT19", "MUC5B", "EPCAM", "VIM",
    "S100A8", "S100A9", "TLR4", "SPP1", "CD44", "LIF", "LIFR",
    "OSM", "OSMR", "AREG", "EGFR", "HMGB1",
]

gene_table = pd.read_csv(IN_DIR / "GSE136831_AllCells.GeneIDs.txt.gz", sep="\t", compression="gzip")
gene_symbol_col = "HGNC_EnsemblAlt_GeneID"
symbols = gene_table[gene_symbol_col].astype(str).tolist()

# Matrix Market uses 1-based row indices. If a symbol were duplicated, retain all
# matching rows and sum them into the same output symbol.
row_to_target = {}
present = []
for target_pos, gene in enumerate(TARGETS):
    matched = [i + 1 for i, symbol in enumerate(symbols) if symbol == gene]
    if matched:
        present.append(gene)
        for row in matched:
            row_to_target[row] = target_pos

missing = sorted(set(TARGETS) - set(present))
if missing:
    print("Missing target symbols:", ", ".join(missing), flush=True)

barcodes = pd.read_csv(
    IN_DIR / "GSE136831_AllCells.cellBarcodes.txt.gz", header=None,
    compression="gzip", dtype=str,
)[0].to_numpy()

# The companion optimized C++ extractor writes the 24 x 312,928 uint32 block.
# If present, package it directly and avoid rescanning the 692M-entry matrix.
fast_binary = IN_DIR / "GSE136831_targeted_expression.uint32.bin"
expected_binary_bytes = len(TARGETS) * len(barcodes) * np.dtype(np.uint32).itemsize
if fast_binary.exists():
    if fast_binary.stat().st_size != expected_binary_bytes:
        raise ValueError(
            f"Fast binary has {fast_binary.stat().st_size} bytes; expected {expected_binary_bytes}"
        )
    selected = np.fromfile(fast_binary, dtype=np.uint32).reshape(len(TARGETS), len(barcodes))
    np.savez_compressed(
        OUT,
        genes=np.asarray(TARGETS, dtype=str),
        expression=selected,
        cells=barcodes.astype(str),
    )
    report = {
        "source": "GSE136831_RawCounts_Sparse.mtx.gz",
        "matrix_shape": [45947, len(barcodes)],
        "declared_nonzero_entries": 692789348,
        "retained_target_entries": int((selected > 0).sum()),
        "genes_requested": TARGETS,
        "genes_present": TARGETS,
        "genes_missing": [],
        "extractor": "extract_gse136831_fast.cpp",
        "output": str(OUT),
    }
    (IN_DIR / "GSE136831_targeted_extraction_report.json").write_text(
        json.dumps(report, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, indent=2), flush=True)
    sys.exit(0)

matrix_path = IN_DIR / "GSE136831_RawCounts_Sparse.mtx.gz"
with gzip.open(matrix_path, "rb") as fh:
    header = fh.readline().decode().strip()
    if not header.startswith("%%MatrixMarket matrix coordinate"):
        raise ValueError(f"Unexpected Matrix Market header: {header}")
    line = fh.readline()
    while line.startswith(b"%"):
        line = fh.readline()
    n_rows, n_cols, n_entries = map(int, line.split())
    if n_cols != len(barcodes):
        raise ValueError(f"Matrix columns {n_cols} != barcodes {len(barcodes)}")

    selected = np.zeros((len(TARGETS), n_cols), dtype=np.uint32)
    target_rows = np.asarray(sorted(row_to_target), dtype=np.int32)
    row_target_map = np.full(n_rows + 1, -1, dtype=np.int32)
    for row, target_pos in row_to_target.items():
        row_target_map[row] = target_pos
    total_parsed = 0
    retained = 0
    remainder = b""
    started = time.time()
    block_number = 0

    while True:
        chunk = fh.read(BLOCK_BYTES)
        if not chunk:
            block = remainder
            remainder = b""
        else:
            block = remainder + chunk
            last_newline = block.rfind(b"\n")
            if last_newline < 0:
                remainder = block
                continue
            remainder = block[last_newline + 1 :]
            block = block[: last_newline + 1]

        if block:
            values = np.fromstring(block.decode("ascii"), dtype=np.int64, sep=" ")
            if values.size % 3:
                raise ValueError(f"Parsed value count is not divisible by 3 in block {block_number}")
            triples = values.reshape(-1, 3)
            mask = np.isin(triples[:, 0], target_rows, assume_unique=False)
            hits = triples[mask]
            if len(hits):
                np.add.at(
                    selected,
                    (row_target_map[hits[:, 0]], hits[:, 1] - 1),
                    hits[:, 2].astype(np.uint32),
                )
            total_parsed += len(triples)
            retained += len(hits)
            block_number += 1
            if block_number % 20 == 0:
                elapsed = time.time() - started
                print(
                    f"blocks={block_number}, entries={total_parsed:,}/{n_entries:,} "
                    f"({total_parsed/n_entries:.1%}), retained={retained:,}, elapsed={elapsed:.1f}s",
                    flush=True,
                )
        if not chunk:
            break

if total_parsed != n_entries:
    raise ValueError(f"Parsed {total_parsed:,} entries; Matrix Market declares {n_entries:,}")

keep_positions = [TARGETS.index(g) for g in present]
selected = selected[keep_positions]
np.savez_compressed(
    OUT,
    genes=np.asarray(present, dtype=str),
    expression=selected,
    cells=barcodes.astype(str),
)

report = {
    "source": "GSE136831_RawCounts_Sparse.mtx.gz",
    "matrix_shape": [n_rows, n_cols],
    "declared_nonzero_entries": n_entries,
    "parsed_entries": total_parsed,
    "retained_target_entries": retained,
    "genes_requested": TARGETS,
    "genes_present": present,
    "genes_missing": missing,
    "output": str(OUT),
    "elapsed_seconds": time.time() - started,
}
(IN_DIR / "GSE136831_targeted_extraction_report.json").write_text(
    json.dumps(report, indent=2), encoding="utf-8"
)
print(json.dumps(report, indent=2), flush=True)
