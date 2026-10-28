from __future__ import annotations

import json
import struct
import time
import zlib
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data"
GEO = DATA / "GSE135893"
UCSC = DATA / "GSE135893_UCSC"
INDEX_FILE = UCSC / "exprMatrix.json"
BIN_URL = "https://cells.ucsc.edu/lung-pf-control/exprMatrix.bin"
N_CELLS = 114_396

GENE_FILES = [
    DATA / "ferro_aging_genes.txt",
    DATA / "candidate_genes_8.txt",
    DATA / "signature_genes_5.txt",
    DATA / "mechanism_extension_genes.txt",
    DATA / "analysis_extension_genes.txt",
]


def read_genes() -> list[str]:
    genes: list[str] = []
    for path in GENE_FILES:
        genes.extend(
            line.strip()
            for line in path.read_text(encoding="utf-8-sig").splitlines()
            if line.strip()
        )
    return sorted(set(genes))


def make_session() -> requests.Session:
    retry = Retry(
        total=8,
        connect=8,
        read=8,
        status=8,
        backoff_factor=1.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=("GET",),
    )
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retry, pool_connections=8, pool_maxsize=8))
    session.headers.update({"Accept-Encoding": "identity", "User-Agent": "IPF-ferro-aging-analysis/1.0"})
    return session


def fetch_gene(gene: str, offset: int, length: int) -> tuple[str, np.ndarray]:
    end = offset + length - 1
    last_error: Exception | None = None
    for attempt in range(5):
        try:
            session = make_session()
            response = session.get(
                BIN_URL,
                headers={"Range": f"bytes={offset}-{end}"},
                timeout=(30, 120),
            )
            response.raise_for_status()
            payload = response.content
            if len(payload) != length:
                raise ValueError(f"{gene}: expected {length} compressed bytes, got {len(payload)}")

            raw = zlib.decompress(payload)
            name_len = struct.unpack_from("<H", raw, 0)[0]
            stored_gene = raw[2 : 2 + name_len].decode("utf-8")
            if stored_gene != gene:
                raise ValueError(f"index mismatch: requested {gene}, decoded {stored_gene}")

            values = np.frombuffer(raw, dtype="<u4", offset=2 + name_len).copy()
            if values.size != N_CELLS:
                raise ValueError(f"{gene}: expected {N_CELLS} values, got {values.size}")
            return gene, values
        except Exception as exc:  # retry a short range safely
            last_error = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"Failed to fetch {gene}: {last_error}")


def main() -> None:
    index = json.loads(INDEX_FILE.read_text(encoding="utf-8"))
    requested = read_genes()
    present = [gene for gene in requested if gene in index]
    missing = [gene for gene in requested if gene not in index]

    results: dict[str, np.ndarray] = {}
    failures: dict[str, str] = {}
    with ThreadPoolExecutor(max_workers=6) as executor:
        futures = {
            executor.submit(fetch_gene, gene, int(index[gene][0]), int(index[gene][1])): gene
            for gene in present
        }
        for future in as_completed(futures):
            gene = futures[future]
            try:
                name, values = future.result()
                results[name] = values
                print(f"fetched {name}: detected in {(values > 0).sum()} cells", flush=True)
            except Exception as exc:
                failures[gene] = str(exc)
                print(f"FAILED {gene}: {exc}", flush=True)

    genes = sorted(results)
    matrix = np.vstack([results[gene] for gene in genes]).astype(np.uint32, copy=False)
    metadata = pd.read_csv(GEO / "GSE135893_IPF_metadata.csv.gz", compression="gzip")
    if metadata.shape[0] != N_CELLS:
        raise ValueError(f"metadata has {metadata.shape[0]} rows, expected {N_CELLS}")

    np.savez_compressed(
        UCSC / "selected_expression.npz",
        genes=np.asarray(genes, dtype=str),
        expression=matrix,
        cells=np.asarray(metadata.iloc[:, 0].astype(str).tolist(), dtype=str),
    )
    pd.DataFrame(
        {
            "gene": genes,
            "detected_cells": [(results[g] > 0).sum() for g in genes],
            "total_counts": [results[g].sum(dtype=np.uint64) for g in genes],
        }
    ).to_csv(UCSC / "selected_expression_qc.csv", index=False, encoding="utf-8-sig")

    report = {
        "requested_gene_count": len(requested),
        "downloaded_gene_count": len(genes),
        "missing_from_index": missing,
        "failures": failures,
        "cell_count": N_CELLS,
        "source": BIN_URL,
        "cell_order": "GSE135893_IPF_metadata.csv.gz first column",
    }
    (UCSC / "selected_expression_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
