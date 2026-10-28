from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import requests


ROOT = Path(__file__).resolve().parents[2]
KO_ROOT = ROOT / "03_virtual_knockout" / "results"
OUT = KO_ROOT / "pathway_enrichment"
OUT.mkdir(parents=True, exist_ok=True)

sources = ["GO:BP", "GO:CC", "REAC", "KEGG"]
jobs = []
for path in KO_ROOT.glob("*__*_KO/significant_responses_FDR05.csv"):
    tab = pd.read_csv(path)
    clean = tab.loc[
        ~tab["gene"].astype(str).str.match(r"^(RPL|RPS|MT-|LINC|MIR|AC\d)", flags=re.I),
        "gene",
    ].astype(str).tolist()
    jobs.append((path.parent.name, clean))

summaries = []
for name, genes in jobs:
    if len(genes) < 3:
        pd.DataFrame(columns=["source", "native", "name", "p_value", "intersection_size", "intersection_genes"]).to_csv(
            OUT / f"{name}_gProfiler.csv", index=False
        )
        summaries.append({"analysis": name, "input_genes": len(genes), "significant_terms": 0})
        continue

    payload = {
        "organism": "hsapiens",
        "query": genes,
        "sources": sources,
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
    }
    response = requests.post("https://biit.cs.ut.ee/gprofiler/api/gost/profile/", json=payload, timeout=90)
    response.raise_for_status()
    raw = response.json().get("result", [])
    rows = []
    for item in raw:
        evidences = item.get("intersections", [])
        hits = [gene for gene, evidence in zip(genes, evidences) if evidence]
        rows.append(
            {
                "source": item["source"],
                "native": item["native"],
                "name": item["name"],
                "p_value": item["p_value"],
                "intersection_size": item["intersection_size"],
                "term_size": item["term_size"],
                "precision": item["precision"],
                "intersection_genes": ";".join(hits),
            }
        )
    result = pd.DataFrame(rows).sort_values("p_value") if rows else pd.DataFrame()
    result.to_csv(OUT / f"{name}_gProfiler.csv", index=False, encoding="utf-8-sig")
    summaries.append({"analysis": name, "input_genes": len(genes), "significant_terms": len(result)})

    if len(result):
        top = result.head(15).sort_values("p_value", ascending=False).copy()
        top["score"] = -top["p_value"].clip(lower=1e-300).map(__import__("math").log10)
        top["label"] = top["name"].str.slice(0, 65)
        colors = {"GO:BP": "#4C78A8", "GO:CC": "#F58518", "REAC": "#54A24B", "KEGG": "#E45756"}
        fig, ax = plt.subplots(figsize=(9, 6.5))
        ax.barh(top["label"], top["score"], color=[colors.get(x, "#999999") for x in top["source"]])
        ax.set_xlabel("-log10(g:Profiler FDR)")
        ax.set_ylabel("")
        ax.set_title(f"{name}: non-ribosomal knockout-response enrichment")
        ax.grid(axis="x", alpha=0.25)
        fig.tight_layout()
        fig.savefig(OUT / f"{name}_gProfiler.png", dpi=300)
        fig.savefig(OUT / f"{name}_gProfiler.pdf")
        plt.close(fig)

pd.DataFrame(summaries).to_csv(OUT / "enrichment_summary.csv", index=False)
(OUT / "run_info.json").write_text(
    json.dumps(
        {
            "service": "g:Profiler gost API",
            "organism": "hsapiens",
            "sources": sources,
            "correction": "g:Profiler FDR",
            "filter": "FDR<0.05 knockout responses excluding ribosomal, mitochondrial, lncRNA and placeholder symbols",
        },
        indent=2,
    ),
    encoding="utf-8",
)
print(pd.DataFrame(summaries).to_string(index=False))

