#!/usr/bin/env python3
"""Donor-aware targeted myeloid-to-epithelial ligand-receptor analysis.

Interactions are obtained from the literature-supported OmniPath core network. Expression
support and IPF-control differences are evaluated at the biological-donor level.
This is deliberately described as targeted communication evidence rather than a
comprehensive CellChat/CellPhoneDB analysis because the source matrix contains a
preselected panel of genes.
"""

from io import BytesIO
from pathlib import Path
from datetime import datetime, timezone

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data"
OUT = ROOT / "06_mechanism_axis" / "myeloid_epithelial_communication"
OUT.mkdir(parents=True, exist_ok=True)

OMNIPATH_URL = "https://omnipathdb.org/interactions"
OMNIPATH_PARAMS = {
    "datasets": "omnipath",
    "format": "tsv",
    "genesymbols": "1",
    "fields": "sources,references",
}
MIN_CELLS = 10

MYELOID_TYPES = ["Macrophages", "Monocytes", "Proliferating Macrophages", "cDCs", "pDCs"]
EPITHELIAL_TYPES = [
    "AT1", "AT2", "Transitional AT2", "SCGB3A2+", "SCGB3A2+ SCGB1A1+",
    "MUC5B+", "MUC5AC+ High", "KRT5-/KRT17+", "Basal",
    "Proliferating Epithelial Cells",
]

# Mechanistically relevant interactions whose two genes are in the targeted panel.
FOCUS_PAIRS = [
    ("S100A8", "TLR4"),
    ("S100A9", "TLR4"),
    ("HMGB1", "TLR4"),
    ("SPP1", "CD44"),
    ("LIF", "LIFR"),
    ("OSM", "LIFR"),
    ("OSM", "OSMR"),
    ("AREG", "EGFR"),
]


def add_fdr(frame: pd.DataFrame, pcol: str = "p_value") -> pd.DataFrame:
    frame = frame.copy()
    frame["FDR_BH"] = np.nan
    keep = frame[pcol].notna()
    if keep.any():
        frame.loc[keep, "FDR_BH"] = multipletests(frame.loc[keep, pcol], method="fdr_bh")[1]
    return frame


def compare_groups(frame: pd.DataFrame, keys: list[str]) -> pd.DataFrame:
    rows = []
    for key, dat in frame.groupby(keys, observed=True, sort=False):
        key = key if isinstance(key, tuple) else (key,)
        ipf = dat.loc[dat["Diagnosis"].eq("IPF"), "interaction_score"].to_numpy()
        ctrl = dat.loc[dat["Diagnosis"].eq("Control"), "interaction_score"].to_numpy()
        if len(ipf) < 5 or len(ctrl) < 5:
            continue
        u, p = stats.mannwhitneyu(ipf, ctrl, alternative="two-sided")
        cliff = 2 * u / (len(ipf) * len(ctrl)) - 1
        rows.append({
            **dict(zip(keys, key)),
            "n_IPF_donors": len(ipf), "n_Control_donors": len(ctrl),
            "median_IPF": np.median(ipf), "median_Control": np.median(ctrl),
            "IPF_minus_Control": np.median(ipf) - np.median(ctrl),
            "cliffs_delta": cliff, "p_value": p,
        })
    result = pd.DataFrame(rows)
    return add_fdr(result) if not result.empty else result


# Download and preserve the exact interaction resource used.
response = requests.get(OMNIPATH_URL, params=OMNIPATH_PARAMS, timeout=120)
response.raise_for_status()
(OUT / "OmniPath_core_interactions_used.tsv").write_bytes(response.content)
omni = pd.read_csv(BytesIO(response.content), sep="\t")

dat = np.load(DATA / "GSE135893_UCSC" / "selected_expression.npz", allow_pickle=False)
genes = dat["genes"].astype(str)
expr = dat["expression"].astype(np.float64, copy=False)
cells = dat["cells"].astype(str)
gene_idx = {gene: i for i, gene in enumerate(genes)}

meta = pd.read_csv(DATA / "GSE135893" / "GSE135893_IPF_metadata.csv.gz", compression="gzip")
if "Cell" not in meta.columns:
    meta = meta.rename(columns={meta.columns[0]: "Cell"})
if not np.array_equal(cells, meta["Cell"].astype(str).to_numpy()):
    raise ValueError("Expression and metadata cell order do not match")
keep = meta["Diagnosis"].isin(["IPF", "Control"]).to_numpy()
meta = meta.loc[keep].reset_index(drop=True)
expr = expr[:, keep]
depth = meta["nCount_RNA"].to_numpy(dtype=float)

supported = omni[
    omni["source_genesymbol"].isin(genes) & omni["target_genesymbol"].isin(genes)
].copy()
supported.to_csv(OUT / "OmniPath_interactions_supported_by_targeted_panel.csv", index=False, encoding="utf-8-sig")

focus = pd.DataFrame(FOCUS_PAIRS, columns=["ligand", "receptor"])
annotations = (
    omni.rename(columns={"source_genesymbol": "ligand", "target_genesymbol": "receptor"})
    .merge(focus, on=["ligand", "receptor"], how="inner")
    .drop_duplicates(["ligand", "receptor"])
)
if len(annotations) != len(FOCUS_PAIRS):
    absent = set(FOCUS_PAIRS) - set(map(tuple, annotations[["ligand", "receptor"]].to_numpy()))
    raise RuntimeError(f"Focus interactions absent from the OmniPath core network: {sorted(absent)}")
annotations.to_csv(OUT / "focused_OmniPath_interactions.csv", index=False, encoding="utf-8-sig")


def mean_lognorm(gene: str, indices: np.ndarray) -> tuple[float, float]:
    values = expr[gene_idx[gene], indices]
    lognorm = np.log1p(values / depth[indices] * 1e4)
    return float(lognorm.mean()), float((values > 0).mean() * 100)


# Primary, robust analysis: pooled myeloid sender and epithelial receiver within
# each donor. It tests expression-supported communication, not physical binding.
broad_rows = []
for donor, donor_meta in meta.groupby("Sample_Name", sort=False):
    sender = donor_meta.index[donor_meta["celltype"].isin(MYELOID_TYPES)].to_numpy()
    receiver = donor_meta.index[donor_meta["celltype"].isin(EPITHELIAL_TYPES)].to_numpy()
    if len(sender) < MIN_CELLS or len(receiver) < MIN_CELLS:
        continue
    diagnosis = donor_meta["Diagnosis"].iloc[0]
    for ligand, receptor in FOCUS_PAIRS:
        ligand_mean, ligand_pct = mean_lognorm(ligand, sender)
        receptor_mean, receptor_pct = mean_lognorm(receptor, receiver)
        broad_rows.append({
            "Sample_Name": donor, "Diagnosis": diagnosis,
            "sender": "Myeloid", "receiver": "Epithelial",
            "ligand": ligand, "receptor": receptor,
            "n_sender_cells": len(sender), "n_receiver_cells": len(receiver),
            "ligand_mean_log1pCP10K": ligand_mean,
            "receptor_mean_log1pCP10K": receptor_mean,
            "ligand_pct_expressing": ligand_pct,
            "receptor_pct_expressing": receptor_pct,
            "interaction_score": np.sqrt(ligand_mean * receptor_mean),
        })
broad = pd.DataFrame(broad_rows)
broad["interaction"] = broad["ligand"] + "–" + broad["receptor"]
broad_stats = compare_groups(broad, ["ligand", "receptor"])
broad_stats["interaction"] = broad_stats["ligand"] + "–" + broad_stats["receptor"]
broad.to_csv(OUT / "donor_level_broad_interaction_scores.csv", index=False, encoding="utf-8-sig")
broad_stats.sort_values("p_value").to_csv(OUT / "broad_IPF_vs_Control_interaction_tests.csv", index=False, encoding="utf-8-sig")

# Secondary localization: specific sender and receiver cell types. Donors must
# contribute at least MIN_CELLS to both cell types; all tests share one BH family.
detailed_rows = []
for donor, donor_meta in meta.groupby("Sample_Name", sort=False):
    diagnosis = donor_meta["Diagnosis"].iloc[0]
    for sender_type in MYELOID_TYPES:
        sender = donor_meta.index[donor_meta["celltype"].eq(sender_type)].to_numpy()
        if len(sender) < MIN_CELLS:
            continue
        for receiver_type in EPITHELIAL_TYPES:
            receiver = donor_meta.index[donor_meta["celltype"].eq(receiver_type)].to_numpy()
            if len(receiver) < MIN_CELLS:
                continue
            for ligand, receptor in FOCUS_PAIRS:
                ligand_mean, ligand_pct = mean_lognorm(ligand, sender)
                receptor_mean, receptor_pct = mean_lognorm(receptor, receiver)
                detailed_rows.append({
                    "Sample_Name": donor, "Diagnosis": diagnosis,
                    "sender": sender_type, "receiver": receiver_type,
                    "ligand": ligand, "receptor": receptor,
                    "n_sender_cells": len(sender), "n_receiver_cells": len(receiver),
                    "ligand_mean_log1pCP10K": ligand_mean,
                    "receptor_mean_log1pCP10K": receptor_mean,
                    "ligand_pct_expressing": ligand_pct,
                    "receptor_pct_expressing": receptor_pct,
                    "interaction_score": np.sqrt(ligand_mean * receptor_mean),
                })
detailed = pd.DataFrame(detailed_rows)
detailed["interaction"] = detailed["ligand"] + "–" + detailed["receptor"]
detailed_stats = compare_groups(detailed, ["sender", "receiver", "ligand", "receptor"])
detailed_stats["interaction"] = detailed_stats["ligand"] + "–" + detailed_stats["receptor"]
detailed.to_csv(OUT / "donor_celltype_interaction_scores.csv.gz", index=False, compression="gzip")
detailed_stats.sort_values("p_value").to_csv(OUT / "celltype_resolved_IPF_vs_Control_interaction_tests.csv", index=False, encoding="utf-8-sig")

sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)

# Figure 1: donor-level broad-compartment comparison.
order = broad_stats.sort_values("cliffs_delta", ascending=False)["interaction"].tolist()
fig, ax = plt.subplots(figsize=(11, 5.5))
sns.boxplot(data=broad, x="interaction", y="interaction_score", hue="Diagnosis",
            order=order, hue_order=["Control", "IPF"], palette=["#4C78A8", "#E45756"],
            showfliers=False, width=.65, ax=ax)
sns.stripplot(data=broad, x="interaction", y="interaction_score", hue="Diagnosis",
              order=order, hue_order=["Control", "IPF"], dodge=True,
              palette=["#4C78A8", "#E45756"], alpha=.7, size=4, linewidth=.25, ax=ax)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="", frameon=False)
ax.set_xlabel(""); ax.set_ylabel("sqrt(mean ligand × mean receptor)")
ax.tick_params(axis="x", rotation=40)
ax.set_title("Donor-level myeloid → epithelial interaction support", fontweight="bold")
for xpos, interaction in enumerate(order):
    row = broad_stats.loc[broad_stats["interaction"].eq(interaction)].iloc[0]
    ymax = broad.loc[broad["interaction"].eq(interaction), "interaction_score"].max()
    label = f"FDR={row['FDR_BH']:.3g}"
    ax.text(xpos, ymax * 1.08 + .005, label, ha="center", va="bottom", fontsize=7, rotation=45)
ax.margins(y=.16)
fig.tight_layout()
fig.savefig(OUT / "Figure_1_broad_myeloid_epithelial_interactions.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_1_broad_myeloid_epithelial_interactions.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 2: most strongly altered cell-type-resolved combinations.
top = detailed_stats.sort_values("p_value").head(30).copy()
top["route"] = top["sender"] + " → " + top["receiver"]
top["label"] = top["interaction"] + " | " + top["route"]
top = top.sort_values("cliffs_delta")
fig, ax = plt.subplots(figsize=(10, max(7, .31 * len(top))))
sizes = np.clip(-np.log10(top["FDR_BH"].clip(lower=1e-12)) * 45, 20, 260)
points = ax.scatter(top["cliffs_delta"], top["label"], s=sizes,
                    c=top["IPF_minus_Control"], cmap="RdBu_r", linewidths=.4, edgecolors="black")
ax.axvline(0, color="grey", ls="--", lw=1)
ax.set_xlabel("Cliff's delta (IPF vs Control)"); ax.set_ylabel("")
ax.set_title("Cell-type localization of targeted interactions", fontweight="bold")
cbar = fig.colorbar(points, ax=ax, pad=.02)
cbar.set_label("Median score difference (IPF − Control)")
fig.tight_layout()
fig.savefig(OUT / "Figure_2_celltype_resolved_interactions.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_2_celltype_resolved_interactions.pdf", bbox_inches="tight")
plt.close(fig)

with open(OUT / "README.txt", "w", encoding="utf-8") as fh:
    fh.write(
        "Targeted myeloid-to-epithelial ligand-receptor analysis\n"
        f"Run UTC: {datetime.now(timezone.utc).isoformat()}\n"
        f"OmniPath request: {response.url}\n"
        f"Targeted panel genes: {len(genes)}; retained strict-cohort cells: {len(meta)}\n"
        f"Donors: {meta.loc[meta.Diagnosis.eq('IPF'), 'Sample_Name'].nunique()} IPF and "
        f"{meta.loc[meta.Diagnosis.eq('Control'), 'Sample_Name'].nunique()} controls.\n"
        "Expression: per-cell log1p(count/nCount_RNA*10,000); donor interaction score is the\n"
        "geometric mean of sender-ligand and receiver-receptor mean expression.\n"
        "Inference: two-sided Mann-Whitney tests across biological donors, with BH correction.\n"
        "The cell-type-resolved analysis requires >=10 cells of both types in a donor and is\n"
        "secondary localization. This targeted panel cannot establish global communication\n"
        "probabilities or causality.\n"
    )

print("Broad myeloid-to-epithelial tests:")
print(broad_stats.sort_values("p_value")[["interaction", "n_IPF_donors", "n_Control_donors",
      "median_IPF", "median_Control", "cliffs_delta", "p_value", "FDR_BH"]].to_string(index=False))
print("\nTop cell-type-resolved tests:")
print(detailed_stats.sort_values("p_value")[["sender", "receiver", "interaction",
      "n_IPF_donors", "n_Control_donors", "cliffs_delta", "p_value", "FDR_BH"]].head(20).to_string(index=False))
print(f"\nOutputs: {OUT}")
