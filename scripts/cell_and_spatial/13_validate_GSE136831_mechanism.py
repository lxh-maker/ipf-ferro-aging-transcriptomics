#!/usr/bin/env python3
"""Independent donor-aware validation of the IPF epithelial mechanism axis.

Dataset: GSE136831 (Adams et al.), 32 IPF and 28 control subjects. The primary
validation is a paired within-IPF comparison of secretory (Club/Goblet) versus
ATII states. A secondary analysis tests broad myeloid-to-epithelial interactions.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data" / "GSE136831"
OUT = ROOT / "06_mechanism_axis" / "external_validation_GSE136831"
OUT.mkdir(parents=True, exist_ok=True)

MIN_STATE_CELLS = 5
MIN_COMMUNICATION_CELLS = 20
STATE_GENES = ["LCN2", "ACSL4", "SCGB1A1", "SCGB3A1", "SFTPC", "SFTPA1", "KRT17", "KRT19", "MUC5B"]
EPITHELIAL_TYPES = ["ATII", "ATI", "Club", "Goblet", "Basal", "Aberrant_Basaloid", "Ciliated"]
FOCUS_PAIRS = [
    ("S100A8", "TLR4"), ("S100A9", "TLR4"), ("HMGB1", "TLR4"),
    ("SPP1", "CD44"), ("LIF", "LIFR"), ("OSM", "LIFR"),
    ("OSM", "OSMR"), ("AREG", "EGFR"),
]


def add_fdr(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["FDR_BH"] = multipletests(frame["p_value"], method="fdr_bh")[1]
    return frame


archive = np.load(DATA / "GSE136831_targeted_expression.npz", allow_pickle=False)
genes = archive["genes"].astype(str)
expr_all = archive["expression"].astype(np.float64, copy=False)
cells = archive["cells"].astype(str)
gene_idx = {gene: i for i, gene in enumerate(genes)}

meta_all = pd.read_csv(
    DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz",
    sep="\t", compression="gzip",
)
if not np.array_equal(cells, meta_all["CellBarcode_Identity"].astype(str).to_numpy()):
    raise ValueError("Expression and metadata cell order do not match")

keep = meta_all["Disease_Identity"].isin(["IPF", "Control"]).to_numpy()
meta = meta_all.loc[keep].reset_index(drop=True)
expr = expr_all[:, keep]
depth = meta["nUMI"].to_numpy(dtype=float)

# Primary external validation: within each IPF donor, aggregate mature ATII and
# secretory Club/Goblet cells separately and form a paired state contrast.
state_rows = []
ipf_meta = meta.loc[meta["Disease_Identity"].eq("IPF")]
for donor, donor_meta in ipf_meta.groupby("Subject_Identity", sort=False):
    for state, celltypes in {"ATII": ["ATII"], "Secretory": ["Club", "Goblet"]}.items():
        idx = donor_meta.index[donor_meta["Manuscript_Identity"].isin(celltypes)].to_numpy()
        if len(idx) < MIN_STATE_CELLS:
            continue
        total_umis = depth[idx].sum()
        sums = expr[:, idx].sum(axis=1)
        for gene in STATE_GENES:
            state_rows.append({
                "Subject_Identity": donor, "state": state, "gene": gene,
                "n_cells": len(idx),
                "log2CPM": np.log2(sums[gene_idx[gene]] / total_umis * 1e6 + 1),
            })
state_pb = pd.DataFrame(state_rows)
state_wide = (
    state_pb.pivot_table(index=["Subject_Identity", "gene"], columns="state", values="log2CPM")
    .dropna().reset_index()
)
state_wide["Secretory_minus_ATII"] = state_wide["Secretory"] - state_wide["ATII"]

state_test_rows = []
for gene, frame in state_wide.groupby("gene", sort=False):
    stat, pval = stats.wilcoxon(frame["Secretory"], frame["ATII"], alternative="two-sided")
    state_test_rows.append({
        "gene": gene, "n_paired_IPF_donors": len(frame),
        "median_ATII_log2CPM": frame["ATII"].median(),
        "median_Secretory_log2CPM": frame["Secretory"].median(),
        "median_Secretory_minus_ATII": frame["Secretory_minus_ATII"].median(),
        "wilcoxon_statistic": stat, "p_value": pval,
    })
state_tests = add_fdr(pd.DataFrame(state_test_rows)).sort_values("p_value")
state_pb.to_csv(OUT / "IPF_donor_state_pseudobulk.csv", index=False, encoding="utf-8-sig")
state_wide.to_csv(OUT / "IPF_paired_secretory_vs_ATII_values.csv", index=False, encoding="utf-8-sig")
state_tests.to_csv(OUT / "IPF_paired_secretory_vs_ATII_tests.csv", index=False, encoding="utf-8-sig")

# Descriptive IPF epithelial map across author-annotated states.
celltype_rows = []
for (donor, celltype), celltype_frame in ipf_meta.groupby(["Subject_Identity", "Manuscript_Identity"]):
    idx = celltype_frame.index.to_numpy()
    if celltype not in EPITHELIAL_TYPES or len(idx) < MIN_STATE_CELLS:
        continue
    total_umis = depth[idx].sum()
    sums = expr[:, idx].sum(axis=1)
    for gene in STATE_GENES:
        celltype_rows.append({
            "Subject_Identity": donor, "celltype": celltype, "gene": gene,
            "n_cells": len(idx), "log2CPM": np.log2(sums[gene_idx[gene]] / total_umis * 1e6 + 1),
        })
celltype_pb = pd.DataFrame(celltype_rows)
celltype_summary = (
    celltype_pb.groupby(["celltype", "gene"], as_index=False)
    .agg(n_donors=("Subject_Identity", "nunique"), median_log2CPM=("log2CPM", "median"))
)
celltype_pb.to_csv(OUT / "IPF_epithelial_celltype_pseudobulk.csv", index=False, encoding="utf-8-sig")
celltype_summary.to_csv(OUT / "IPF_epithelial_celltype_summary.csv", index=False, encoding="utf-8-sig")

# Secondary external validation of broad myeloid-to-epithelial LR support.
comm_rows = []
for donor, donor_meta in meta.groupby("Subject_Identity", sort=False):
    sender = donor_meta.index[donor_meta["CellType_Category"].eq("Myeloid")].to_numpy()
    receiver = donor_meta.index[donor_meta["CellType_Category"].eq("Epithelial")].to_numpy()
    if len(sender) < MIN_COMMUNICATION_CELLS or len(receiver) < MIN_COMMUNICATION_CELLS:
        continue
    diagnosis = donor_meta["Disease_Identity"].iloc[0]
    for ligand, receptor in FOCUS_PAIRS:
        ligand_values = expr[gene_idx[ligand], sender]
        receptor_values = expr[gene_idx[receptor], receiver]
        ligand_mean = np.log1p(ligand_values / depth[sender] * 1e4).mean()
        receptor_mean = np.log1p(receptor_values / depth[receiver] * 1e4).mean()
        comm_rows.append({
            "Subject_Identity": donor, "Diagnosis": diagnosis,
            "ligand": ligand, "receptor": receptor,
            "n_sender_cells": len(sender), "n_receiver_cells": len(receiver),
            "ligand_mean_log1pCP10K": ligand_mean,
            "receptor_mean_log1pCP10K": receptor_mean,
            "interaction_score": np.sqrt(ligand_mean * receptor_mean),
        })
communication = pd.DataFrame(comm_rows)
communication["interaction"] = communication["ligand"] + "–" + communication["receptor"]

comm_test_rows = []
for (ligand, receptor), frame in communication.groupby(["ligand", "receptor"], sort=False):
    ipf = frame.loc[frame["Diagnosis"].eq("IPF"), "interaction_score"].to_numpy()
    ctrl = frame.loc[frame["Diagnosis"].eq("Control"), "interaction_score"].to_numpy()
    u, pval = stats.mannwhitneyu(ipf, ctrl, alternative="two-sided")
    comm_test_rows.append({
        "ligand": ligand, "receptor": receptor,
        "n_IPF_donors": len(ipf), "n_Control_donors": len(ctrl),
        "median_IPF": np.median(ipf), "median_Control": np.median(ctrl),
        "IPF_minus_Control": np.median(ipf) - np.median(ctrl),
        "cliffs_delta": 2 * u / (len(ipf) * len(ctrl)) - 1,
        "p_value": pval,
    })
comm_tests = add_fdr(pd.DataFrame(comm_test_rows)).sort_values("p_value")
comm_tests["interaction"] = comm_tests["ligand"] + "–" + comm_tests["receptor"]
communication.to_csv(OUT / "donor_level_myeloid_epithelial_interactions.csv", index=False, encoding="utf-8-sig")
comm_tests.to_csv(OUT / "myeloid_epithelial_interaction_tests.csv", index=False, encoding="utf-8-sig")

sns.set_theme(style="whitegrid", context="paper", font_scale=1.05)

# Figure 1: paired external validation; each dot is an independent IPF donor.
order = state_tests.sort_values("median_Secretory_minus_ATII")["gene"].tolist()
plot_state = state_wide.copy()
plot_state["gene"] = pd.Categorical(plot_state["gene"], categories=order, ordered=True)
fig, ax = plt.subplots(figsize=(8.8, 5.4))
sns.stripplot(data=plot_state, x="Secretory_minus_ATII", y="gene", order=order,
              color="#4C78A8", alpha=.65, size=5, jitter=.17, ax=ax)
for ypos, gene in enumerate(order):
    row = state_tests.set_index("gene").loc[gene]
    ax.scatter(row["median_Secretory_minus_ATII"], ypos, marker="D", s=35,
               color="#D62728", zorder=5)
    ax.text(ax.get_xlim()[1], ypos, f" FDR={row['FDR_BH']:.2g}", va="center", fontsize=8)
ax.axvline(0, color="black", ls="--", lw=1)
ax.set_xlabel("Secretory (Club/Goblet) − ATII, log2CPM")
ax.set_ylabel("")
ax.set_title("GSE136831: paired epithelial-state validation in 15 IPF donors", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "Figure_1_external_paired_state_validation.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_1_external_paired_state_validation.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 2: descriptive external cell-state expression landscape.
heat = celltype_summary.pivot(index="celltype", columns="gene", values="median_log2CPM")
heat = heat.reindex(index=EPITHELIAL_TYPES, columns=STATE_GENES)
fig, ax = plt.subplots(figsize=(11, 5.8))
sns.heatmap(heat, cmap="viridis", annot=True, fmt=".1f", linewidths=.4,
            cbar_kws={"label": "Median donor log2CPM"}, ax=ax)
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title("GSE136831 IPF epithelial state landscape", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT / "Figure_2_external_epithelial_state_heatmap.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_2_external_epithelial_state_heatmap.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 3: independent replication of LR evidence.
comm_order = comm_tests.sort_values("cliffs_delta", ascending=False)["interaction"].tolist()
fig, ax = plt.subplots(figsize=(11, 5.5))
sns.boxplot(data=communication, x="interaction", y="interaction_score", hue="Diagnosis",
            order=comm_order, hue_order=["Control", "IPF"], palette=["#4C78A8", "#E45756"],
            showfliers=False, width=.65, ax=ax)
sns.stripplot(data=communication, x="interaction", y="interaction_score", hue="Diagnosis",
              order=comm_order, hue_order=["Control", "IPF"], dodge=True,
              palette=["#4C78A8", "#E45756"], alpha=.6, size=3.5, linewidth=.2, ax=ax)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:2], labels[:2], title="", frameon=False)
ax.set_xlabel(""); ax.set_ylabel("sqrt(mean ligand × mean receptor)")
ax.tick_params(axis="x", rotation=40)
ax.set_title("GSE136831: external myeloid → epithelial interaction support", fontweight="bold")
for xpos, interaction in enumerate(comm_order):
    row = comm_tests.loc[comm_tests["interaction"].eq(interaction)].iloc[0]
    ymax = communication.loc[communication["interaction"].eq(interaction), "interaction_score"].max()
    ax.text(xpos, ymax * 1.08 + .004, f"FDR={row['FDR_BH']:.2g}",
            ha="center", va="bottom", fontsize=7, rotation=45)
ax.margins(y=.17)
fig.tight_layout()
fig.savefig(OUT / "Figure_3_external_interaction_validation.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_3_external_interaction_validation.pdf", bbox_inches="tight")
plt.close(fig)

with open(OUT / "README.txt", "w", encoding="utf-8") as fh:
    fh.write(
        "Independent validation in GSE136831\n"
        "Cohort in metadata: 32 IPF, 28 Control, and 18 COPD subjects; COPD was excluded.\n"
        "Primary test: paired donor-level pseudobulk contrast between Club/Goblet secretory cells\n"
        "and ATII within IPF; >=5 cells per state were required, yielding 15 paired donors.\n"
        "Secondary LR test: broad Myeloid sender and Epithelial receiver; >=20 cells per compartment;\n"
        "two-sided Mann-Whitney tests across subjects and BH adjustment over eight prespecified pairs.\n"
        "These results independently reproduce a cell-state association and expression-supported\n"
        "communication signal; they do not establish lineage or causal signaling.\n"
    )

print("Paired IPF secretory-vs-ATII validation:")
print(state_tests[["gene", "n_paired_IPF_donors", "median_Secretory_minus_ATII", "p_value", "FDR_BH"]].to_string(index=False))
print("\nExternal myeloid-to-epithelial tests:")
print(comm_tests[["interaction", "n_IPF_donors", "n_Control_donors", "median_IPF", "median_Control",
                  "cliffs_delta", "p_value", "FDR_BH"]].to_string(index=False))
print(f"\nOutputs: {OUT}")
