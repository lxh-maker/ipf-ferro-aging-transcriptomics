#!/usr/bin/env python3
"""Independent fixed-module validation of IPF myeloid states in GSE136831.

The four gene modules are copied unchanged from the GSE135893 discovery
analysis.  The reduced expression panel is used only as a fixed-module state
projection; diffusion pseudotime is therefore descriptive and not evidence of
lineage.  Statistical summaries use patients as the unit of replication.
"""

from pathlib import Path
import json

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data" / "GSE136831"
OUT = ROOT / "06_mechanism_axis" / "myeloid_state_coupling" / "validation_GSE136831"
OUT.mkdir(parents=True, exist_ok=True)
SEED = 20261028
rng = np.random.default_rng(SEED)

GENES = [
    "FABP4", "INHBA", "MARCO", "MCEMP1", "LPL", "MRC1", "C1QA", "C1QB", "C1QC",
    "S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "TIMP1", "LYZ",
    "GPNMB", "LIPA", "FTL", "CTSD", "ASAH1", "HEXB", "LAMP1", "CYP27A1",
    "SPP1", "MMP9", "LGMN", "PLA2G7", "CHIT1", "APOC1", "TREM2", "CD9", "CHI3L1",
]
MODULES = {
    "Resident/alveolar": GENES[0:9],
    "S100A8 transitional": GENES[9:16],
    "Lysosome/lipid/iron": GENES[16:24],
    "SPP1 effector": GENES[24:33],
}
MODULE_ORDER = list(MODULES)
MYELOID_TYPES = ["cMonocyte", "ncMonocyte", "Macrophage_Alveolar", "Macrophage"]


def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    keep = np.isfinite(p)
    if keep.any():
        out[keep] = multipletests(p[keep], method="fdr_bh")[1]
    return out


meta = pd.read_csv(DATA / "GSE136831_AllCells.Samples.CellType.MetadataTable.txt.gz", sep="\t")
meta.columns = meta.columns.str.strip('"')
for col in meta.select_dtypes(include="object"):
    meta[col] = meta[col].str.strip('"')
assert len(meta) == 312928

raw = np.memmap(DATA / "GSE136831_myeloid_modules.uint32.bin", dtype=np.uint32,
                mode="r", shape=(len(GENES), len(meta)))
mask = (meta["Disease_Identity"] == "IPF") & meta["Manuscript_Identity"].isin(MYELOID_TYPES)
idx_all = np.flatnonzero(mask.to_numpy())
ipf = meta.iloc[idx_all].reset_index(drop=False).rename(columns={"index": "matrix_col"})
counts = np.asarray(raw[:, idx_all], dtype=np.float32).T
lib = ipf["nUMI"].to_numpy(dtype=np.float32)
expr = np.log1p(counts / np.maximum(lib[:, None], 1.0) * 1e4)

# Standardize genes within the independent IPF myeloid compartment, then apply
# the exact discovery modules without gene selection or weight refitting.
z = stats.zscore(expr, axis=0, nan_policy="omit")
z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
gene_to_i = {g: i for i, g in enumerate(GENES)}
for module, genes in MODULES.items():
    ipf[module] = z[:, [gene_to_i[g] for g in genes]].mean(axis=1)

# Author-annotated cell identities provide an independent, transparent anchor.
pseudobulk = ipf.groupby(["Subject_Identity", "Manuscript_Identity"], observed=True)[MODULE_ORDER].mean().reset_index()
pseudobulk.to_csv(OUT / "fixed_module_scores_by_patient_and_author_celltype.csv", index=False)

type_test_rows = []
for module in MODULE_ORDER:
    wide = pseudobulk.pivot(index="Subject_Identity", columns="Manuscript_Identity", values=module)
    complete = wide.dropna(subset=MYELOID_TYPES)
    if len(complete) >= 3:
        stat, p = stats.friedmanchisquare(*[complete[t] for t in MYELOID_TYPES])
    else:
        stat, p = np.nan, np.nan
    type_test_rows.append({"module": module, "n_complete_patients": len(complete),
                           "friedman_statistic": stat, "p_value": p})
type_tests = pd.DataFrame(type_test_rows)
type_tests["FDR_BH"] = bh(type_tests["p_value"])
type_tests.to_csv(OUT / "author_celltype_module_tests.csv", index=False)

# Up to 100 cells per patient and author myeloid identity avoids domination by
# large libraries while preserving the four reported identities.
sampled = []
for (_, _), dat in ipf.groupby(["Subject_Identity", "Manuscript_Identity"], observed=True):
    n = min(100, len(dat))
    sampled.extend(rng.choice(dat.index.to_numpy(), size=n, replace=False).tolist())
sampled = np.asarray(sorted(sampled), dtype=int)
smeta = ipf.loc[sampled].copy().reset_index(drop=True)
sx = expr[sampled]

adata = ad.AnnData(sx, obs=smeta.set_index("CellBarcode_Identity"), var=pd.DataFrame(index=GENES))
adata.layers["lognorm"] = adata.X.copy()
for module in MODULE_ORDER:
    adata.obs[module] = smeta[module].to_numpy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=20, svd_solver="arpack", random_state=SEED)
sc.external.pp.harmony_integrate(adata, key="Subject_Identity", basis="X_pca",
                                 adjusted_basis="X_pca_harmony", random_state=SEED,
                                 max_iter_harmony=30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=18, use_rep="X_pca_harmony", random_state=SEED)
sc.tl.umap(adata, min_dist=.35, random_state=SEED)
sc.tl.diffmap(adata, n_comps=15, random_state=SEED)

root_score = (adata.obs["Resident/alveolar"]
              - .5 * adata.obs["S100A8 transitional"]
              - .5 * adata.obs["SPP1 effector"]).to_numpy()
is_alv = adata.obs["Manuscript_Identity"].to_numpy() == "Macrophage_Alveolar"
candidate_pool = np.flatnonzero(is_alv)
cut = np.quantile(root_score[candidate_pool], .90)
candidates = candidate_pool[root_score[candidate_pool] >= cut]
space = adata.obsm["X_pca_harmony"][candidates, :15]
centroid = space.mean(axis=0)
root = int(candidates[np.argmin(((space - centroid) ** 2).sum(axis=1))])
adata.uns["iroot"] = root
sc.tl.dpt(adata, n_dcs=10)

cell = adata.obs[["Subject_Identity", "Manuscript_Identity", "dpt_pseudotime"] + MODULE_ORDER].copy()
cell["UMAP1"] = adata.obsm["X_umap"][:, 0]
cell["UMAP2"] = adata.obsm["X_umap"][:, 1]
cell.to_csv(OUT / "balanced_IPF_myeloid_fixed_module_projection.csv.gz", compression="gzip")

positions = []
correlations = []
for donor, dat in cell.groupby("Subject_Identity", observed=True):
    for module in MODULE_ORDER:
        threshold = dat[module].quantile(.75)
        top = dat.loc[dat[module] >= threshold, "dpt_pseudotime"]
        positions.append({"donor": donor, "module": module, "n_cells": len(dat),
                          "median_top_quartile_pseudotime": top.median()})
        rho, p = stats.spearmanr(dat["dpt_pseudotime"], dat[module])
        correlations.append({"donor": donor, "module": module, "rho": rho, "cell_level_p": p})
positions = pd.DataFrame(positions)
correlations = pd.DataFrame(correlations)
positions.to_csv(OUT / "fixed_module_positions_by_patient.csv", index=False)
correlations.to_csv(OUT / "fixed_module_pseudotime_correlations_by_patient.csv", index=False)

wide = positions.pivot(index="donor", columns="module", values="median_top_quartile_pseudotime")[MODULE_ORDER]
complete = wide.dropna()
friedman_stat, friedman_p = stats.friedmanchisquare(*[complete[m] for m in MODULE_ORDER])
pair_rows = []
for left, right in zip(MODULE_ORDER[:-1], MODULE_ORDER[1:]):
    d = complete[right] - complete[left]
    try:
        stat, p = stats.wilcoxon(d, zero_method="wilcox")
    except ValueError:
        stat, p = np.nan, 1.0
    pair_rows.append({"left": left, "right": right, "n_patients": len(d),
                      "median_right_minus_left": d.median(),
                      "n_patients_right_later": int((d > 0).sum()),
                      "wilcoxon_statistic": stat, "p_value": p})
pair_tests = pd.DataFrame(pair_rows)
pair_tests["FDR_BH"] = bh(pair_tests["p_value"])
pair_tests.to_csv(OUT / "adjacent_fixed_module_position_tests.csv", index=False)

full_order = ((complete[MODULE_ORDER[0]] < complete[MODULE_ORDER[1]]) &
              (complete[MODULE_ORDER[1]] < complete[MODULE_ORDER[2]]) &
              (complete[MODULE_ORDER[2]] < complete[MODULE_ORDER[3]]))

corr_summary = []
for module, dat in correlations.groupby("module", sort=False):
    vals = dat["rho"].dropna().to_numpy()
    try:
        stat, p = stats.wilcoxon(vals, zero_method="wilcox")
    except ValueError:
        stat, p = np.nan, 1.0
    corr_summary.append({"module": module, "n_patients": len(vals),
                         "median_patient_rho": np.median(vals), "p_value": p})
corr_summary = pd.DataFrame(corr_summary)
corr_summary["FDR_BH"] = bh(corr_summary["p_value"])
corr_summary.to_csv(OUT / "fixed_module_pseudotime_trend_tests.csv", index=False)

sns.set_theme(style="white", context="paper", font_scale=1.05)
palette = {"Resident/alveolar": "#4C78A8", "S100A8 transitional": "#E45756",
           "Lysosome/lipid/iron": "#F2A104", "SPP1 effector": "#7A5195"}

fig, axes = plt.subplots(2, 3, figsize=(14.5, 9.0))
p = axes[0, 0].scatter(cell.UMAP1, cell.UMAP2, c=cell.dpt_pseudotime,
                       cmap="viridis", s=5, linewidths=0)
axes[0, 0].scatter(cell.iloc[root].UMAP1, cell.iloc[root].UMAP2, marker="*", s=120,
                   c="white", edgecolors="black")
axes[0, 0].set_title("Fixed-module diffusion ordering")
fig.colorbar(p, ax=axes[0, 0], label="DPT")
for ax, module in zip(axes.flat[1:5], MODULE_ORDER):
    p = ax.scatter(cell.UMAP1, cell.UMAP2, c=cell[module], cmap="coolwarm",
                   vmin=-2, vmax=2, s=5, linewidths=0)
    ax.set_title(module)
    fig.colorbar(p, ax=ax, label="module z-score")
type_colors = dict(zip(MYELOID_TYPES, sns.color_palette("Set2", len(MYELOID_TYPES))))
for typ in MYELOID_TYPES:
    q = cell.Manuscript_Identity == typ
    axes.flat[5].scatter(cell.loc[q, "UMAP1"], cell.loc[q, "UMAP2"], s=5,
                         linewidths=0, color=type_colors[typ], label=typ)
axes.flat[5].set_title("Author annotation")
axes.flat[5].legend(frameon=False, fontsize=8, markerscale=2)
for ax in axes.flat:
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
fig.suptitle("GSE136831 independent fixed-module projection (IPF myeloid cells)", fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(OUT / "Figure_4_validation_fixed_module_projection.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_4_validation_fixed_module_projection.pdf", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
sns.boxplot(data=pseudobulk, x="Manuscript_Identity", y="value", hue="module",
            order=MYELOID_TYPES, hue_order=MODULE_ORDER, palette=palette,
            showfliers=False, ax=axes[0]) if False else None
plotdat = pseudobulk.melt(id_vars=["Subject_Identity", "Manuscript_Identity"],
                          value_vars=MODULE_ORDER, var_name="module", value_name="score")
sns.pointplot(data=plotdat, x="Manuscript_Identity", y="score", hue="module",
              order=MYELOID_TYPES, hue_order=MODULE_ORDER, palette=palette,
              estimator=np.median, errorbar=("pi", 50), dodge=.18, markers="o", ax=axes[0])
axes[0].axhline(0, color="grey", lw=.7)
axes[0].tick_params(axis="x", rotation=20)
axes[0].set_xlabel(""); axes[0].set_ylabel("Patient pseudobulk fixed-module score")
axes[0].set_title("Prespecified modules across author cell identities")
axes[0].legend(frameon=False, fontsize=8)

sns.stripplot(data=positions, x="module", y="median_top_quartile_pseudotime",
              order=MODULE_ORDER, hue="module", palette=palette, legend=False,
              jitter=.14, size=4, alpha=.65, ax=axes[1])
sns.boxplot(data=positions, x="module", y="median_top_quartile_pseudotime",
            order=MODULE_ORDER, color="white", showfliers=False, width=.5, ax=axes[1])
axes[1].tick_params(axis="x", rotation=25)
axes[1].set_xlabel(""); axes[1].set_ylabel("Median DPT of top-quartile cells")
axes[1].set_title(f"Patient-level state positions\nFriedman P={friedman_p:.3g}; complete order {full_order.sum()}/{len(full_order)}")
fig.tight_layout()
fig.savefig(OUT / "Figure_5_validation_patient_level_results.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_5_validation_patient_level_results.pdf", bbox_inches="tight")
plt.close(fig)

summary = {
    "dataset": "GSE136831", "n_IPF_myeloid_cells_all": int(len(ipf)),
    "n_balanced_projection_cells": int(len(cell)),
    "n_IPF_patients": int(ipf.Subject_Identity.nunique()),
    "friedman_state_position_statistic": float(friedman_stat),
    "friedman_state_position_p": float(friedman_p),
    "n_patients_complete_order": int(full_order.sum()),
    "n_patients_with_all_positions": int(len(full_order)),
    "root_cell": str(adata.obs_names[root]),
    "interpretation": "Independent fixed-module state projection, not lineage tracing.",
}
(OUT / "validation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print("Selected IPF myeloid cells", len(ipf), "patients", ipf.Subject_Identity.nunique())
print("Balanced projection cells", len(cell))
print("\nAdjacent state positions")
print(pair_tests.to_string(index=False))
print("\nPseudotime trends")
print(corr_summary.to_string(index=False))
print("\nSummary", summary)
print("Outputs", OUT)
