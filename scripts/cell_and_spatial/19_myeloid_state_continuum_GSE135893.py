#!/usr/bin/env python3
"""Donor-aware IPF macrophage state-continuum analysis in GSE135893.

This analysis tests a prespecified four-module model:
resident/alveolar -> S100A8/S100A9 transitional -> lysosomal/lipid/ferritin
processing -> SPP1-rich effector. Pseudotime is interpreted as a state ordering,
not lineage tracing. All inferential summaries use biological donors.
"""

from pathlib import Path
import json
import warnings

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import io, stats
from statsmodels.nonparametric.smoothers_lowess import lowess
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
IN_DIR = ROOT / "03_virtual_knockout" / "inputs" / "IPF_Macrophages"
OUT = ROOT / "06_mechanism_axis" / "myeloid_state_coupling" / "discovery_GSE135893"
OUT.mkdir(parents=True, exist_ok=True)

SEED = 20261028
np.random.seed(SEED)
sc.settings.verbosity = 2

MODULES = {
    "Resident/alveolar": ["FABP4", "INHBA", "MARCO", "MCEMP1", "LPL", "MRC1", "C1QA", "C1QB", "C1QC"],
    "S100A8 transitional": ["S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "TIMP1", "LYZ"],
    "Lysosome/lipid/iron": ["GPNMB", "LIPA", "FTL", "CTSD", "ASAH1", "HEXB", "LAMP1", "CYP27A1"],
    "SPP1 effector": ["SPP1", "MMP9", "LGMN", "PLA2G7", "CHIT1", "APOC1", "TREM2", "CD9", "CHI3L1"],
}
MODULE_ORDER = list(MODULES)


def read_vector(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", header=None, dtype=str)[0].tolist()


def bh(values) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    keep = np.isfinite(values)
    if keep.any():
        out[keep] = multipletests(values[keep], method="fdr_bh")[1]
    return out


def dense_col(adata: ad.AnnData, gene: str) -> np.ndarray:
    arr = adata.raw[:, gene].X
    return np.asarray(arr.toarray() if hasattr(arr, "toarray") else arr).ravel()


genes = read_vector(IN_DIR / "genes.tsv")
cells = read_vector(IN_DIR / "cells.tsv")
meta = pd.read_csv(IN_DIR / "cell_metadata.csv").set_index("Cell").loc[cells].copy()
x = io.mmread(IN_DIR / "counts.mtx").tocsr().T
adata = ad.AnnData(x, obs=meta, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# Fixed module scores are mean standardized log-normalized expression. The same
# gene sets are carried unchanged into the independent cohort.
module_gene_report = []
for module, requested in MODULES.items():
    present = [g for g in requested if g in adata.var_names]
    missing = sorted(set(requested) - set(present))
    mat = np.column_stack([dense_col(adata, g) for g in present])
    z = stats.zscore(mat, axis=0, nan_policy="omit")
    z = np.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
    adata.obs[module] = np.clip(z, -5, 5).mean(axis=1)
    for gene in requested:
        module_gene_report.append({"module": module, "gene": gene, "detected": gene in present})
pd.DataFrame(module_gene_report).to_csv(OUT / "module_gene_definitions.csv", index=False)

# Manifold uses the broader 982-gene matrix and excludes technical families.
technical = adata.var_names.str.match(r"^(MT-|RPL|RPS)")
sc.pp.highly_variable_genes(adata, n_top_genes=min(600, int((~technical).sum())), flavor="seurat")
adata.var.loc[technical, "highly_variable"] = False
sc.pp.scale(adata, max_value=10)
n_pcs = min(30, int(adata.var["highly_variable"].sum()) - 1, adata.n_obs - 1)
sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack", random_state=SEED)
sc.external.pp.harmony_integrate(
    adata, key="Sample_Name", basis="X_pca", adjusted_basis="X_pca_harmony",
    random_state=SEED, max_iter_harmony=30,
)
sc.pp.neighbors(adata, n_neighbors=22, n_pcs=min(25, n_pcs), use_rep="X_pca_harmony", random_state=SEED)
sc.tl.umap(adata, min_dist=.35, spread=1.0, random_state=SEED)
sc.tl.diffmap(adata, n_comps=15, random_state=SEED)

# Root at the resident/alveolar end, avoiding a donor-specific expression extreme.
root_score = (
    adata.obs["Resident/alveolar"]
    - .5 * adata.obs["S100A8 transitional"]
    - .5 * adata.obs["SPP1 effector"]
).to_numpy()
cut = np.quantile(root_score, .90)
candidates = np.flatnonzero(root_score >= cut)
space = adata.obsm["X_pca_harmony"][candidates, :min(15, n_pcs)]
centroid = space.mean(axis=0)
root = int(candidates[np.argmin(((space - centroid) ** 2).sum(axis=1))])
adata.uns["iroot"] = root
sc.tl.dpt(adata, n_dcs=10)

# Module-dominant labels are descriptive only.
module_matrix = adata.obs[MODULE_ORDER].to_numpy()
adata.obs["dominant_module"] = pd.Categorical(
    np.asarray(MODULE_ORDER)[np.argmax(module_matrix, axis=1)], categories=MODULE_ORDER, ordered=True
)

cell_out = adata.obs[["Sample_Name", "orig.ident", "celltype", "dpt_pseudotime", "dominant_module"] + MODULE_ORDER].copy()
cell_out["UMAP1"] = adata.obsm["X_umap"][:, 0]
cell_out["UMAP2"] = adata.obsm["X_umap"][:, 1]
cell_out.to_csv(OUT / "IPF_macrophage_cell_state_scores.csv.gz", compression="gzip")
# HDF5 paths cannot contain forward slashes in DataFrame column names. Keep the
# publication-facing labels in CSV/figures and sanitize only the serialized copy.
adata_h5 = adata.copy()
adata_h5.obs.rename(columns={c: c.replace("/", "_") for c in adata_h5.obs.columns}, inplace=True)
adata_h5.write_h5ad(OUT / "IPF_macrophage_state_continuum.h5ad", compression="gzip")

# Donor-level module correlations with pseudotime.
donor_corr_rows = []
for donor, dat in cell_out.groupby("Sample_Name", observed=True):
    for module in MODULE_ORDER:
        rho, p = stats.spearmanr(dat["dpt_pseudotime"], dat[module])
        donor_corr_rows.append({"donor": donor, "module": module, "n_cells": len(dat), "rho": rho, "cell_level_p": p})
donor_corr = pd.DataFrame(donor_corr_rows)
donor_corr.to_csv(OUT / "module_pseudotime_donor_correlations.csv", index=False)

trend_rows = []
for module, dat in donor_corr.groupby("module", sort=False):
    rhos = dat["rho"].dropna().to_numpy()
    try:
        stat, p = stats.wilcoxon(rhos, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        stat, p = np.nan, 1.0
    trend_rows.append({
        "module": module, "n_donors": len(rhos), "median_donor_rho": np.median(rhos),
        "q25_rho": np.quantile(rhos, .25), "q75_rho": np.quantile(rhos, .75),
        "wilcoxon_statistic": stat, "p_value": p,
    })
module_trends = pd.DataFrame(trend_rows)
module_trends["FDR_BH"] = bh(module_trends["p_value"])
module_trends.to_csv(OUT / "module_pseudotime_trend_tests.csv", index=False)

# State order: within each donor, locate the cells in the top quartile for each
# fixed module and summarize their median pseudotime.
position_rows = []
for donor, dat in cell_out.groupby("Sample_Name", observed=True):
    for module in MODULE_ORDER:
        threshold = dat[module].quantile(.75)
        top = dat.loc[dat[module] >= threshold, "dpt_pseudotime"]
        position_rows.append({
            "donor": donor, "module": module, "n_top_quartile_cells": len(top),
            "median_top_quartile_pseudotime": top.median(),
        })
module_positions = pd.DataFrame(position_rows)
module_positions.to_csv(OUT / "module_state_positions_by_donor.csv", index=False)
position_wide = module_positions.pivot(index="donor", columns="module", values="median_top_quartile_pseudotime")[MODULE_ORDER]

friedman_stat, friedman_p = stats.friedmanchisquare(*[position_wide[m] for m in MODULE_ORDER])
pair_rows = []
for left, right in zip(MODULE_ORDER[:-1], MODULE_ORDER[1:]):
    diff = position_wide[right] - position_wide[left]
    try:
        stat, p = stats.wilcoxon(diff, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        stat, p = np.nan, 1.0
    pair_rows.append({
        "left": left, "right": right, "n_donors": len(diff),
        "median_right_minus_left": diff.median(),
        "n_donors_right_later": int((diff > 0).sum()),
        "wilcoxon_statistic": stat, "p_value": p,
    })
position_tests = pd.DataFrame(pair_rows)
position_tests["FDR_BH"] = bh(position_tests["p_value"])
position_tests.to_csv(OUT / "adjacent_module_position_tests.csv", index=False)

full_order = (
    (position_wide["Resident/alveolar"] < position_wide["S100A8 transitional"])
    & (position_wide["S100A8 transitional"] < position_wide["Lysosome/lipid/iron"])
    & (position_wide["Lysosome/lipid/iron"] < position_wide["SPP1 effector"])
)
ordering_summary = {
    "friedman_statistic": float(friedman_stat), "friedman_p": float(friedman_p),
    "n_donors_complete_order": int(full_order.sum()), "n_donors_total": int(len(full_order)),
    "complete_order_donors": position_wide.index[full_order].tolist(),
}
(OUT / "ordering_summary.json").write_text(json.dumps(ordering_summary, indent=2), encoding="utf-8")

# Gene-level donor correlations for transparent inspection.
all_module_genes = [g for geneset in MODULES.values() for g in geneset if g in adata.raw.var_names]
gene_expr = pd.DataFrame({g: dense_col(adata, g) for g in all_module_genes}, index=adata.obs_names)
gene_expr = pd.concat([adata.obs[["Sample_Name", "dpt_pseudotime"]], gene_expr], axis=1)
gene_rows = []
for gene in all_module_genes:
    module = next(m for m, gs in MODULES.items() if gene in gs)
    for donor, dat in gene_expr.groupby("Sample_Name", observed=True):
        if dat[gene].nunique() < 2:
            continue
        rho, p = stats.spearmanr(dat["dpt_pseudotime"], dat[gene])
        gene_rows.append({"gene": gene, "module": module, "donor": donor, "rho": rho, "cell_level_p": p})
gene_corr = pd.DataFrame(gene_rows)
gene_corr.to_csv(OUT / "gene_pseudotime_donor_correlations.csv", index=False)
gene_summary = []
for (module, gene), dat in gene_corr.groupby(["module", "gene"], sort=False):
    rhos = dat["rho"].dropna().to_numpy()
    try:
        stat, p = stats.wilcoxon(rhos, alternative="two-sided", zero_method="wilcox")
    except ValueError:
        stat, p = np.nan, 1.0
    gene_summary.append({"module": module, "gene": gene, "n_donors": len(rhos), "median_donor_rho": np.median(rhos), "p_value": p})
gene_summary = pd.DataFrame(gene_summary)
gene_summary["FDR_BH_global"] = bh(gene_summary["p_value"])
gene_summary.to_csv(OUT / "gene_pseudotime_trend_tests.csv", index=False)

# Pseudotime-bin trajectories, summarized per donor before aggregation.
cell_out["pseudotime_bin"] = pd.cut(cell_out["dpt_pseudotime"], bins=np.linspace(0, 1, 11), include_lowest=True, labels=False)
bin_donor = cell_out.groupby(["Sample_Name", "pseudotime_bin"], observed=True)[MODULE_ORDER].mean().reset_index()
bin_donor.to_csv(OUT / "module_scores_by_donor_pseudotime_bin.csv", index=False)
bin_summary = bin_donor.groupby("pseudotime_bin", observed=True)[MODULE_ORDER].agg(["median", lambda x: x.quantile(.25), lambda x: x.quantile(.75)])
bin_summary.to_csv(OUT / "module_scores_pseudotime_bin_summary.csv")

sns.set_theme(style="white", context="paper", font_scale=1.05)
palette = {
    "Resident/alveolar": "#4C78A8",
    "S100A8 transitional": "#E45756",
    "Lysosome/lipid/iron": "#F2A104",
    "SPP1 effector": "#7A5195",
}

# Figure 1: state manifold.
fig, axes = plt.subplots(2, 3, figsize=(14.5, 9.0))
pts = axes[0, 0].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                         c=adata.obs["dpt_pseudotime"], cmap="viridis", s=11, linewidths=0)
axes[0, 0].scatter(adata.obsm["X_umap"][root, 0], adata.obsm["X_umap"][root, 1], marker="*", s=120,
                   c="white", edgecolors="black", linewidths=.8)
axes[0, 0].set_title("Diffusion pseudotime")
fig.colorbar(pts, ax=axes[0, 0], label="DPT")
for ax, module in zip(axes.flat[1:5], MODULE_ORDER):
    p = ax.scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1], c=adata.obs[module],
                   cmap="coolwarm", s=11, linewidths=0, vmin=-2, vmax=2)
    ax.set_title(module)
    fig.colorbar(p, ax=ax, label="fixed module z-score")
ax = axes.flat[5]
for state in MODULE_ORDER:
    idx = adata.obs["dominant_module"].astype(str).to_numpy() == state
    ax.scatter(adata.obsm["X_umap"][idx, 0], adata.obsm["X_umap"][idx, 1], s=11, linewidths=0,
               c=palette[state], label=state, alpha=.8)
ax.set_title("Dominant module (descriptive)")
ax.legend(frameon=False, fontsize=8, loc="best")
for ax in axes.flat:
    ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
fig.suptitle("GSE135893 IPF macrophage state manifold (12 donor-balanced patients)", fontweight="bold", y=1.01)
fig.tight_layout()
fig.savefig(OUT / "Figure_1_discovery_myeloid_state_manifold.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_1_discovery_myeloid_state_manifold.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 2: donor-first trajectory summary and state positions.
fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.0))
for module in MODULE_ORDER:
    summ = bin_donor.groupby("pseudotime_bin", observed=True)[module].agg(
        median="median", q25=lambda x: x.quantile(.25), q75=lambda x: x.quantile(.75)
    ).reset_index()
    xbin = (summ["pseudotime_bin"].to_numpy() + .5) / 10
    axes[0].plot(xbin, summ["median"], marker="o", lw=2.2, color=palette[module], label=module)
    axes[0].fill_between(xbin, summ["q25"], summ["q75"], color=palette[module], alpha=.14)
axes[0].axhline(0, color="grey", lw=.7)
axes[0].set_xlabel("Diffusion pseudotime bin")
axes[0].set_ylabel("Median donor module z-score")
axes[0].set_title("Fixed modules across the state continuum")
axes[0].legend(frameon=False, fontsize=8)

sns.stripplot(data=module_positions, x="module", y="median_top_quartile_pseudotime",
              order=MODULE_ORDER, palette=palette, jitter=.15, size=5, alpha=.75, ax=axes[1])
sns.boxplot(data=module_positions, x="module", y="median_top_quartile_pseudotime",
            order=MODULE_ORDER, color="white", showfliers=False, width=.5, ax=axes[1])
axes[1].set_xlabel(""); axes[1].set_ylabel("Median pseudotime of top-quartile cells")
axes[1].tick_params(axis="x", rotation=25)
axes[1].set_title(f"Donor-level module positions\nFriedman P={friedman_p:.3g}; complete order {full_order.sum()}/{len(full_order)} donors")
fig.tight_layout()
fig.savefig(OUT / "Figure_2_discovery_module_trajectories.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_2_discovery_module_trajectories.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 3: every point is one donor.
fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.1))
sns.stripplot(data=donor_corr, x="rho", y="module", order=MODULE_ORDER, palette=palette,
              jitter=.16, size=5, alpha=.75, ax=axes[0])
for ypos, module in enumerate(MODULE_ORDER):
    row = module_trends.set_index("module").loc[module]
    axes[0].plot([row.q25_rho, row.q75_rho], [ypos, ypos], color="black", lw=2.2)
    axes[0].scatter(row.median_donor_rho, ypos, color="black", marker="D", s=32, zorder=5)
    axes[0].text(.98, ypos, f"FDR={row.FDR_BH:.3g}", transform=axes[0].get_yaxis_transform(),
                 ha="right", va="center", fontsize=8)
axes[0].axvline(0, color="grey", ls="--", lw=1)
axes[0].set_xlabel("Within-donor Spearman rho with pseudotime")
axes[0].set_ylabel("")
axes[0].set_title("Donor-consistent module trends")

heat = gene_summary.pivot(index="gene", columns="module", values="median_donor_rho")
ordered_genes = [g for m in MODULE_ORDER for g in MODULES[m] if g in heat.index]
heat = heat.reindex(ordered_genes)[MODULE_ORDER]
mask = heat.isna()
sns.heatmap(heat, mask=mask, cmap="vlag", center=0, vmin=-.7, vmax=.7, linewidths=.3,
            cbar_kws={"label": "Median donor rho"}, ax=axes[1])
axes[1].set_xlabel(""); axes[1].set_ylabel("")
axes[1].set_title("Gene-level state trends")
fig.tight_layout()
fig.savefig(OUT / "Figure_3_discovery_donor_consistency.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_3_discovery_donor_consistency.pdf", bbox_inches="tight")
plt.close(fig)

run_info = {
    "dataset": "GSE135893", "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
    "n_donors": int(adata.obs["Sample_Name"].nunique()), "cells_per_donor": adata.obs["Sample_Name"].value_counts().to_dict(),
    "root_cell": str(adata.obs_names[root]), "root_donor": str(adata.obs.iloc[root]["Sample_Name"]),
    "modules": MODULES, "ordering_summary": ordering_summary,
    "interpretation": "Computational macrophage state ordering; not lineage tracing.",
}
(OUT / "run_info.json").write_text(json.dumps(run_info, indent=2), encoding="utf-8")

print("Module trends")
print(module_trends.to_string(index=False))
print("\nAdjacent module positions")
print(position_tests.to_string(index=False))
print("\nOrdering summary", ordering_summary)
print("\nOutputs", OUT)
