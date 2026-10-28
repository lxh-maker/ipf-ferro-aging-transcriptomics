#!/usr/bin/env python3
"""Targeted epithelial-state pseudotime analysis for the IPF scRNA-seq cohort.

The input is the donor-balanced epithelial matrix prepared for scTenifoldKnk.
Pseudotime is interpreted as a computational state ordering, not lineage tracing.
Inference is summarized within each biological donor to avoid cell-level
pseudoreplication.
"""

from pathlib import Path
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
IN_DIR = ROOT / "03_virtual_knockout" / "inputs" / "IPF_Epithelial"
OUT_DIR = ROOT / "06_mechanism_axis" / "epithelial_pseudotime"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RANDOM_SEED = 20261027
np.random.seed(RANDOM_SEED)
sc.settings.verbosity = 2


def read_vector(path: Path) -> list[str]:
    return pd.read_csv(path, sep="\t", header=None, dtype=str)[0].tolist()


def bh(values: pd.Series) -> np.ndarray:
    out = np.full(len(values), np.nan)
    keep = values.notna().to_numpy()
    if keep.any():
        out[keep] = multipletests(values.to_numpy()[keep], method="fdr_bh")[1]
    return out


genes = read_vector(IN_DIR / "genes.tsv")
cells = read_vector(IN_DIR / "cells.tsv")
meta = pd.read_csv(IN_DIR / "cell_metadata.csv").set_index("Cell")
meta = meta.loc[cells].copy()

# Matrix Market file is genes x cells; AnnData expects cells x genes.
x = io.mmread(IN_DIR / "counts.mtx").tocsr().T
adata = ad.AnnData(x, obs=meta, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()
adata.obs_names_make_unique()

adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# The matrix already consists of variable/network genes. Remove technical gene
# families before selecting the most variable features for the state manifold.
technical = adata.var_names.str.match(r"^(MT-|RPL|RPS)")
sc.pp.highly_variable_genes(adata, n_top_genes=min(600, int((~technical).sum())), flavor="seurat")
adata.var.loc[technical, "highly_variable"] = False
if int(adata.var["highly_variable"].sum()) < 100:
    adata.var["highly_variable"] = ~technical

sc.pp.scale(adata, max_value=10)
n_pcs = min(30, int(adata.var["highly_variable"].sum()) - 1, adata.n_obs - 1)
sc.tl.pca(adata, n_comps=n_pcs, use_highly_variable=True, svd_solver="arpack", random_state=RANDOM_SEED)
# Multiple library preparations can otherwise form donor-specific branches. Harmony
# is used only to construct the state manifold; all expression values remain raw
# library-size-normalized observations and all inference remains donor-aware.
sc.external.pp.harmony_integrate(adata, key="Sample_Name", basis="X_pca", adjusted_basis="X_pca_harmony",
                                 random_state=RANDOM_SEED, max_iter_harmony=30)
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=min(25, n_pcs), use_rep="X_pca_harmony",
                random_state=RANDOM_SEED)
sc.tl.umap(adata, min_dist=0.35, spread=1.0, random_state=RANDOM_SEED)
sc.tl.diffmap(adata, n_comps=15, random_state=RANDOM_SEED)

# Root: the mature AT2 cell with the strongest alveolar program and the weakest
# injury/basal program. This makes the biological assumption explicit and auditable.
root_markers_pos = [g for g in ["SFTPC", "SFTPA1", "SFTPB", "ABCA3"] if g in adata.raw.var_names]
root_markers_neg = [g for g in ["KRT17", "KRT8", "VIM", "LCN2"] if g in adata.raw.var_names]
at2_idx = np.flatnonzero(adata.obs["celltype"].astype(str).to_numpy() == "AT2")
if len(at2_idx) == 0:
    raise RuntimeError("No AT2 cells are available for choosing a pseudotime root.")

def raw_gene(gene: str) -> np.ndarray:
    arr = adata.raw[:, gene].X
    return np.asarray(arr.toarray() if hasattr(arr, "toarray") else arr).ravel()

root_score = np.zeros(adata.n_obs)
for gene in root_markers_pos:
    root_score += stats.zscore(raw_gene(gene), nan_policy="omit")
for gene in root_markers_neg:
    root_score -= stats.zscore(raw_gene(gene), nan_policy="omit")
# Avoid choosing a donor-specific extreme: among the top quartile of mature-AT2
# scores, use the cell nearest their centroid in the donor-corrected embedding.
cutoff = np.nanquantile(root_score[at2_idx], .75)
root_candidates = at2_idx[root_score[at2_idx] >= cutoff]
root_space = adata.obsm["X_pca_harmony"][root_candidates, :min(15, n_pcs)]
root_centroid = root_space.mean(axis=0)
root_local = int(root_candidates[np.argmin(((root_space - root_centroid) ** 2).sum(axis=1))])
adata.uns["iroot"] = root_local
sc.tl.dpt(adata, n_dcs=10)

adata.write_h5ad(OUT_DIR / "IPF_epithelial_targeted_pseudotime.h5ad", compression="gzip")

cell_out = adata.obs[["Sample_Name", "orig.ident", "celltype", "dpt_pseudotime"]].copy()
cell_out["UMAP1"] = adata.obsm["X_umap"][:, 0]
cell_out["UMAP2"] = adata.obsm["X_umap"][:, 1]
cell_out.to_csv(OUT_DIR / "IPF_epithelial_cell_pseudotime.csv")

state_summary = (
    cell_out.groupby("celltype", observed=True)
    .agg(n_cells=("dpt_pseudotime", "size"), n_donors=("Sample_Name", "nunique"),
         median_pseudotime=("dpt_pseudotime", "median"),
         q25=("dpt_pseudotime", lambda z: z.quantile(.25)),
         q75=("dpt_pseudotime", lambda z: z.quantile(.75)))
    .sort_values("median_pseudotime")
)
state_summary.to_csv(OUT_DIR / "IPF_epithelial_state_pseudotime_summary.csv")

genes_of_interest = [
    g for g in ["SFTPC", "SFTPA1", "SCGB1A1", "SCGB3A1", "LCN2", "ACSL4",
                   "KRT8", "KRT17", "KRT19", "MUC5B", "VIM", "EPCAM"]
    if g in adata.raw.var_names
]
expression = pd.DataFrame(
    {g: raw_gene(g) for g in genes_of_interest},
    index=adata.obs_names,
)
expression = pd.concat([adata.obs[["Sample_Name", "celltype", "dpt_pseudotime"]], expression], axis=1)

# Donor-aware trend inference: calculate a rank correlation within each donor,
# then test whether donor correlations consistently differ from zero.
donor_rows = []
for gene in genes_of_interest:
    for donor, dat in expression.groupby("Sample_Name"):
        if len(dat) < 20 or dat[gene].nunique() < 2:
            continue
        rho, pval = stats.spearmanr(dat["dpt_pseudotime"], dat[gene])
        donor_rows.append({"gene": gene, "donor": donor, "n_cells": len(dat), "rho": rho, "cell_level_p": pval})
donor_corr = pd.DataFrame(donor_rows)
donor_corr.to_csv(OUT_DIR / "IPF_epithelial_pseudotime_donor_correlations.csv", index=False)

summary_rows = []
for gene, dat in donor_corr.groupby("gene"):
    rhos = dat["rho"].dropna().to_numpy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            stat, pval = stats.wilcoxon(rhos, alternative="two-sided", zero_method="wilcox")
        except ValueError:
            stat, pval = np.nan, 1.0
    summary_rows.append({
        "gene": gene, "n_donors": len(rhos), "median_donor_rho": np.median(rhos),
        "q25_donor_rho": np.quantile(rhos, .25), "q75_donor_rho": np.quantile(rhos, .75),
        "wilcoxon_statistic": stat, "p_value": pval,
    })
trend_summary = pd.DataFrame(summary_rows).sort_values("p_value")
trend_summary["FDR_BH"] = bh(trend_summary["p_value"])
trend_summary.to_csv(OUT_DIR / "IPF_epithelial_pseudotime_gene_trends.csv", index=False)

sns.set_theme(style="white", context="paper", font_scale=1.05)
palette = dict(zip(
    sorted(adata.obs["celltype"].astype(str).unique()),
    sns.color_palette("tab10", n_colors=adata.obs["celltype"].nunique()),
))

# Figure 1: state manifold and transparent root definition.
fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.7))
for state in sorted(palette):
    idx = adata.obs["celltype"].astype(str).to_numpy() == state
    axes[0].scatter(adata.obsm["X_umap"][idx, 0], adata.obsm["X_umap"][idx, 1], s=8, alpha=.75,
                    c=[palette[state]], label=state, linewidths=0)
axes[0].scatter(adata.obsm["X_umap"][root_local, 0], adata.obsm["X_umap"][root_local, 1],
                s=90, marker="*", c="black", edgecolors="white", linewidths=.6, label="AT2 root")
axes[0].set_title("IPF epithelial states")
axes[0].set_xlabel("UMAP1"); axes[0].set_ylabel("UMAP2")
axes[0].legend(frameon=False, fontsize=7, markerscale=1.7, bbox_to_anchor=(1.01, 1), loc="upper left")

pts = axes[1].scatter(adata.obsm["X_umap"][:, 0], adata.obsm["X_umap"][:, 1],
                      c=adata.obs["dpt_pseudotime"], cmap="viridis", s=9, linewidths=0)
axes[1].scatter(adata.obsm["X_umap"][root_local, 0], adata.obsm["X_umap"][root_local, 1],
                s=90, marker="*", c="white", edgecolors="black", linewidths=.7)
axes[1].set_title("Diffusion pseudotime")
axes[1].set_xlabel("UMAP1"); axes[1].set_ylabel("UMAP2")
fig.colorbar(pts, ax=axes[1], label="DPT")

order = state_summary.index.tolist()
sns.violinplot(data=cell_out, x="celltype", y="dpt_pseudotime", order=order,
               palette=palette, inner="quartile", cut=0, linewidth=.7, ax=axes[2])
axes[2].set_title("State ordering (descriptive)")
axes[2].set_xlabel(""); axes[2].set_ylabel("Diffusion pseudotime")
axes[2].tick_params(axis="x", rotation=55, labelsize=7)
fig.suptitle("Targeted pseudotime of donor-balanced IPF epithelial cells", y=1.03, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT_DIR / "Figure_1_IPF_epithelial_pseudotime.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT_DIR / "Figure_1_IPF_epithelial_pseudotime.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 2: LOWESS trajectories. Curves show expression-state relationships;
# donor-level consistency is annotated from the inferential table.
plot_genes = [g for g in ["SFTPC", "SCGB1A1", "SCGB3A1", "LCN2", "ACSL4", "KRT17"] if g in genes_of_interest]
fig, axes = plt.subplots(2, 3, figsize=(13, 7.2), sharex=True)
for ax, gene in zip(axes.ravel(), plot_genes):
    dat = expression[["dpt_pseudotime", gene]].sort_values("dpt_pseudotime")
    if len(dat) > 500:
        shown = dat.sample(500, random_state=RANDOM_SEED)
    else:
        shown = dat
    ax.scatter(shown["dpt_pseudotime"], shown[gene], s=5, alpha=.12, color="#4C78A8", linewidths=0)
    smooth = lowess(dat[gene], dat["dpt_pseudotime"], frac=.25, it=1, return_sorted=True)
    ax.plot(smooth[:, 0], smooth[:, 1], color="#D62728", lw=2.2)
    row = trend_summary.loc[trend_summary["gene"] == gene].iloc[0]
    ax.set_title(f"{gene}  median donor rho={row['median_donor_rho']:.2f}\nFDR={row['FDR_BH']:.3g}")
    ax.set_ylabel("log1p normalized expression")
    ax.set_xlabel("Diffusion pseudotime")
fig.suptitle("Expression remodeling along the IPF epithelial state axis", y=1.01, fontweight="bold")
fig.tight_layout()
fig.savefig(OUT_DIR / "Figure_2_IPF_epithelial_gene_trajectories.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT_DIR / "Figure_2_IPF_epithelial_gene_trajectories.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 3: each point is a biological donor.
ordered_genes = trend_summary.sort_values("median_donor_rho")["gene"].tolist()
fig, ax = plt.subplots(figsize=(8.2, max(4.5, .43 * len(ordered_genes))))
sns.stripplot(data=donor_corr, x="rho", y="gene", order=ordered_genes, color="#4C78A8",
              alpha=.65, size=4.5, jitter=.17, ax=ax)
for y, gene in enumerate(ordered_genes):
    row = trend_summary.set_index("gene").loc[gene]
    ax.plot([row["q25_donor_rho"], row["q75_donor_rho"]], [y, y], color="black", lw=2)
    ax.scatter(row["median_donor_rho"], y, color="#D62728", marker="D", s=26, zorder=5)
ax.axvline(0, color="grey", ls="--", lw=1)
ax.set_xlabel("Within-donor Spearman rho with pseudotime")
ax.set_ylabel("")
ax.set_title("Donor-consistent epithelial gene trends (12 IPF patients)", fontweight="bold")
fig.tight_layout()
fig.savefig(OUT_DIR / "Figure_3_IPF_epithelial_donor_trends.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT_DIR / "Figure_3_IPF_epithelial_donor_trends.pdf", bbox_inches="tight")
plt.close(fig)

with open(OUT_DIR / "README.txt", "w", encoding="utf-8") as fh:
    fh.write(
        "IPF epithelial targeted pseudotime\n"
        f"Cells: {adata.n_obs}; genes: {adata.n_vars}; biological donors: {adata.obs['Sample_Name'].nunique()}\n"
        f"Root cell: {adata.obs_names[root_local]} (author annotation: {adata.obs.iloc[root_local]['celltype']})\n"
        "State manifold: PCA followed by Harmony integration across Sample_Name; expression is not batch-adjusted.\n"
        "Root criterion: centroid-nearest cell among the top-quartile mature AT2 scores (high SFTPC/SFTPA1/\n"
        "SFTPB/ABCA3 and low KRT17/KRT8/VIM/LCN2).\n"
        "Statistics: Spearman correlations are computed within donors; Wilcoxon signed-rank tests assess\n"
        "consistency across donors and are BH-adjusted across genes. Pseudotime is a computational\n"
        "state ordering and must not be described as direct lineage evidence.\n"
    )

print("State order:")
print(state_summary[["n_cells", "n_donors", "median_pseudotime"]].to_string())
print("\nDonor-aware gene trends:")
print(trend_summary[["gene", "n_donors", "median_donor_rho", "p_value", "FDR_BH"]].to_string(index=False))
print(f"\nOutputs: {OUT_DIR}")
