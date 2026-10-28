#!/usr/bin/env python3
"""Spatial validation of branched IPF myeloid programs in E-MTAB-14121.

Tests fixed resident, S100A8-transitional, lysosome/lipid/iron-processing and
SPP1-effector modules against fibrosis severity and the previously derived
LCN2-high/ACSL4-low epithelial transition score.  Replicates are collapsed
before donor-level inference. Spatial associations are not causal evidence.
"""

from pathlib import Path
import itertools
import json

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data" / "E-MTAB-14121_spatial"
BASE = ROOT / "06_mechanism_axis" / "spatial_ETAB14121"
OUT = ROOT / "06_mechanism_axis" / "myeloid_state_coupling" / "spatial_ETAB14121"
OUT.mkdir(parents=True, exist_ok=True)

MODULES = {
    "Resident/alveolar": ["FABP4", "INHBA", "MARCO", "MCEMP1", "LPL", "MRC1", "C1QA", "C1QB", "C1QC"],
    "S100A8 transitional": ["S100A8", "S100A9", "FCN1", "VCAN", "CTSS", "TIMP1", "LYZ"],
    "Lysosome/lipid/iron": ["GPNMB", "LIPA", "FTL", "CTSD", "ASAH1", "HEXB", "LAMP1", "CYP27A1"],
    "SPP1 effector": ["SPP1", "MMP9", "LGMN", "PLA2G7", "CHIT1", "APOC1", "TREM2", "CD9", "CHI3L1"],
}
MODULE_ORDER = list(MODULES)
K_VALUES = [7, 19, 37]  # index includes the focal spot; effective neighbors: 6, 18, 36


def decode(x):
    return np.array([v.decode() if isinstance(v, bytes) else str(v) for v in x])


def robust_z(x):
    x = np.asarray(x, float)
    med = np.nanmedian(x)
    scale = 1.4826 * np.nanmedian(np.abs(x - med))
    if not np.isfinite(scale) or scale < 1e-8:
        scale = np.nanstd(x)
    return np.zeros_like(x) if not np.isfinite(scale) or scale < 1e-8 else (x - med) / scale


def bh(p):
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    keep = np.isfinite(p)
    if keep.any():
        out[keep] = multipletests(p[keep], method="fdr_bh")[1]
    return out


def exact_sign_flip(values):
    v = np.asarray(values, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return np.nan
    observed = abs(v.mean())
    null = [abs(np.mean(v * np.asarray(s))) for s in itertools.product([-1, 1], repeat=len(v))]
    return float(np.mean(np.asarray(null) >= observed - 1e-12))


def partial_spearman(x, y, cov):
    good = np.isfinite(x) & np.isfinite(y) & np.isfinite(cov)
    if good.sum() < 30:
        return np.nan
    rx = stats.rankdata(np.asarray(x)[good]); ry = stats.rankdata(np.asarray(y)[good]); rc = stats.rankdata(np.asarray(cov)[good])
    design = np.column_stack([np.ones(len(rc)), rc])
    ex = rx - design @ np.linalg.lstsq(design, rx, rcond=None)[0]
    ey = ry - design @ np.linalg.lstsq(design, ry, rcond=None)[0]
    return float(stats.pearsonr(ex, ey).statistic)


def read_expr(path, requested):
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        mat = sparse.csc_matrix((g["data"][:], g["indices"][:], g["indptr"][:]), shape=tuple(g["shape"][:]))
        barcodes = decode(g["barcodes"][:])
        fg = g["features"]
        names = decode(fg["name"][:] if "name" in fg else fg["id"][:])
    total = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    total[total <= 0] = 1
    index = {}
    for i, gene in enumerate(names):
        if gene in requested:
            index.setdefault(gene, []).append(i)
    values = {}
    for gene in requested:
        rows = index.get(gene, [])
        if rows:
            cnt = np.asarray(mat[rows, :].sum(axis=0)).ravel().astype(float)
            values[gene] = np.log1p(cnt / total * 1e4)
        else:
            values[gene] = np.zeros(len(barcodes))
    return barcodes, values, set(index)


meta = pd.read_csv(DATA / "hs_visium_metadata.tsv", sep="\t")
meta["severity"] = meta["B_tissue_selection"].str.extract(r"B(\d)").astype(int)
meta["disease"] = np.where(meta["condition"].str.lower().eq("ipf"), "IPF", "Control")
base = pd.read_csv(BASE / "spot_scores_and_neighborhoods.csv.gz")
requested = sorted(set().union(*MODULES.values()))

spot_parts = []
section_rows = []
coupling_rows = []
coloc_rows = []
gene_rows = []

for _, mr in meta.iterrows():
    sample = mr["sample_id"]
    old = base[base.sample_id == sample].copy().reset_index(drop=True)
    barcodes, expr, present = read_expr(DATA / f"{sample}__filtered_feature_bc_matrix.h5", requested)
    fresh = pd.DataFrame({"barcode": barcodes})
    for module, genes in MODULES.items():
        used = [g for g in genes if g in present]
        # Equal gene weights after within-section robust standardization reduce
        # domination by abundant housekeeping-like transcripts such as FTL.
        fresh[module] = np.mean(np.vstack([robust_z(expr[g]) for g in used]), axis=0)
        gene_rows.append({"sample_id": sample, "module": module,
                          "genes_requested": ";".join(genes), "genes_used": ";".join(used),
                          "n_requested": len(genes), "n_used": len(used)})
    df = old.merge(fresh, on="barcode", how="left", validate="one_to_one")
    if df[MODULE_ORDER].isna().any().any():
        raise ValueError(f"module merge failed: {sample}")

    myeloid = df.myeloid_fraction.to_numpy(float)
    epi = df.epithelial_fraction.to_numpy(float)
    receiver = df.epithelial_receiver.astype(bool).to_numpy()
    coords = df[["array_row", "array_col"]].to_numpy(float)
    tree = cKDTree(coords)

    for module in MODULE_ORDER:
        sender = robust_z(df[module].to_numpy(float)) * myeloid
        df[f"sender_{module}"] = robust_z(sender)
        section_rows.append({
            "sample_id": sample, "donor": mr["subject_alias"], "disease": mr["disease"],
            "severity": int(mr["severity"]), "replicate": int(mr["replicate"]), "module": module,
            "myeloid_weighted_score": float(np.average(df[module], weights=myeloid + 1e-8)),
        })

    for k in K_VALUES:
        kk = min(k, len(df))
        _, ix = tree.query(coords, k=kk)
        neigh = ix[:, 1:] if kk > 1 else ix[:, :0]
        for module in MODULE_ORDER:
            exposure = df[f"sender_{module}"].to_numpy()[neigh].mean(axis=1)
            rho = stats.spearmanr(exposure[receiver], df.loc[receiver, "transition_score"]).statistic
            coupling_rows.append({
                "sample_id": sample, "donor": mr["subject_alias"], "disease": mr["disease"],
                "severity": int(mr["severity"]), "replicate": int(mr["replicate"]),
                "neighbors": kk - 1, "module": module, "rho_with_epithelial_transition": rho,
                "n_receiver_spots": int(receiver.sum()),
            })

    rich = myeloid >= max(.20, np.quantile(myeloid, .75))
    pairs = [("Resident/alveolar", "S100A8 transitional"),
             ("S100A8 transitional", "Lysosome/lipid/iron"),
             ("Lysosome/lipid/iron", "SPP1 effector"),
             ("S100A8 transitional", "SPP1 effector")]
    for left, right in pairs:
        rho = partial_spearman(df.loc[rich, left].to_numpy(), df.loc[rich, right].to_numpy(), myeloid[rich])
        coloc_rows.append({"sample_id": sample, "donor": mr["subject_alias"], "disease": mr["disease"],
                           "severity": int(mr["severity"]), "replicate": int(mr["replicate"]),
                           "left": left, "right": right, "partial_spearman": rho,
                           "n_myeloid_rich_spots": int(rich.sum())})

    keep = ["sample_id", "donor", "disease", "severity", "replicate", "barcode",
            "array_row", "array_col", "myeloid_fraction", "epithelial_fraction",
            "transition_score", "epithelial_receiver"] + MODULE_ORDER + [f"sender_{m}" for m in MODULE_ORDER]
    spot_parts.append(df[keep])

spots = pd.concat(spot_parts, ignore_index=True)
sections = pd.DataFrame(section_rows)
coupling = pd.DataFrame(coupling_rows)
coloc = pd.DataFrame(coloc_rows)
pd.DataFrame(gene_rows).to_csv(OUT / "spatial_module_gene_detection.csv", index=False)
spots.to_csv(OUT / "spatial_fixed_module_spot_scores.csv.gz", index=False, compression="gzip")
sections.to_csv(OUT / "section_fixed_module_scores.csv", index=False)
coupling.to_csv(OUT / "section_multiscale_epithelial_coupling.csv", index=False)
coloc.to_csv(OUT / "section_myeloid_module_colocalization.csv", index=False)

# Average technical replicates, then perform donor-level inference.
ds = sections.groupby(["donor", "disease", "severity", "module"], as_index=False)["myeloid_weighted_score"].mean()
dc = coupling.groupby(["donor", "disease", "severity", "neighbors", "module"], as_index=False)["rho_with_epithelial_transition"].mean()
dcol = coloc.groupby(["donor", "disease", "severity", "left", "right"], as_index=False)["partial_spearman"].mean()
ds.to_csv(OUT / "donor_severity_fixed_module_scores.csv", index=False)
dc.to_csv(OUT / "donor_severity_multiscale_coupling.csv", index=False)
dcol.to_csv(OUT / "donor_severity_module_colocalization.csv", index=False)

severity_tests = []
for module in MODULE_ORDER:
    vals = []
    sub = ds[(ds.disease == "IPF") & (ds.module == module)]
    for donor, dat in sub.groupby("donor"):
        if dat.severity.nunique() >= 3:
            vals.append((donor, stats.linregress(dat.severity, dat.myeloid_weighted_score).slope))
    x = np.asarray([v for _, v in vals])
    severity_tests.append({"module": module, "n_donors": len(x), "mean_slope": x.mean(),
                           "median_slope": np.median(x), "n_positive": int((x > 0).sum()),
                           "exact_sign_flip_p": exact_sign_flip(x),
                           "donor_slopes": ";".join(f"{d}:{v:.5g}" for d, v in vals)})
severity_tests = pd.DataFrame(severity_tests)
severity_tests["FDR_BH"] = bh(severity_tests.exact_sign_flip_p)
severity_tests.to_csv(OUT / "fixed_module_severity_trend_tests.csv", index=False)

coupling_tests = []
for (neighbors, module), sub in dc.groupby(["neighbors", "module"], observed=True):
    b3 = sub[(sub.disease == "IPF") & (sub.severity == 3)].rho_with_epithelial_transition.to_numpy()
    ctrl = sub[sub.disease == "Control"].rho_with_epithelial_transition.to_numpy()
    if len(b3) and len(ctrl):
        mw = stats.mannwhitneyu(b3, ctrl, method="exact")
        delta = np.mean(b3[:, None] > ctrl[None, :]) - np.mean(b3[:, None] < ctrl[None, :])
        coupling_tests.append({"neighbors": neighbors, "module": module, "n_IPF_B3": len(b3),
                               "n_control": len(ctrl), "median_rho_IPF_B3": np.median(b3),
                               "median_rho_control": np.median(ctrl), "cliffs_delta": delta,
                               "mann_whitney_exact_p": mw.pvalue, "B3_positive_donors": int((b3 > 0).sum())})
coupling_tests = pd.DataFrame(coupling_tests)
coupling_tests["FDR_BH"] = bh(coupling_tests.mann_whitney_exact_p)
coupling_tests.to_csv(OUT / "B3_vs_control_multiscale_coupling_tests.csv", index=False)

coloc_tests = []
for (left, right), sub in dcol.groupby(["left", "right"], observed=True):
    b3 = sub[(sub.disease == "IPF") & (sub.severity == 3)].partial_spearman.to_numpy()
    ctrl = sub[sub.disease == "Control"].partial_spearman.to_numpy()
    mw = stats.mannwhitneyu(b3, ctrl, method="exact")
    coloc_tests.append({"left": left, "right": right, "median_rho_IPF_B3": np.median(b3),
                        "median_rho_control": np.median(ctrl), "mann_whitney_exact_p": mw.pvalue,
                        "B3_positive_donors": int((b3 > 0).sum())})
coloc_tests = pd.DataFrame(coloc_tests)
coloc_tests["FDR_BH"] = bh(coloc_tests.mann_whitney_exact_p)
coloc_tests.to_csv(OUT / "B3_vs_control_module_colocalization_tests.csv", index=False)

sns.set_theme(style="whitegrid", context="paper", font_scale=1.1)
palette = {"Resident/alveolar": "#4C78A8", "S100A8 transitional": "#E45756",
           "Lysosome/lipid/iron": "#F2A104", "SPP1 effector": "#7A5195"}
fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.2))
ipf = ds[ds.disease == "IPF"]
for module in MODULE_ORDER:
    for donor, dat in ipf[ipf.module == module].groupby("donor"):
        dat = dat.sort_values("severity")
        axes[0].plot(dat.severity, dat.myeloid_weighted_score, color=palette[module], alpha=.32, lw=1)
    med = ipf[ipf.module == module].groupby("severity").myeloid_weighted_score.median()
    axes[0].plot(med.index, med.values, marker="o", color=palette[module], lw=2.7, label=module)
axes[0].set_xticks([1, 2, 3], ["B1 mild", "B2 moderate", "B3 severe"])
axes[0].set_ylabel("Myeloid-weighted spatial module score"); axes[0].set_xlabel("")
axes[0].set_title("Within-patient severity trajectories")
axes[0].legend(frameon=False, fontsize=8)

primary = dc[dc.neighbors == 18].copy()
primary["group"] = np.where(primary.disease == "Control", "Control", primary.severity.map({1: "IPF B1", 2: "IPF B2", 3: "IPF B3"}))
sns.pointplot(data=primary, x="group", y="rho_with_epithelial_transition", hue="module",
              order=["Control", "IPF B1", "IPF B2", "IPF B3"], hue_order=MODULE_ORDER,
              palette=palette, estimator=np.median, errorbar=("pi", 50), dodge=.25, ax=axes[1])
axes[1].axhline(0, color="black", lw=.8)
axes[1].set_xlabel(""); axes[1].set_ylabel("Neighbor module vs epithelial transition (rho)")
axes[1].set_title("Epithelial–myeloid spatial coupling (18 neighbors)")
axes[1].legend(frameon=False, fontsize=8)
fig.tight_layout()
fig.savefig(OUT / "Figure_6_spatial_branch_severity_and_coupling.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_6_spatial_branch_severity_and_coupling.pdf", bbox_inches="tight")
plt.close(fig)

summary = {"dataset": "E-MTAB-14121", "n_sections": int(meta.shape[0]),
           "n_IPF_donors": int(meta.loc[meta.disease == "IPF", "subject_alias"].nunique()),
           "n_control_donors": int(meta.loc[meta.disease == "Control", "subject_alias"].nunique()),
           "modules_fixed_before_spatial_analysis": MODULES,
           "interpretation": "Spatial association and deconvolution; no lineage or ligand-target causality."}
(OUT / "spatial_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print("Severity tests\n", severity_tests.to_string(index=False))
print("\n18-neighbor coupling\n", coupling_tests[coupling_tests.neighbors == 18].to_string(index=False))
print("\nColocalization\n", coloc_tests.to_string(index=False))
print("\nOutputs", OUT)
