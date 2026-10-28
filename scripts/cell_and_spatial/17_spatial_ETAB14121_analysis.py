from __future__ import annotations

import hashlib
import itertools
import json
import math
import os
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse, stats
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "00_data" / "E-MTAB-14121_spatial"
OUT_DIR = ROOT / "06_mechanism_axis" / "spatial_ETAB14121"
C2L_DIR = DATA_DIR / "cell2location_habermann2020" / "cell2location_habermann2020"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20261027)
N_PERM = 200

MODULES = {
    "myeloid_ferroaging": [
        "S100A8", "S100A9", "LIPA", "FTL", "FTH1", "ASAH1", "CTSD",
        "HEXB", "LAMP1", "CYP27A1", "GPNMB", "CD68", "PLA2G7",
    ],
    "myeloid_identity": [
        "LST1", "FCER1G", "TYROBP", "CTSS", "AIF1", "LILRB1", "CD68",
    ],
    "secretory_remodeling": [
        "LCN2", "SCGB1A1", "SCGB3A1", "MUC5B", "KRT17", "KRT19",
        "CLDN1", "KRT8",
    ],
    "alveolar_acsl4": [
        "ACSL4", "SFTPC", "SFTPA1", "SFTPA2", "ABCA3", "NAPSA", "AGER", "CAV1",
    ],
    "epithelial_identity": ["EPCAM", "KRT8", "KRT18", "KRT19", "CLDN4", "KRT17"],
    "fibrotic_stroma": ["COL1A1", "COL1A2", "COL3A1", "CTHRC1", "COL14A1", "POSTN"],
}

FOCUS_GENES = [
    "TGFB1", "S100A8", "S100A9", "SPP1", "LCN2", "ACSL4", "PTGS2",
    "EGR1", "HIF1A", "CD44", "TGFBR1", "TGFBR2", "ITGB6", "ITGB8",
]


def decode_array(x):
    return np.array([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in x])


def robust_z(x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med)) * 1.4826
    if not np.isfinite(mad) or mad < 1e-8:
        mad = np.nanstd(x)
    if not np.isfinite(mad) or mad < 1e-8:
        return np.zeros_like(x)
    return (x - med) / mad


def weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    good = np.isfinite(x) & np.isfinite(w) & (w >= 0)
    if not np.any(good) or np.sum(w[good]) <= 0:
        return float("nan")
    return float(np.average(x[good], weights=w[good] + 1e-8))


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    good = np.isfinite(x) & np.isfinite(y)
    if np.sum(good) < 20 or np.nanstd(x[good]) < 1e-10 or np.nanstd(y[good]) < 1e-10:
        return float("nan")
    return float(stats.spearmanr(x[good], y[good]).statistic)


def bh_adjust(pvalues):
    p = np.asarray(pvalues, dtype=float)
    out = np.full(p.shape, np.nan)
    good = np.isfinite(p)
    if not np.any(good):
        return out
    vals = p[good]
    order = np.argsort(vals)
    ranked = vals[order] * len(vals) / np.arange(1, len(vals) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adj = np.empty_like(ranked)
    adj[order] = np.minimum(ranked, 1.0)
    out[good] = adj
    return out


def exact_signed_mean_p(slopes: np.ndarray) -> float:
    slopes = np.asarray(slopes, dtype=float)
    slopes = slopes[np.isfinite(slopes)]
    if len(slopes) == 0:
        return float("nan")
    obs = abs(np.mean(slopes))
    null = []
    for signs in itertools.product([-1, 1], repeat=len(slopes)):
        null.append(abs(np.mean(slopes * np.asarray(signs))))
    return float(np.mean(np.asarray(null) >= obs - 1e-12))


def read_10x_h5(path: Path, genes_needed: list[str]):
    with h5py.File(path, "r") as h:
        g = h["matrix"]
        shape = tuple(int(v) for v in g["shape"][:])
        mat = sparse.csc_matrix(
            (g["data"][:], g["indices"][:], g["indptr"][:]), shape=shape
        )
        barcodes = decode_array(g["barcodes"][:])
        feature_group = g["features"]
        if "name" in feature_group:
            feature_names = decode_array(feature_group["name"][:])
        else:
            feature_names = decode_array(feature_group["id"][:])

    total = np.asarray(mat.sum(axis=0)).ravel().astype(float)
    total[total <= 0] = 1.0
    gene_to_rows: dict[str, list[int]] = {}
    needed = set(genes_needed)
    for i, gene in enumerate(feature_names):
        if gene in needed:
            gene_to_rows.setdefault(gene, []).append(i)

    expr = {}
    for gene in genes_needed:
        rows = gene_to_rows.get(gene, [])
        if rows:
            counts = np.asarray(mat[rows, :].sum(axis=0)).ravel().astype(float)
            expr[gene] = np.log1p(counts / total * 10000.0)
        else:
            expr[gene] = np.zeros(mat.shape[1], dtype=float)
    return barcodes, total, expr, sorted(gene_to_rows)


def module_mean(expr: dict[str, np.ndarray], genes: list[str], present: set[str]):
    used = [g for g in genes if g in present]
    if not used:
        return np.zeros_like(next(iter(expr.values()))), used
    return np.mean(np.vstack([expr[g] for g in used]), axis=0), used


meta = pd.read_csv(DATA_DIR / "hs_visium_metadata.tsv", sep="\t")
meta["severity"] = meta["B_tissue_selection"].str.extract(r"B(\d)").astype(int)
meta["disease"] = np.where(meta["condition"].str.lower().eq("ipf"), "IPF", "Control")
meta.to_csv(OUT_DIR / "sample_metadata.csv", index=False)

all_needed = sorted(set(FOCUS_GENES).union(*MODULES.values()))
all_spots = []
section_rows = []
neighbor_rows = []
module_gene_rows = []

# Visium barcodes have fixed array coordinates. BioStudies mistakenly serves the
# V10T03-280-D1 positions path as an HDF5 file, so retain a barcode-to-array map
# from a valid section as a documented fallback for neighborhood construction.
master_pos_path = DATA_DIR / "V10T03-279-B1_spatial_tissue_positions_list.csv"
master_pos = pd.read_csv(
    master_pos_path,
    header=None,
    names=["barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"],
)[["barcode", "array_row", "array_col"]].drop_duplicates("barcode")

for _, mrow in meta.iterrows():
    sample = mrow["sample_id"]
    h5_path = DATA_DIR / f"{sample}__filtered_feature_bc_matrix.h5"
    pos_path = DATA_DIR / f"{sample}_spatial_tissue_positions_list.csv"
    c2l_path = C2L_DIR / f"{sample}_spot_cell_abundances_5pc.csv"
    if not (h5_path.exists() and pos_path.exists() and c2l_path.exists()):
        raise FileNotFoundError(f"Missing input for {sample}")

    barcodes, total_counts, expr, present_genes = read_10x_h5(h5_path, all_needed)
    present = set(present_genes)

    position_source = "sample_specific"
    try:
        pos = pd.read_csv(
            pos_path,
            header=None,
            names=["barcode", "in_tissue", "array_row", "array_col", "pixel_row", "pixel_col"],
        )
        if pos.shape[1] != 6 or not pd.api.types.is_numeric_dtype(pos["array_row"]):
            raise ValueError("not a valid positions CSV")
    except Exception:
        pos = master_pos.copy()
        pos["in_tissue"] = np.nan
        pos["pixel_row"] = np.nan
        pos["pixel_col"] = np.nan
        position_source = "fixed_visium_array_fallback"

    c2l = pd.read_csv(c2l_path)
    prefix = f"{sample}_"
    c2l["barcode"] = c2l["spot_id"].astype(str).str.replace(
        f"^{prefix}", "", regex=True
    )
    c2l_cols = [c for c in c2l.columns if c.startswith("q05cell_abundance_w_sf_")]
    rename = {c: c.replace("q05cell_abundance_w_sf_", "") for c in c2l_cols}
    c2l = c2l[["barcode"] + c2l_cols].rename(columns=rename)

    df = pd.DataFrame({"barcode": barcodes, "total_counts": total_counts})
    df = df.merge(pos, on="barcode", how="left", validate="one_to_one")
    df = df.merge(c2l, on="barcode", how="left", validate="one_to_one")
    if df["array_row"].isna().any():
        raise ValueError(f"Spatial coordinates failed to merge for {sample}")
    c2l_names = list(rename.values())
    if df[c2l_names].isna().any().any():
        raise ValueError(f"cell2location abundance failed to merge for {sample}")

    for gene in FOCUS_GENES:
        df[gene] = expr[gene]

    for module, genes in MODULES.items():
        score, used = module_mean(expr, genes, present)
        df[module] = score
        module_gene_rows.append(
            {
                "sample_id": sample,
                "module": module,
                "genes_requested": ";".join(genes),
                "genes_used": ";".join(used),
                "n_requested": len(genes),
                "n_used": len(used),
            }
        )

    def csum(names):
        cols = [c for c in names if c in df.columns]
        return df[cols].sum(axis=1).to_numpy(float)

    myeloid_ab = csum(["Macrophages", "Monocytes", "Proliferating Macrophages", "cDCs", "pDCs"])
    epithelial_ab = csum([
        "AT1", "AT2", "Basal", "Ciliated", "Differentiating Ciliated", "KRT5-/KRT17+",
        "MUC5AC+ High", "MUC5B+", "Proliferating Epithelial Cells", "SCGB3A2+",
        "SCGB3A2+ SCGB1A1+", "Transitional AT2",
    ])
    abnormal_epi_ab = csum([
        "KRT5-/KRT17+", "MUC5AC+ High", "MUC5B+", "SCGB3A2+",
        "SCGB3A2+ SCGB1A1+", "Transitional AT2",
    ])
    at2_ab = csum(["AT2"])
    total_ab = df[c2l_names].sum(axis=1).to_numpy(float)
    total_ab[total_ab <= 0] = 1.0
    df["myeloid_fraction"] = myeloid_ab / total_ab
    df["epithelial_fraction"] = epithelial_ab / total_ab
    df["abnormal_epithelial_fraction"] = abnormal_epi_ab / total_ab
    df["AT2_fraction"] = at2_ab / total_ab

    df["secretory_z"] = robust_z(df["secretory_remodeling"].to_numpy())
    df["alveolar_z"] = robust_z(df["alveolar_acsl4"].to_numpy())
    df["transition_score"] = 0.5 * (df["secretory_z"] - df["alveolar_z"])
    df["LCN2_ACSL4_balance"] = 0.5 * (
        robust_z(df["LCN2"].to_numpy()) - robust_z(df["ACSL4"].to_numpy())
    )

    tgfb1_sender_raw = df["TGFB1"].to_numpy() * df["myeloid_fraction"].to_numpy()
    s100a8_sender_raw = df["myeloid_ferroaging"].to_numpy() * df["myeloid_fraction"].to_numpy()
    spp1_sender_raw = df["SPP1"].to_numpy() * df["myeloid_fraction"].to_numpy()
    df["TGFB1_myeloid_sender"] = robust_z(tgfb1_sender_raw)
    df["S100A8_ferro_sender"] = robust_z(s100a8_sender_raw)
    df["SPP1_myeloid_sender"] = robust_z(spp1_sender_raw)

    coords = df[["array_row", "array_col"]].to_numpy(float)
    k = min(7, len(df))
    _, idx = cKDTree(coords).query(coords, k=k)
    neighbors = idx[:, 1:] if k > 1 else idx[:, :0]
    exposure_columns = {
        "TGFB1": "TGFB1_myeloid_sender",
        "S100A8_ferro": "S100A8_ferro_sender",
        "SPP1": "SPP1_myeloid_sender",
    }
    for label, col in exposure_columns.items():
        vals = df[col].to_numpy(float)
        df[f"neighbor_{label}"] = vals[neighbors].mean(axis=1)

    receiver = (
        (df["epithelial_fraction"].to_numpy() >= 0.20)
        & (df["epithelial_fraction"].to_numpy() > df["myeloid_fraction"].to_numpy())
    )
    if np.sum(receiver) < 100:
        receiver = df["epithelial_fraction"].to_numpy() >= np.quantile(
            df["epithelial_fraction"].to_numpy(), 0.60
        )
    df["epithelial_receiver"] = receiver

    responses = {
        "secretory_transition": "transition_score",
        "LCN2_minus_ACSL4": "LCN2_ACSL4_balance",
        "secretory_program": "secretory_z",
        "alveolar_ACSL4_program": "alveolar_z",
    }
    for exposure, sender_col in exposure_columns.items():
        observed_exposure = df[f"neighbor_{exposure}"].to_numpy(float)
        sender = df[sender_col].to_numpy(float)
        for response, response_col in responses.items():
            y = df[response_col].to_numpy(float)
            obs = safe_spearman(observed_exposure[receiver], y[receiver])
            if response in {"secretory_transition", "LCN2_minus_ACSL4"}:
                null = np.empty(N_PERM, dtype=float)
                for b in range(N_PERM):
                    perm_sender = RNG.permutation(sender)
                    perm_exposure = perm_sender[neighbors].mean(axis=1)
                    null[b] = safe_spearman(perm_exposure[receiver], y[receiver])
                p = (1 + np.sum(np.abs(null) >= abs(obs))) / (N_PERM + 1)
                null_mean = float(np.nanmean(null))
                null_sd = float(np.nanstd(null))
            else:
                p = np.nan
                null_mean = np.nan
                null_sd = np.nan
            neighbor_rows.append({
                "sample_id": sample,
                "donor": mrow["subject_alias"],
                "disease": mrow["disease"],
                "severity": int(mrow["severity"]),
                "replicate": int(mrow["replicate"]),
                "exposure": exposure,
                "response": response,
                "n_receiver_spots": int(np.sum(receiver)),
                "spearman_rho": obs,
                "permutation_p": p,
                "null_mean": null_mean,
                "null_sd": null_sd,
            })

    epi_w = df["epithelial_fraction"].to_numpy(float)
    my_w = df["myeloid_fraction"].to_numpy(float)
    section_rows.append({
        "sample_id": sample,
        "donor": mrow["subject_alias"],
        "disease": mrow["disease"],
        "severity": int(mrow["severity"]),
        "replicate": int(mrow["replicate"]),
        "n_spots": len(df),
        "n_epithelial_receivers": int(np.sum(receiver)),
        "position_source": position_source,
        "has_image_coordinates": bool(df["pixel_row"].notna().all()),
        "mean_epithelial_fraction": float(np.mean(epi_w)),
        "mean_myeloid_fraction": float(np.mean(my_w)),
        "epithelial_secretory_score": weighted_mean(df["secretory_remodeling"].to_numpy(), epi_w),
        "epithelial_alveolar_ACSL4_score": weighted_mean(df["alveolar_acsl4"].to_numpy(), epi_w),
        "epithelial_transition_score": weighted_mean(df["transition_score"].to_numpy(), epi_w),
        "epithelial_LCN2": weighted_mean(df["LCN2"].to_numpy(), epi_w),
        "epithelial_ACSL4": weighted_mean(df["ACSL4"].to_numpy(), epi_w),
        "myeloid_ferroaging_score": weighted_mean(df["myeloid_ferroaging"].to_numpy(), my_w),
        "myeloid_TGFB1": weighted_mean(df["TGFB1"].to_numpy(), my_w),
        "myeloid_S100A8": weighted_mean(df["S100A8"].to_numpy(), my_w),
        "myeloid_SPP1": weighted_mean(df["SPP1"].to_numpy(), my_w),
        "fibrotic_stroma_score": float(np.mean(df["fibrotic_stroma"])),
    })

    df.insert(0, "sample_id", sample)
    df.insert(1, "donor", mrow["subject_alias"])
    df.insert(2, "disease", mrow["disease"])
    df.insert(3, "severity", int(mrow["severity"]))
    df.insert(4, "replicate", int(mrow["replicate"]))
    df["position_source"] = position_source
    keep = [
        "sample_id", "donor", "disease", "severity", "replicate", "barcode", "position_source",
        "array_row", "array_col", "pixel_row", "pixel_col", "total_counts",
        "myeloid_fraction", "epithelial_fraction", "abnormal_epithelial_fraction", "AT2_fraction",
        "epithelial_receiver", "myeloid_ferroaging", "secretory_remodeling",
        "alveolar_acsl4", "fibrotic_stroma", "transition_score", "LCN2_ACSL4_balance",
        "TGFB1_myeloid_sender", "S100A8_ferro_sender", "SPP1_myeloid_sender",
        "neighbor_TGFB1", "neighbor_S100A8_ferro", "neighbor_SPP1",
    ] + FOCUS_GENES
    all_spots.append(df[keep])

spots = pd.concat(all_spots, ignore_index=True)
sections = pd.DataFrame(section_rows)
neighbors_df = pd.DataFrame(neighbor_rows)
neighbors_df["permutation_FDR"] = bh_adjust(neighbors_df["permutation_p"])
module_genes = pd.DataFrame(module_gene_rows)

spots.to_csv(OUT_DIR / "spot_scores_and_neighborhoods.csv.gz", index=False, compression="gzip")
sections.to_csv(OUT_DIR / "section_level_spatial_scores.csv", index=False)
neighbors_df.to_csv(OUT_DIR / "section_level_neighbor_correlations.csv", index=False)
module_genes.to_csv(OUT_DIR / "module_genes_detected_by_section.csv", index=False)

# Collapse technical replicates before any donor-level inference.
metric_cols = [c for c in sections.columns if c not in {
    "sample_id", "donor", "disease", "severity", "replicate", "n_spots", "n_epithelial_receivers",
    "position_source", "has_image_coordinates"
}]
donor_severity = sections.groupby(["donor", "disease", "severity"], as_index=False)[metric_cols].mean()
donor_severity.to_csv(OUT_DIR / "donor_severity_spatial_scores.csv", index=False)

neighbor_donor_severity = (
    neighbors_df.groupby(["donor", "disease", "severity", "exposure", "response"], as_index=False)
    ["spearman_rho"].mean()
)
neighbor_donor_severity.to_csv(OUT_DIR / "donor_severity_neighbor_correlations.csv", index=False)

# Exact donor-sign-flip tests for the B1->B3 trend (n=4 donors).
trend_rows = []
ipf_ds = donor_severity[donor_severity["disease"] == "IPF"]
for metric in metric_cols:
    slopes = []
    for donor, sub in ipf_ds.groupby("donor"):
        sub = sub.sort_values("severity")
        if sub["severity"].nunique() >= 3:
            slope = stats.linregress(sub["severity"], sub[metric]).slope
            slopes.append((donor, slope))
    vals = np.asarray([v for _, v in slopes], dtype=float)
    trend_rows.append({
        "metric": metric,
        "n_donors": len(vals),
        "mean_within_donor_slope": float(np.mean(vals)) if len(vals) else np.nan,
        "median_within_donor_slope": float(np.median(vals)) if len(vals) else np.nan,
        "n_positive": int(np.sum(vals > 0)),
        "n_negative": int(np.sum(vals < 0)),
        "exact_sign_flip_p_two_sided": exact_signed_mean_p(vals),
        "donor_slopes": ";".join(f"{d}:{v:.6g}" for d, v in slopes),
    })
trend_tests = pd.DataFrame(trend_rows)
trend_tests["FDR"] = bh_adjust(trend_tests["exact_sign_flip_p_two_sided"])
trend_tests.to_csv(OUT_DIR / "severity_trend_exact_tests.csv", index=False)

# Exact B3 IPF versus control comparison at donor level.
control = donor_severity[donor_severity["disease"] == "Control"].copy()
b3 = donor_severity[(donor_severity["disease"] == "IPF") & (donor_severity["severity"] == 3)].copy()
contrast_rows = []
for metric in metric_cols:
    x = b3[metric].dropna().to_numpy(float)
    y = control[metric].dropna().to_numpy(float)
    if len(x) and len(y):
        mw = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact")
        delta = np.mean(x[:, None] > y[None, :]) - np.mean(x[:, None] < y[None, :])
        contrast_rows.append({
            "metric": metric,
            "n_IPF_B3_donors": len(x),
            "n_control_donors": len(y),
            "median_IPF_B3": float(np.median(x)),
            "median_control": float(np.median(y)),
            "IPF_B3_minus_control": float(np.median(x) - np.median(y)),
            "cliffs_delta": float(delta),
            "mann_whitney_exact_p": float(mw.pvalue),
        })
contrast = pd.DataFrame(contrast_rows)
contrast["FDR"] = bh_adjust(contrast["mann_whitney_exact_p"])
contrast.to_csv(OUT_DIR / "IPF_B3_vs_control_donor_tests.csv", index=False)

# Donor-level neighborhood evidence in B3, contrasted with controls.
nb_b3 = neighbor_donor_severity[
    (neighbor_donor_severity["disease"] == "IPF")
    & (neighbor_donor_severity["severity"] == 3)
]
nb_ctrl = neighbor_donor_severity[neighbor_donor_severity["disease"] == "Control"]
nb_tests = []
for exposure in sorted(neighbor_donor_severity["exposure"].unique()):
    for response in sorted(neighbor_donor_severity["response"].unique()):
        x = nb_b3[(nb_b3["exposure"] == exposure) & (nb_b3["response"] == response)]["spearman_rho"].dropna().to_numpy()
        y = nb_ctrl[(nb_ctrl["exposure"] == exposure) & (nb_ctrl["response"] == response)]["spearman_rho"].dropna().to_numpy()
        if len(x) and len(y):
            mw = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact")
            delta = np.mean(x[:, None] > y[None, :]) - np.mean(x[:, None] < y[None, :])
            nb_tests.append({
                "exposure": exposure,
                "response": response,
                "median_rho_IPF_B3": float(np.median(x)),
                "median_rho_control": float(np.median(y)),
                "n_IPF_B3": len(x),
                "n_control": len(y),
                "cliffs_delta": float(delta),
                "mann_whitney_exact_p": float(mw.pvalue),
                "IPF_B3_positive_donors": int(np.sum(x > 0)),
            })
nb_tests = pd.DataFrame(nb_tests)
nb_tests["FDR"] = bh_adjust(nb_tests["mann_whitney_exact_p"])
nb_tests.to_csv(OUT_DIR / "neighbor_B3_vs_control_donor_tests.csv", index=False)

# Figure 1: donor-paired severity trajectories.
sns.set_theme(style="whitegrid", context="talk")
fig_metrics = [
    ("epithelial_secretory_score", "Epithelial secretory remodeling"),
    ("epithelial_alveolar_ACSL4_score", "Epithelial alveolar/ACSL4"),
    ("myeloid_ferroaging_score", "Myeloid ferro-aging module"),
    ("myeloid_TGFB1", "Myeloid-weighted TGFB1"),
]
fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)
for ax, (metric, title) in zip(axes.flat, fig_metrics):
    for donor, sub in ipf_ds.groupby("donor"):
        sub = sub.sort_values("severity")
        ax.plot(sub["severity"], sub[metric], marker="o", linewidth=2, label=donor)
    ax.set_xticks([1, 2, 3], ["B1 mild", "B2 moderate", "B3 severe"])
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Weighted log-normalized score")
axes.flat[0].legend(loc="upper left", ncol=2, frameon=False, fontsize=9)
fig.suptitle("Within-patient spatial programs across IPF fibrosis severity", fontsize=20, y=1.02)
fig.savefig(OUT_DIR / "Figure_1_spatial_severity_trajectories.png", dpi=350, bbox_inches="tight")
fig.savefig(OUT_DIR / "Figure_1_spatial_severity_trajectories.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 2: donor-level neighborhood associations with the primary transition response.
primary_nb = neighbor_donor_severity[neighbor_donor_severity["response"] == "secretory_transition"].copy()
primary_nb["group"] = np.where(
    primary_nb["disease"] == "Control", "Control",
    primary_nb["severity"].map({1: "IPF B1", 2: "IPF B2", 3: "IPF B3"}),
)
order = ["Control", "IPF B1", "IPF B2", "IPF B3"]
fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True, constrained_layout=True)
for ax, exposure in zip(axes, ["TGFB1", "S100A8_ferro", "SPP1"]):
    sub = primary_nb[primary_nb["exposure"] == exposure]
    sns.boxplot(data=sub, x="group", y="spearman_rho", order=order, color="#e7eef8", ax=ax, fliersize=0)
    sns.stripplot(data=sub, x="group", y="spearman_rho", order=order, hue="donor", ax=ax,
                  size=7, jitter=0.08, palette="tab10", legend=False)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_title(exposure.replace("_", " "))
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=35)
    ax.set_ylabel("Neighbor exposure vs epithelial transition\nSpearman rho" if exposure == "TGFB1" else "")
fig.suptitle("Donor-level spatial neighborhood evidence", fontsize=20, y=1.04)
fig.savefig(OUT_DIR / "Figure_2_neighbor_associations.png", dpi=350, bbox_inches="tight")
fig.savefig(OUT_DIR / "Figure_2_neighbor_associations.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 3: representative severe-IPF spatial maps.
b3_sections = sections[
    (sections["disease"] == "IPF") & (sections["severity"] == 3)
    & sections["has_image_coordinates"]
]
representative = b3_sections.sort_values("epithelial_secretory_score", ascending=False).iloc[0]["sample_id"]
plot_df = spots[spots["sample_id"] == representative].copy()
with open(DATA_DIR / f"{representative}_spatial_scalefactors_json.json", "r", encoding="utf-8") as fh:
    scale = json.load(fh)["tissue_lowres_scalef"]
image = plt.imread(DATA_DIR / f"{representative}_spatial_tissue_lowres_image.png")
plot_df["x"] = plot_df["pixel_col"] * scale
plot_df["y"] = plot_df["pixel_row"] * scale
map_cols = [
    ("myeloid_fraction", "Cell2location myeloid fraction", "magma"),
    ("TGFB1_myeloid_sender", "TGFB1 × myeloid sender score", "magma"),
    ("transition_score", "LCN2/SCGB-high vs ACSL4/alveolar-low", "coolwarm"),
    ("neighbor_TGFB1", "Six-neighbor TGFB1 exposure", "magma"),
]
fig, axes = plt.subplots(2, 2, figsize=(14, 12), constrained_layout=True)
for ax, (col, title, cmap) in zip(axes.flat, map_cols):
    ax.imshow(image)
    values = plot_df[col].to_numpy(float)
    lo, hi = np.nanquantile(values, [0.02, 0.98])
    sc = ax.scatter(plot_df["x"], plot_df["y"], c=values, s=8, cmap=cmap,
                    vmin=lo, vmax=hi, linewidths=0, alpha=0.82)
    ax.set_title(title)
    ax.set_axis_off()
    fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.01)
fig.suptitle(f"Representative severe-IPF section: {representative}", fontsize=20, y=1.02)
fig.savefig(OUT_DIR / "Figure_3_representative_B3_spatial_maps.png", dpi=350, bbox_inches="tight")
fig.savefig(OUT_DIR / "Figure_3_representative_B3_spatial_maps.pdf", bbox_inches="tight")
plt.close(fig)

# Provenance and concise machine-readable summary.
def md5(path: Path):
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

summary = {
    "dataset": "E-MTAB-14121 / S-BSST1410",
    "source_matrix": "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/121/E-MTAB-14121",
    "source_cell2location": "https://ftp.ebi.ac.uk/biostudies/fire/S-BSST/410/S-BSST1410",
    "n_sections": int(sections.shape[0]),
    "n_spots": int(spots.shape[0]),
    "n_ipf_donors": int(meta.loc[meta["disease"] == "IPF", "subject_alias"].nunique()),
    "n_control_donors": int(meta.loc[meta["disease"] == "Control", "subject_alias"].nunique()),
    "representative_B3_section": representative,
    "cell2location_zip_md5": md5(DATA_DIR / "cell2location_habermann2020.zip"),
    "analysis_seed": 20261027,
    "permutations_per_section_test": N_PERM,
    "interpretation_limit": "Spatial association and deconvolution do not establish ligand-target causality.",
}
with open(OUT_DIR / "run_info.json", "w", encoding="utf-8") as fh:
    json.dump(summary, fh, ensure_ascii=False, indent=2)

print(json.dumps(summary, indent=2))
print("\nTop B3 vs control donor-level contrasts:")
print(contrast.sort_values("mann_whitney_exact_p").head(12).to_string(index=False))
print("\nPrimary neighborhood tests:")
print(nb_tests[nb_tests["response"] == "secretory_transition"].to_string(index=False))
