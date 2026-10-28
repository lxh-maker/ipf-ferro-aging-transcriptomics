from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from scipy.spatial import cKDTree


ROOT = Path(__file__).resolve().parents[2]
INFILE = ROOT / "06_mechanism_axis" / "spatial_ETAB14121" / "spot_scores_and_neighborhoods.csv.gz"
OUTDIR = INFILE.parent
SCALES = [0, 6, 12, 24, 48, 96]


def rz(x):
    x = np.asarray(x, float)
    med = np.nanmedian(x)
    mad = 1.4826 * np.nanmedian(np.abs(x - med))
    if not np.isfinite(mad) or mad < 1e-8:
        mad = np.nanstd(x)
    return np.zeros_like(x) if (not np.isfinite(mad) or mad < 1e-8) else (x - med) / mad


def rho(x, y):
    good = np.isfinite(x) & np.isfinite(y)
    if good.sum() < 20 or np.std(x[good]) < 1e-10 or np.std(y[good]) < 1e-10:
        return np.nan
    return stats.spearmanr(x[good], y[good]).statistic


spots = pd.read_csv(INFILE)
rows = []
for sample, d in spots.groupby("sample_id", sort=False):
    d = d.reset_index(drop=True)
    coords = d[["array_row", "array_col"]].to_numpy(float)
    max_k = min(max(SCALES) + 1, len(d))
    _, nn = cKDTree(coords).query(coords, k=max_k)
    nn = nn[:, 1:]
    receiver = d["epithelial_receiver"].astype(bool).to_numpy()

    exposures = {
        "TGFB1_myeloid": d["TGFB1_myeloid_sender"].to_numpy(float),
        "TGFB1_all_spots": rz(d["TGFB1"].to_numpy(float)),
        "S100A8_ferro_myeloid": d["S100A8_ferro_sender"].to_numpy(float),
        "S100A8_all_spots": rz(d["S100A8"].to_numpy(float)),
        "SPP1_myeloid": d["SPP1_myeloid_sender"].to_numpy(float),
        "myeloid_fraction": rz(d["myeloid_fraction"].to_numpy(float)),
    }
    responses = {
        "secretory_transition": d["transition_score"].to_numpy(float),
        "LCN2_minus_ACSL4": d["LCN2_ACSL4_balance"].to_numpy(float),
        "abnormal_epithelial_fraction": rz(d["abnormal_epithelial_fraction"].to_numpy(float)),
    }

    for k in SCALES:
        for ename, sender in exposures.items():
            exposure = sender if k == 0 else sender[nn[:, :k]].mean(axis=1)
            for rname, response in responses.items():
                rows.append({
                    "sample_id": sample,
                    "donor": d.loc[0, "donor"],
                    "disease": d.loc[0, "disease"],
                    "severity": int(d.loc[0, "severity"]),
                    "scale_k_neighbors": k,
                    "exposure": ename,
                    "response": rname,
                    "n_receiver_spots": int(receiver.sum()),
                    "spearman_rho": rho(exposure[receiver], response[receiver]),
                })

section = pd.DataFrame(rows)
section.to_csv(OUTDIR / "spatial_scale_sensitivity_section.csv", index=False)

donor = (
    section.groupby(["donor", "disease", "severity", "scale_k_neighbors", "exposure", "response"], as_index=False)
    ["spearman_rho"].mean()
)
donor.to_csv(OUTDIR / "spatial_scale_sensitivity_donor.csv", index=False)

tests = []
for (k, exposure, response), sub in donor.groupby(["scale_k_neighbors", "exposure", "response"]):
    x = sub[(sub["disease"] == "IPF") & (sub["severity"] == 3)]["spearman_rho"].dropna().to_numpy()
    y = sub[sub["disease"] == "Control"]["spearman_rho"].dropna().to_numpy()
    if len(x) and len(y):
        mw = stats.mannwhitneyu(x, y, alternative="two-sided", method="exact")
        delta = np.mean(x[:, None] > y[None, :]) - np.mean(x[:, None] < y[None, :])
        tests.append({
            "scale_k_neighbors": k,
            "exposure": exposure,
            "response": response,
            "median_rho_IPF_B3": np.median(x),
            "median_rho_control": np.median(y),
            "IPF_B3_positive_donors": int(np.sum(x > 0)),
            "cliffs_delta": delta,
            "mann_whitney_exact_p": mw.pvalue,
        })
tests = pd.DataFrame(tests)
tests.to_csv(OUTDIR / "spatial_scale_sensitivity_B3_vs_control.csv", index=False)

sns.set_theme(style="whitegrid", context="talk")
plot_exposures = ["TGFB1_myeloid", "TGFB1_all_spots", "S100A8_ferro_myeloid", "SPP1_myeloid"]
plot = donor[
    donor["exposure"].isin(plot_exposures)
    & (donor["response"] == "secretory_transition")
    & (((donor["disease"] == "IPF") & (donor["severity"] == 3)) | (donor["disease"] == "Control"))
].copy()
plot["group"] = np.where(plot["disease"] == "Control", "Control", "IPF B3")
summary = plot.groupby(["exposure", "group", "scale_k_neighbors"], as_index=False)["spearman_rho"].agg(
    median="median", q25=lambda x: np.quantile(x, 0.25), q75=lambda x: np.quantile(x, 0.75)
)

fig, axes = plt.subplots(2, 2, figsize=(14, 10), sharex=True, sharey=True, constrained_layout=True)
palette = {"Control": "#4c78a8", "IPF B3": "#e45756"}
for ax, exposure in zip(axes.flat, plot_exposures):
    sub = summary[summary["exposure"] == exposure]
    for group in ["Control", "IPF B3"]:
        g = sub[sub["group"] == group].sort_values("scale_k_neighbors")
        x = g["scale_k_neighbors"].to_numpy(float)
        ax.plot(x, g["median"], marker="o", linewidth=2.2, color=palette[group], label=group)
        ax.fill_between(x, g["q25"].to_numpy(float), g["q75"].to_numpy(float),
                        color=palette[group], alpha=0.15)
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_title(exposure.replace("_", " "))
    ax.set_xlabel("Neighborhood size (nearest spots; 0 = same spot)")
    ax.set_ylabel("Median donor Spearman rho")
axes.flat[0].legend(loc="upper right", ncol=1, frameon=False, fontsize=9)
fig.suptitle("Spatial-scale sensitivity of ligand/source associations", fontsize=20, y=1.02)
fig.savefig(OUTDIR / "Figure_4_spatial_scale_sensitivity.png", dpi=350, bbox_inches="tight")
fig.savefig(OUTDIR / "Figure_4_spatial_scale_sensitivity.pdf", bbox_inches="tight")
plt.close(fig)

# Figure 5: focused SPP1 evidence for the cell2location-defined abnormal epithelium.
spp1 = donor[
    (donor["exposure"] == "SPP1_myeloid")
    & (donor["response"] == "abnormal_epithelial_fraction")
].copy()
b3_spp1 = spp1[(spp1["disease"] == "IPF") & (spp1["severity"] == 3)]
fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), constrained_layout=True)
for donor_id, sub in b3_spp1.groupby("donor"):
    sub = sub.sort_values("scale_k_neighbors")
    axes[0].plot(sub["scale_k_neighbors"], sub["spearman_rho"], marker="o", alpha=0.72, label=donor_id)
med = b3_spp1.groupby("scale_k_neighbors", as_index=False)["spearman_rho"].median()
axes[0].plot(med["scale_k_neighbors"], med["spearman_rho"], color="black", marker="o",
             linewidth=3, label="Median")
axes[0].axhline(0, color="black", linewidth=0.8)
axes[0].set_title("SPP1-rich myeloid neighborhood vs abnormal epithelium")
axes[0].set_xlabel("Neighborhood size (nearest spots)")
axes[0].set_ylabel("Spearman rho")
axes[0].legend(frameon=False, ncol=2, fontsize=9)

k96 = spp1[spp1["scale_k_neighbors"] == 96].copy()
k96["group"] = np.where(k96["disease"] == "Control", "Control",
                         np.where(k96["severity"] == 3, "IPF B3", "Other IPF"))
k96 = k96[k96["group"].isin(["Control", "IPF B3"])]
sns.boxplot(data=k96, x="group", y="spearman_rho", order=["Control", "IPF B3"],
            color="#e7eef8", fliersize=0, ax=axes[1])
sns.stripplot(data=k96, x="group", y="spearman_rho", order=["Control", "IPF B3"],
              hue="donor", palette="tab10", size=8, jitter=0.08, legend=False, ax=axes[1])
axes[1].axhline(0, color="black", linewidth=0.8)
axes[1].set_title("Broad neighborhood (96 nearest spots)")
axes[1].set_xlabel("")
axes[1].set_ylabel("SPP1 neighborhood vs abnormal epithelium rho")
axes[1].text(0.5, 0.98, "Exact P=0.057; Cliff's delta=0.875", transform=axes[1].transAxes,
             ha="center", va="top", fontsize=11)
fig.suptitle("Spatial support for an SPP1-rich myeloid–abnormal epithelial niche", fontsize=18, y=1.10)
fig.savefig(OUTDIR / "Figure_5_SPP1_abnormal_epithelial_niche.png", dpi=350, bbox_inches="tight")
fig.savefig(OUTDIR / "Figure_5_SPP1_abnormal_epithelial_niche.pdf", bbox_inches="tight")
plt.close(fig)

print(tests[
    (tests["response"] == "secretory_transition")
    & tests["exposure"].isin(plot_exposures)
].to_string(index=False))
