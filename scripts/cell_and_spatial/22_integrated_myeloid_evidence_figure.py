#!/usr/bin/env python3
"""Create an integrated, evidence-calibrated summary figure."""

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "06_mechanism_axis" / "myeloid_state_coupling"
DISC = BASE / "discovery_GSE135893"
VAL = BASE / "validation_GSE136831"
SPA = BASE / "spatial_ETAB14121"
OUT = BASE

ORDER = ["Resident/alveolar", "S100A8 transitional", "Lysosome/lipid/iron", "SPP1 effector"]
COLORS = dict(zip(ORDER, ["#4C78A8", "#E45756", "#F2A104", "#7A5195"]))

disc = pd.read_csv(DISC / "module_state_positions_by_donor.csv")
valpos = pd.read_csv(VAL / "fixed_module_positions_by_patient.csv")
pbulk = pd.read_csv(VAL / "fixed_module_scores_by_patient_and_author_celltype.csv")
coloc = pd.read_csv(SPA / "donor_severity_module_colocalization.csv")

sns.set_theme(style="white", context="paper", font_scale=1.05)
fig = plt.figure(figsize=(16, 10))
gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.18], height_ratios=[1, 1], hspace=.38, wspace=.35)

ax = fig.add_subplot(gs[0, 0])
sns.stripplot(data=disc, x="module", y="median_top_quartile_pseudotime", order=ORDER,
              hue="module", palette=COLORS, legend=False, jitter=.14, size=5, alpha=.75, ax=ax)
sns.boxplot(data=disc, x="module", y="median_top_quartile_pseudotime", order=ORDER,
            color="white", width=.5, showfliers=False, ax=ax)
ax.tick_params(axis="x", rotation=28); ax.set_xlabel("")
ax.set_ylabel("Top-quartile cell DPT")
ax.set_title("A  Discovery: partial ordering", loc="left", fontweight="bold")
ax.text(.02, .98, "Resident earlier; transition≈processing\nSPP1 later in 10/12 donors",
        transform=ax.transAxes, va="top", fontsize=8)

ax = fig.add_subplot(gs[0, 1])
types = ["cMonocyte", "ncMonocyte", "Macrophage_Alveolar", "Macrophage"]
heat = pbulk.groupby("Manuscript_Identity")[ORDER].median().reindex(types)
sns.heatmap(heat, cmap="vlag", center=0, annot=True, fmt=".2f", linewidths=.5,
            cbar_kws={"label": "Median patient module score"}, ax=ax)
ax.set_xlabel(""); ax.set_ylabel("")
ax.set_title("B  Validation: branch identity", loc="left", fontweight="bold")

ax = fig.add_subplot(gs[0, 2])
sns.stripplot(data=valpos, x="module", y="median_top_quartile_pseudotime", order=ORDER,
              hue="module", palette=COLORS, legend=False, jitter=.14, size=3.7, alpha=.62, ax=ax)
sns.boxplot(data=valpos, x="module", y="median_top_quartile_pseudotime", order=ORDER,
            color="white", width=.5, showfliers=False, ax=ax)
ax.tick_params(axis="x", rotation=28); ax.set_xlabel("")
ax.set_ylabel("Top-quartile cell DPT")
ax.set_title("C  Independent fixed-module projection", loc="left", fontweight="bold")
ax.text(.02, .98, "Complete four-state order: 0/32 patients",
        transform=ax.transAxes, va="top", fontsize=8)

ax = fig.add_subplot(gs[1, 0:2])
sub = coloc[(coloc.left == "Lysosome/lipid/iron") & (coloc.right == "SPP1 effector")].copy()
sub["group"] = np.where(sub.disease == "Control", "Control", sub.severity.map({1: "IPF B1", 2: "IPF B2", 3: "IPF B3"}))
group_order = ["Control", "IPF B1", "IPF B2", "IPF B3"]
for donor, dat in sub[sub.disease == "IPF"].groupby("donor"):
    dat = dat.sort_values("severity")
    ax.plot(dat.group.map({g: i for i, g in enumerate(group_order)}), dat.partial_spearman,
            color="#7A5195", alpha=.35, lw=1.2)
sns.stripplot(data=sub, x="group", y="partial_spearman", order=group_order,
              color="#7A5195", jitter=.08, size=6, alpha=.8, ax=ax)
sns.boxplot(data=sub, x="group", y="partial_spearman", order=group_order,
            color="white", width=.5, showfliers=False, ax=ax)
ax.axhline(0, color="grey", lw=.8)
ax.set_xlabel(""); ax.set_ylabel("Partial Spearman rho in myeloid-rich spots")
ax.set_title("D  Processing–SPP1 spatial colocalization", loc="left", fontweight="bold")
ax.text(.02, .97, "B3 vs control exact P=0.0286; FDR=0.114; 4/4 B3 donors positive",
        transform=ax.transAxes, va="top", fontsize=9)

ax = fig.add_subplot(gs[1, 2]); ax.axis("off")
ax.set_title("E  Evidence-calibrated model", loc="left", fontweight="bold")
box = dict(boxstyle="round,pad=.45", ec="none")
ax.text(.05, .76, "S100A8-high\ninflammatory monocytes", ha="center", va="center",
        transform=ax.transAxes, bbox={**box, "fc": "#F7C7C6"}, fontsize=10)
ax.text(.05, .25, "Resident/alveolar\nmacrophages", ha="center", va="center",
        transform=ax.transAxes, bbox={**box, "fc": "#C8D9EC"}, fontsize=10)
ax.text(.50, .25, "Lysosome–lipid–iron\nprocessing macrophages", ha="center", va="center",
        transform=ax.transAxes, bbox={**box, "fc": "#F9D995"}, fontsize=10)
ax.text(.88, .25, "SPP1-rich\nmacrophages", ha="center", va="center",
        transform=ax.transAxes, bbox={**box, "fc": "#D8C7E5"}, fontsize=10)
ax.text(.50, .76, "LCN2-high / ACSL4-low\nabnormal epithelium", ha="center", va="center",
        transform=ax.transAxes, bbox={**box, "fc": "#D5D5D5"}, fontsize=10)
ax.annotate("", xy=(.36, .25), xytext=(.18, .25), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="->", lw=2, color="#555"))
ax.annotate("", xy=(.77, .25), xytext=(.65, .25), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="->", lw=2, color="#555"))
ax.annotate("not resolved as\na linear precursor", xy=(.38, .39), xytext=(.12, .65), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", lw=1.5, ls="--", color="#E45756"),
            color="#B23A38", ha="center", fontsize=8)
ax.annotate("spatial link\nnot significant", xy=(.67, .66), xytext=(.78, .42), xycoords="axes fraction",
            arrowprops=dict(arrowstyle="-|>", lw=1.5, ls="--", color="#777"),
            color="#666", ha="center", fontsize=8)
ax.text(.5, .03, "Solid: supported association  |  Dashed: unproven connection",
        transform=ax.transAxes, ha="center", fontsize=8, color="#555")

fig.suptitle("Cross-cohort evaluation of the proposed IPF myeloid state continuum",
             fontsize=17, fontweight="bold", y=.995)
fig.savefig(OUT / "Figure_7_integrated_myeloid_evidence.png", dpi=400, bbox_inches="tight")
fig.savefig(OUT / "Figure_7_integrated_myeloid_evidence.pdf", bbox_inches="tight")
plt.close(fig)
print(OUT / "Figure_7_integrated_myeloid_evidence.png")
