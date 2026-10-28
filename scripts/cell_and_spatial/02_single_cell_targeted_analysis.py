from __future__ import annotations

import json
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu, spearmanr
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "00_data"
GEO = DATA / "GSE135893"
UCSC = DATA / "GSE135893_UCSC"
OUT = ROOT / "02_results" / "GSE135893_targeted_single_cell"
FIG = OUT / "figures"
TAB = OUT / "tables"
for directory in (OUT, FIG, TAB):
    directory.mkdir(parents=True, exist_ok=True)

sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 150})

CANDIDATES = ["CDO1", "EGR1", "HIF1A", "LCN2", "LIFR", "PDE4B", "PTGS2", "S100A8"]
SIGNATURE = ["PDE4B", "LIFR", "S100A8", "EGR1", "LCN2"]
MECHANISM = ["S100A8", "S100A9", "LCN2", "ACSL4", "LIF", "LIFR", "SPHK1", "PLA2G7"]
DOT_GENES = CANDIDATES + ["S100A9", "ACSL4", "LIF", "SPHK1", "PLA2G7"]

MYELOID = ["Macrophages", "Monocytes", "Proliferating Macrophages", "cDCs", "pDCs"]
FIBRO = ["Fibroblasts", "PLIN2+ Fibroblasts", "HAS1 High Fibroblasts", "Myofibroblasts"]
KEY_CELLTYPES = [
    "Macrophages",
    "Monocytes",
    "AT2",
    "Transitional AT2",
    "MUC5B+",
    "SCGB3A2+",
    "KRT5-/KRT17+",
    "Fibroblasts",
    "PLIN2+ Fibroblasts",
    "HAS1 High Fibroblasts",
    "Myofibroblasts",
    "Endothelial Cells",
]


def read_gene_file(path: Path) -> list[str]:
    return [x.strip() for x in path.read_text(encoding="utf-8-sig").splitlines() if x.strip()]


def fdr(values: pd.Series) -> np.ndarray:
    out = np.full(len(values), np.nan)
    valid = values.notna().to_numpy()
    if valid.any():
        out[valid] = multipletests(values[valid], method="fdr_bh")[1]
    return out


def safe_mwu(a: np.ndarray, b: np.ndarray) -> float:
    a = a[np.isfinite(a)]
    b = b[np.isfinite(b)]
    if len(a) < 3 or len(b) < 3:
        return np.nan
    if np.all(a == a[0]) and np.all(b == b[0]) and a[0] == b[0]:
        return 1.0
    return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)


def load_data():
    dat = np.load(UCSC / "selected_expression.npz", allow_pickle=False)
    genes = dat["genes"].astype(str).tolist()
    expr = dat["expression"].astype(np.uint32, copy=False)
    cells = dat["cells"].astype(str)

    meta = pd.read_csv(GEO / "GSE135893_IPF_metadata.csv.gz", compression="gzip")
    meta = meta.rename(columns={meta.columns[0]: "Cell"})
    if not np.array_equal(cells, meta["Cell"].astype(str).to_numpy()):
        raise ValueError("Expression and metadata cell order do not match")

    coords = pd.read_csv(
        UCSC / "Seurat_UMAP.coords.tsv.gz",
        sep="\t",
        header=None,
        names=["Cell", "UMAP_1", "UMAP_2"],
    )
    if not np.array_equal(cells, coords["Cell"].astype(str).to_numpy()):
        raise ValueError("UMAP and metadata cell order do not match")

    keep = meta["Diagnosis"].isin(["IPF", "Control"]).to_numpy()
    meta = meta.loc[keep].reset_index(drop=True)
    coords = coords.loc[keep].reset_index(drop=True)
    expr = expr[:, keep]
    return genes, expr, meta, coords


def normalize_selected(expr: np.ndarray, totals: np.ndarray) -> np.ndarray:
    totals = np.maximum(totals.astype(np.float64), 1.0)
    return np.log1p(expr.astype(np.float64) * (10_000.0 / totals[None, :])).astype(np.float32)


def build_targeted_score(lognorm: np.ndarray, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    ferro = read_gene_file(DATA / "ferro_aging_genes.txt")
    aliases = {"CAVIN1": "PTRF"}
    mapped = [aliases.get(g, g) for g in ferro]
    used = sorted(set(g for g in mapped if g in genes))
    idx = [genes.index(g) for g in used]
    block = lognorm[idx].astype(np.float64)
    means = block.mean(axis=1, keepdims=True)
    sds = block.std(axis=1, ddof=1, keepdims=True)
    sds[sds == 0] = 1.0
    z = np.clip((block - means) / sds, -5, 5)
    return z.mean(axis=0).astype(np.float32), used


def aggregate_pseudobulk(
    expr: np.ndarray,
    lognorm: np.ndarray,
    score: np.ndarray,
    genes: list[str],
    meta: pd.DataFrame,
) -> pd.DataFrame:
    records = []
    grouped = meta.groupby(["Sample_Name", "Diagnosis", "celltype"], sort=False).indices
    for (sample, diagnosis, celltype), idx in grouped.items():
        idx = np.asarray(idx, dtype=int)
        total_umis = float(meta.loc[idx, "nCount_RNA"].sum())
        sums = expr[:, idx].sum(axis=1, dtype=np.uint64)
        means = lognorm[:, idx].mean(axis=1)
        pct = (expr[:, idx] > 0).mean(axis=1) * 100.0
        log2cpm = np.log2(1.0 + sums.astype(np.float64) * 1_000_000.0 / max(total_umis, 1.0))
        for j, gene in enumerate(genes):
            records.append(
                {
                    "Sample_Name": sample,
                    "Diagnosis": diagnosis,
                    "celltype": celltype,
                    "n_cells": len(idx),
                    "gene": gene,
                    "sum_counts": int(sums[j]),
                    "log2CPM": float(log2cpm[j]),
                    "mean_log_normalized": float(means[j]),
                    "pct_expressing": float(pct[j]),
                    "ferro_aging_targeted_score": float(score[idx].mean()),
                }
            )
    return pd.DataFrame.from_records(records)


def test_pseudobulk(pb: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (celltype, gene), frame in pb.groupby(["celltype", "gene"], sort=False):
        ipf = frame.loc[frame.Diagnosis.eq("IPF"), "log2CPM"].to_numpy()
        ctrl = frame.loc[frame.Diagnosis.eq("Control"), "log2CPM"].to_numpy()
        rows.append(
            {
                "celltype": celltype,
                "gene": gene,
                "n_IPF_samples": len(ipf),
                "n_Control_samples": len(ctrl),
                "median_IPF_log2CPM": np.median(ipf) if len(ipf) else np.nan,
                "median_Control_log2CPM": np.median(ctrl) if len(ctrl) else np.nan,
                "IPF_minus_Control": np.median(ipf) - np.median(ctrl) if len(ipf) and len(ctrl) else np.nan,
                "p_value": safe_mwu(ipf, ctrl),
            }
        )
    out = pd.DataFrame(rows)
    out["FDR"] = fdr(out["p_value"])
    return out.sort_values(["FDR", "p_value"], na_position="last")


def score_tests(pb: pd.DataFrame) -> pd.DataFrame:
    scores = pb.drop_duplicates(["Sample_Name", "Diagnosis", "celltype"])
    rows = []
    for celltype, frame in scores.groupby("celltype"):
        ipf = frame.loc[frame.Diagnosis.eq("IPF"), "ferro_aging_targeted_score"].to_numpy()
        ctrl = frame.loc[frame.Diagnosis.eq("Control"), "ferro_aging_targeted_score"].to_numpy()
        rows.append(
            {
                "celltype": celltype,
                "n_IPF_samples": len(ipf),
                "n_Control_samples": len(ctrl),
                "median_IPF": np.median(ipf) if len(ipf) else np.nan,
                "median_Control": np.median(ctrl) if len(ctrl) else np.nan,
                "IPF_minus_Control": np.median(ipf) - np.median(ctrl) if len(ipf) and len(ctrl) else np.nan,
                "p_value": safe_mwu(ipf, ctrl),
            }
        )
    out = pd.DataFrame(rows)
    out["FDR"] = fdr(out["p_value"])
    return out.sort_values(["FDR", "p_value"], na_position="last")


def proportion_tests(meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    counts = meta.groupby(["Sample_Name", "Diagnosis", "celltype"]).size().rename("n_cells").reset_index()
    totals = meta.groupby(["Sample_Name", "Diagnosis"]).size().rename("sample_total").reset_index()
    prop = counts.merge(totals, on=["Sample_Name", "Diagnosis"], how="left")
    prop["proportion"] = prop["n_cells"] / prop["sample_total"]
    all_samples = meta[["Sample_Name", "Diagnosis"]].drop_duplicates()
    all_types = pd.DataFrame({"celltype": sorted(meta.celltype.unique())})
    grid = all_samples.merge(all_types, how="cross")
    prop = grid.merge(prop, on=["Sample_Name", "Diagnosis", "celltype"], how="left").fillna(
        {"n_cells": 0, "proportion": 0}
    )
    rows = []
    for celltype, frame in prop.groupby("celltype"):
        ipf = frame.loc[frame.Diagnosis.eq("IPF"), "proportion"].to_numpy()
        ctrl = frame.loc[frame.Diagnosis.eq("Control"), "proportion"].to_numpy()
        rows.append(
            {
                "celltype": celltype,
                "median_IPF": np.median(ipf),
                "median_Control": np.median(ctrl),
                "IPF_minus_Control": np.median(ipf) - np.median(ctrl),
                "p_value": safe_mwu(ipf, ctrl),
            }
        )
    tests = pd.DataFrame(rows)
    tests["FDR"] = fdr(tests["p_value"])
    return prop, tests.sort_values(["FDR", "p_value"])


def plot_umap(meta: pd.DataFrame, coords: pd.DataFrame) -> None:
    rng = np.random.default_rng(20261027)
    n = min(len(meta), 70_000)
    idx = np.sort(rng.choice(len(meta), n, replace=False)) if len(meta) > n else np.arange(len(meta))
    frame = pd.concat([coords.iloc[idx].reset_index(drop=True), meta.iloc[idx].reset_index(drop=True)], axis=1)

    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    sns.scatterplot(
        data=frame,
        x="UMAP_1",
        y="UMAP_2",
        hue="celltype",
        palette="husl",
        s=2,
        linewidth=0,
        alpha=0.65,
        ax=axes[0],
    )
    axes[0].set_title("GSE135893: author-annotated cell types")
    axes[0].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, markerscale=3)
    sns.scatterplot(
        data=frame,
        x="UMAP_1",
        y="UMAP_2",
        hue="Diagnosis",
        hue_order=["Control", "IPF"],
        palette={"Control": "#3B82F6", "IPF": "#DC2626"},
        s=2,
        linewidth=0,
        alpha=0.55,
        ax=axes[1],
    )
    axes[1].set_title("Strict comparison: Control vs IPF")
    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel("")
        ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(FIG / "01_UMAP_celltypes_and_groups.pdf", bbox_inches="tight")
    fig.savefig(FIG / "01_UMAP_celltypes_and_groups.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_dotplot(expr: np.ndarray, lognorm: np.ndarray, genes: list[str], meta: pd.DataFrame) -> pd.DataFrame:
    gene_use = [g for g in DOT_GENES if g in genes]
    celltypes = meta.celltype.value_counts().index.tolist()
    rows = []
    for ct in celltypes:
        idx = np.flatnonzero(meta.celltype.eq(ct).to_numpy())
        for gene in gene_use:
            j = genes.index(gene)
            rows.append(
                {
                    "celltype": ct,
                    "gene": gene,
                    "mean": float(lognorm[j, idx].mean()),
                    "pct": float((expr[j, idx] > 0).mean() * 100),
                }
            )
    dot = pd.DataFrame(rows)
    dot["scaled_mean"] = dot.groupby("gene")["mean"].transform(
        lambda x: (x - x.mean()) / (x.std(ddof=1) if x.std(ddof=1) > 0 else 1)
    )
    xmap = {g: i for i, g in enumerate(gene_use)}
    ymap = {ct: i for i, ct in enumerate(celltypes[::-1])}
    fig, ax = plt.subplots(figsize=(10, max(8, len(celltypes) * 0.32)))
    sc = ax.scatter(
        dot.gene.map(xmap),
        dot.celltype.map(ymap),
        s=np.maximum(dot.pct.to_numpy(), 0.5) * 3.2,
        c=dot.scaled_mean,
        cmap="RdBu_r",
        vmin=-2,
        vmax=2,
        edgecolors="none",
    )
    ax.set_xticks(range(len(gene_use)), gene_use, rotation=45, ha="right")
    ax.set_yticks(range(len(celltypes)), celltypes[::-1])
    ax.set_title("Candidate and mechanism genes across lung cell types")
    ax.set_xlabel("")
    ax.set_ylabel("")
    plt.colorbar(sc, ax=ax, label="Scaled mean expression")
    fig.tight_layout()
    fig.savefig(FIG / "02_candidate_mechanism_dotplot.pdf", bbox_inches="tight")
    fig.savefig(FIG / "02_candidate_mechanism_dotplot.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    return dot


def plot_fc_heatmap(tests: pd.DataFrame) -> None:
    frame = tests[tests.gene.isin(DOT_GENES) & tests.celltype.isin(KEY_CELLTYPES)].copy()
    mat = frame.pivot(index="celltype", columns="gene", values="IPF_minus_Control")
    mat = mat.reindex(index=[x for x in KEY_CELLTYPES if x in mat.index], columns=[x for x in DOT_GENES if x in mat.columns])
    fig, ax = plt.subplots(figsize=(10, 7))
    sns.heatmap(mat, cmap="RdBu_r", center=0, annot=True, fmt=".2f", linewidths=0.3, ax=ax)
    ax.set_title("Patient-level pseudobulk difference (IPF − Control)")
    ax.set_xlabel("")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(FIG / "03_pseudobulk_effect_heatmap.pdf", bbox_inches="tight")
    fig.savefig(FIG / "03_pseudobulk_effect_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def compartment_gene_table(pb: pd.DataFrame) -> pd.DataFrame:
    contexts = {
        "S100A8": MYELOID,
        "PDE4B": MYELOID,
        "LCN2": ["MUC5B+", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "AT2", "Transitional AT2", "KRT5-/KRT17+"],
        "EGR1": ["MUC5B+", "SCGB3A2+", "SCGB3A2+ SCGB1A1+", "AT2", "Transitional AT2", "KRT5-/KRT17+"],
        "LIFR": FIBRO,
    }
    rows = []
    for gene, celltypes in contexts.items():
        frame = pb[(pb.gene == gene) & pb.celltype.isin(celltypes)].copy()
        if frame.empty:
            continue
        weighted = []
        for (sample, diagnosis), sf in frame.groupby(["Sample_Name", "Diagnosis"]):
            value = np.average(sf.log2CPM, weights=np.maximum(sf.n_cells, 1))
            weighted.append({"gene": gene, "Sample_Name": sample, "Diagnosis": diagnosis, "log2CPM": value})
        rows.extend(weighted)
    return pd.DataFrame(rows)


def plot_signature_context(pb: pd.DataFrame) -> pd.DataFrame:
    frame = compartment_gene_table(pb)
    fig, axes = plt.subplots(1, 5, figsize=(17, 4.2), sharey=False)
    stats = []
    for ax, gene in zip(axes, SIGNATURE):
        sub = frame[frame.gene.eq(gene)]
        ipf = sub.loc[sub.Diagnosis.eq("IPF"), "log2CPM"].to_numpy()
        ctrl = sub.loc[sub.Diagnosis.eq("Control"), "log2CPM"].to_numpy()
        p = safe_mwu(ipf, ctrl)
        stats.append({"gene": gene, "p_value": p, "median_IPF": np.median(ipf), "median_Control": np.median(ctrl)})
        sns.boxplot(data=sub, x="Diagnosis", y="log2CPM", order=["Control", "IPF"], color="white", ax=ax)
        sns.stripplot(
            data=sub,
            x="Diagnosis",
            y="log2CPM",
            order=["Control", "IPF"],
            palette={"Control": "#3B82F6", "IPF": "#DC2626"},
            size=5,
            jitter=0.12,
            ax=ax,
        )
        ax.set_title(f"{gene}\np={p:.3g}" if np.isfinite(p) else gene)
        ax.set_xlabel("")
        ax.set_ylabel("log2 CPM" if ax is axes[0] else "")
        if ax.get_legend() is not None:
            ax.get_legend().remove()
    fig.suptitle("Five-gene signature in biologically relevant compartments", y=1.04)
    fig.tight_layout()
    fig.savefig(FIG / "04_five_gene_patient_level_boxplots.pdf", bbox_inches="tight")
    fig.savefig(FIG / "04_five_gene_patient_level_boxplots.png", bbox_inches="tight", dpi=300)
    plt.close(fig)
    out = pd.DataFrame(stats)
    out["FDR"] = fdr(out.p_value)
    return out


def plot_score(score_pb: pd.DataFrame) -> None:
    top = score_pb[
        score_pb.celltype.isin(KEY_CELLTYPES)
        & (score_pb.n_IPF_samples >= 5)
        & (score_pb.n_Control_samples >= 5)
    ].sort_values("IPF_minus_Control")
    if top.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 6))
    colors = np.where(top.IPF_minus_Control >= 0, "#DC2626", "#3B82F6")
    ax.barh(top.celltype, top.IPF_minus_Control, color=colors)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Median targeted ferro-aging score difference (IPF − Control)")
    ax.set_ylabel("")
    ax.set_title("Patient-level ferro-aging score shifts by cell type")
    fig.tight_layout()
    fig.savefig(FIG / "05_ferro_aging_score_effects.pdf", bbox_inches="tight")
    fig.savefig(FIG / "05_ferro_aging_score_effects.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def plot_proportions(prop: pd.DataFrame, prop_tests: pd.DataFrame) -> None:
    use = [x for x in KEY_CELLTYPES if x in prop.celltype.unique()]
    frame = prop[prop.celltype.isin(use)].copy()
    fig, ax = plt.subplots(figsize=(12, 6))
    sns.boxplot(
        data=frame,
        x="celltype",
        y="proportion",
        hue="Diagnosis",
        hue_order=["Control", "IPF"],
        palette={"Control": "#3B82F6", "IPF": "#DC2626"},
        showfliers=False,
        ax=ax,
    )
    ax.tick_params(axis="x", rotation=55)
    ax.set_xlabel("")
    ax.set_ylabel("Fraction of cells per donor")
    ax.set_title("Cell-type composition by donor")
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(FIG / "06_celltype_proportions.pdf", bbox_inches="tight")
    fig.savefig(FIG / "06_celltype_proportions.png", bbox_inches="tight", dpi=300)
    plt.close(fig)


def mechanism_correlations(pb: pd.DataFrame) -> pd.DataFrame:
    contexts = {
        "Myeloid": MYELOID,
        "Fibroblast": FIBRO,
        "Epithelial": ["AT1", "AT2", "Transitional AT2", "MUC5B+", "SCGB3A2+", "KRT5-/KRT17+", "Basal"],
    }
    pairs = {
        "Myeloid": [("S100A8", "LCN2"), ("S100A8", "ACSL4"), ("S100A9", "LCN2")],
        "Epithelial": [("LCN2", "ACSL4"), ("LCN2", "PTGS2")],
        "Fibroblast": [("LIFR", "COL1A1"), ("LIFR", "COL3A1"), ("LIFR", "CTHRC1")],
    }
    rows = []
    for context, celltypes in contexts.items():
        frame = pb[pb.celltype.isin(celltypes)]
        agg = (
            frame.groupby(["Sample_Name", "Diagnosis", "gene"], as_index=False)
            .apply(lambda x: pd.Series({"log2CPM": np.average(x.log2CPM, weights=np.maximum(x.n_cells, 1))}), include_groups=False)
        )
        wide = agg.pivot_table(index=["Sample_Name", "Diagnosis"], columns="gene", values="log2CPM").reset_index()
        for diagnosis in ["IPF", "Control", "All"]:
            use = wide if diagnosis == "All" else wide[wide.Diagnosis.eq(diagnosis)]
            for gene_a, gene_b in pairs[context]:
                if gene_a not in use or gene_b not in use:
                    continue
                valid = use[[gene_a, gene_b]].dropna()
                rho, p = spearmanr(valid[gene_a], valid[gene_b]) if len(valid) >= 4 else (np.nan, np.nan)
                rows.append(
                    {
                        "context": context,
                        "diagnosis": diagnosis,
                        "gene_a": gene_a,
                        "gene_b": gene_b,
                        "n_samples": len(valid),
                        "spearman_rho": rho,
                        "p_value": p,
                    }
                )
    out = pd.DataFrame(rows)
    out["FDR"] = fdr(out.p_value)
    return out.sort_values(["FDR", "p_value"], na_position="last")


def write_summary(meta: pd.DataFrame, gene_tests: pd.DataFrame, score_test: pd.DataFrame, sig_stats: pd.DataFrame, used: list[str]) -> None:
    sample_counts = meta.groupby("Diagnosis").Sample_Name.nunique().to_dict()
    key = gene_tests[gene_tests.gene.isin(DOT_GENES) & gene_tests.celltype.isin(KEY_CELLTYPES)].head(20)
    score_key = score_test[(score_test.n_IPF_samples >= 5) & (score_test.n_Control_samples >= 5)].head(10)
    lines = [
        "# GSE135893 targeted single-cell analysis summary",
        "",
        f"- Strict cohort: {sample_counts.get('IPF', 0)} IPF donors, {sample_counts.get('Control', 0)} controls.",
        f"- Cells retained: {len(meta):,} (IPF and Control only).",
        f"- Author-provided cell types: {meta.celltype.nunique()}.",
        f"- Ferro-aging genes used in targeted composite score: {len(used)}.",
        "- Differential tests are donor-level pseudobulk tests; cells are not treated as independent biological replicates.",
        "- The targeted ferro-aging score is a standardized composite over the supplied gene set, not a full-transcriptome AUCell/UCell score.",
        "",
        "## Five-gene compartment results",
        "",
        "```text",
        sig_stats.to_string(index=False),
        "```",
        "",
        "## Top candidate/mechanism pseudobulk results",
        "",
        "```text",
        key.to_string(index=False),
        "```",
        "",
        "## Top ferro-aging score shifts",
        "",
        "```text",
        score_key.to_string(index=False),
        "```",
    ]
    (OUT / "ANALYSIS_SUMMARY.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    genes, expr, meta, coords = load_data()
    lognorm = normalize_selected(expr, meta.nCount_RNA.to_numpy())
    score, ferro_used = build_targeted_score(lognorm, genes)
    meta["ferro_aging_targeted_score"] = score

    pb = aggregate_pseudobulk(expr, lognorm, score, genes, meta)
    gene_tests = test_pseudobulk(pb)
    score_test = score_tests(pb)
    prop, prop_tests = proportion_tests(meta)
    correlations = mechanism_correlations(pb)

    meta.to_csv(TAB / "strict_IPF_Control_cell_metadata.csv.gz", index=False, compression="gzip")
    pb.to_csv(TAB / "patient_celltype_gene_pseudobulk.csv.gz", index=False, compression="gzip")
    gene_tests.to_csv(TAB / "patient_level_gene_tests.csv", index=False, encoding="utf-8-sig")
    score_test.to_csv(TAB / "ferro_aging_score_tests.csv", index=False, encoding="utf-8-sig")
    prop.to_csv(TAB / "celltype_proportions_by_patient.csv", index=False, encoding="utf-8-sig")
    prop_tests.to_csv(TAB / "celltype_proportion_tests.csv", index=False, encoding="utf-8-sig")
    correlations.to_csv(TAB / "mechanism_correlations.csv", index=False, encoding="utf-8-sig")
    (TAB / "ferro_aging_genes_used.txt").write_text("\n".join(ferro_used) + "\n", encoding="utf-8")

    np.savez_compressed(
        OUT / "GSE135893_IPF_Control_targeted_expression.npz",
        genes=np.asarray(genes),
        expression=expr,
        log_normalized=lognorm,
        ferro_aging_targeted_score=score,
        cells=meta.Cell.astype(str).to_numpy(),
    )

    plot_umap(meta, coords)
    dot = plot_dotplot(expr, lognorm, genes, meta)
    dot.to_csv(TAB / "candidate_mechanism_dotplot_data.csv", index=False, encoding="utf-8-sig")
    plot_fc_heatmap(gene_tests)
    sig_stats = plot_signature_context(pb)
    sig_stats.to_csv(TAB / "five_gene_compartment_tests.csv", index=False, encoding="utf-8-sig")
    plot_score(score_test)
    plot_proportions(prop, prop_tests)
    write_summary(meta, gene_tests, score_test, sig_stats, ferro_used)

    run_info = {
        "dataset": "GSE135893",
        "comparison": "IPF vs Control only",
        "cells": int(len(meta)),
        "IPF_donors": int(meta.loc[meta.Diagnosis.eq("IPF"), "Sample_Name"].nunique()),
        "Control_donors": int(meta.loc[meta.Diagnosis.eq("Control"), "Sample_Name"].nunique()),
        "genes_downloaded": len(genes),
        "ferro_aging_genes_used": len(ferro_used),
        "notes": [
            "Author cell-type annotations were retained.",
            "Gene comparisons use donor-level pseudobulk summaries.",
            "Targeted score is provisional until a full-transcriptome AUCell/UCell score is run.",
        ],
    }
    (OUT / "run_info.json").write_text(json.dumps(run_info, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(run_info, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
