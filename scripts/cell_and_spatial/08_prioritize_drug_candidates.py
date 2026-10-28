from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[2]
DRUG_DIR = ROOT / "05_drug_screen" / "DrugReflector_GSE32537_IPF"
REVERSAL = DRUG_DIR / "GSE32537_IPF_DrugReflector_reversal_top50.csv"
BULK_DE = ROOT / "05_drug_screen" / "GSE32537_signature" / "GSE32537_IPF_vs_Control_limma_all.csv"
SC_TESTS = ROOT / "02_results" / "GSE135893_targeted_single_cell" / "tables" / "patient_level_gene_tests.csv"

drugs = pd.read_csv(REVERSAL)
de = pd.read_csv(BULK_DE).set_index("Gene")
sc = pd.read_csv(SC_TESTS)

# Literature triage is deliberately separated from the model score. A high
# classifier rank is not by itself therapeutic evidence.
evidence = {
    "PF-543": {
        "Evidence_Tier": "A_direct_pulmonary_fibrosis",
        "Recommendation": "Priority 1",
        "Rationale": "SPHK1 inhibitor; post-exposure treatment reduced fibrosis in bleomycin and non-resolving asbestos mouse models, with less epithelial mtDNA damage and fibrogenic monocyte recruitment.",
        "Evidence_URL": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7460639/",
    },
    "darapladib": {
        "Evidence_Tier": "A_direct_pulmonary_fibrosis",
        "Recommendation": "Priority 1",
        "Rationale": "PLA2G7/Lp-PLA2 inhibitor; reduced silica-induced pulmonary fibrosis in mice. PLA2G7 is upregulated in GSE32537 IPF lung and enriched in IPF monocytes/macrophages in GSE135893.",
        "Evidence_URL": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12206922/",
    },
    "NVP-BGJ398": {
        "Evidence_Tier": "B_direct_but_bidirectional",
        "Recommendation": "Priority 2 / caution",
        "Rationale": "Pan-FGFR inhibitor reduced collagen and histologic fibrosis in bleomycin mice, but FGFR blockade had bidirectional effects and increased early mortality in that model.",
        "Evidence_URL": "https://tokushima-u.repo.nii.ac.jp/record/2008014/files/k3430_fulltext.pdf",
    },
    "D-4476": {
        "Evidence_Tier": "C_pathway_rationale",
        "Recommendation": "Exploratory",
        "Rationale": "TGF-beta/CK1-associated inhibitor with pathway relevance, but low DrugReflector rank and no direct IPF validation identified here.",
        "Evidence_URL": "",
    },
    "toremifene": {
        "Evidence_Tier": "D_conflicting_or_risk",
        "Recommendation": "Deprioritize",
        "Rationale": "SERM signal is not IPF-specific; related tamoxifen has clinical evidence of increased radiation-associated pulmonary fibrosis.",
        "Evidence_URL": "https://pubmed.ncbi.nlm.nih.gov/8656444/",
    },
    "tamoxifen": {
        "Evidence_Tier": "D_conflicting_or_risk",
        "Recommendation": "Deprioritize",
        "Rationale": "Clinical studies associate concurrent tamoxifen/radiotherapy with increased pulmonary fibrosis risk.",
        "Evidence_URL": "https://pubmed.ncbi.nlm.nih.gov/8656444/",
    },
    "PHA-665752": {
        "Evidence_Tier": "D_conflicting_or_risk",
        "Recommendation": "Deprioritize",
        "Rationale": "MET inhibitor developed for cancer; no direct IPF benefit found, and HGF-MET signaling can support epithelial repair.",
        "Evidence_URL": "https://pubmed.ncbi.nlm.nih.gov/14612533/",
    },
    "fenticonazole": {
        "Evidence_Tier": "E_no_direct_mechanistic_support",
        "Recommendation": "Deprioritize",
        "Rationale": "Topical antifungal; no direct pulmonary-fibrosis mechanism or validation identified.",
        "Evidence_URL": "",
    },
    "liranaftate": {
        "Evidence_Tier": "E_no_direct_mechanistic_support",
        "Recommendation": "Deprioritize",
        "Rationale": "Topical antifungal; no direct pulmonary-fibrosis mechanism or validation identified.",
        "Evidence_URL": "",
    },
}

rows = []
for _, row in drugs.iterrows():
    name = str(row.get("Compound_Name", ""))
    if not name or name.startswith("BRD-"):
        continue
    record = row.to_dict()
    record.update(
        evidence.get(
            name,
            {
                "Evidence_Tier": "E_not_reviewed",
                "Recommendation": "Not prioritized",
                "Rationale": "Named candidate without a direct pulmonary-fibrosis evidence match in the targeted review.",
                "Evidence_URL": "",
            },
        )
    )
    targets = [x.strip() for x in str(row.get("Target", "")).split("|") if x.strip()]
    target = next((x for x in targets if x in de.index), None)
    record["Matched_Bulk_Target"] = target or ""
    record["Bulk_IPF_logFC"] = de.loc[target, "logFC"] if target else float("nan")
    record["Bulk_IPF_FDR"] = de.loc[target, "adj.P.Val"] if target else float("nan")
    if target and target in set(sc["gene"]):
        cell_hits = sc[(sc["gene"] == target) & sc["p_value"].notna()].sort_values("p_value")
        if len(cell_hits):
            best = cell_hits.iloc[0]
            record["Best_scRNA_celltype"] = best["celltype"]
            record["Best_scRNA_IPF_minus_Control"] = best["IPF_minus_Control"]
            record["Best_scRNA_nominal_P"] = best["p_value"]
            record["Best_scRNA_global_FDR"] = best["FDR"]
    rows.append(record)

priority_order = {
    "Priority 1": 0,
    "Priority 2 / caution": 1,
    "Exploratory": 2,
    "Not prioritized": 3,
    "Deprioritize": 4,
}
prioritized = pd.DataFrame(rows)
prioritized["Priority_Order"] = prioritized["Recommendation"].map(priority_order).fillna(9)
prioritized = prioritized.sort_values(["Priority_Order", "Rank"]).drop(columns="Priority_Order")
prioritized.to_csv(DRUG_DIR / "GSE32537_IPF_named_candidates_evidence_triage.csv", index=False, encoding="utf-8-sig")

plot = prioritized.head(18).copy()
plot["Display"] = plot["Compound_Name"] + "  (rank " + plot["Rank"].astype(int).astype(str) + ")"
order = plot.sort_values("Probability", ascending=True)["Display"]
plot["Display"] = pd.Categorical(plot["Display"], categories=order, ordered=True)
palette = {
    "A_direct_pulmonary_fibrosis": "#D62728",
    "B_direct_but_bidirectional": "#FF7F0E",
    "C_pathway_rationale": "#9467BD",
    "D_conflicting_or_risk": "#7F7F7F",
    "E_no_direct_mechanistic_support": "#BDBDBD",
    "E_not_reviewed": "#9ECAE1",
}
sns.set_theme(style="whitegrid")
fig, ax = plt.subplots(figsize=(9, 7))
sns.barplot(data=plot, x="Probability", y="Display", hue="Evidence_Tier", palette=palette, dodge=False, ax=ax)
ax.set_title("Named DrugReflector reversal candidates: rank plus evidence triage")
ax.set_xlabel("DrugReflector probability (ranking score, not clinical probability)")
ax.set_ylabel("")
ax.legend(title="Evidence tier", loc="lower right", fontsize=8)
fig.tight_layout()
fig.savefig(DRUG_DIR / "GSE32537_IPF_named_candidates_evidence_triage.png", dpi=300)
fig.savefig(DRUG_DIR / "GSE32537_IPF_named_candidates_evidence_triage.pdf")
plt.close(fig)

print(prioritized[["Rank", "Compound_Name", "Target", "Evidence_Tier", "Recommendation", "Bulk_IPF_logFC", "Bulk_IPF_FDR"]].to_string(index=False))

