from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
LEGACY_MODEL_ROOT = Path(
    os.environ.get("DRUGREFLECTOR_ROOT", ROOT / "external" / "drugreflector")
)
SIGNATURE_FILE = ROOT / "05_drug_screen" / "GSE32537_signature" / "GSE32537_IPF_vs_Control_vscores.csv"
OUT = ROOT / "05_drug_screen" / "DrugReflector_GSE32537_IPF"
OUT.mkdir(parents=True, exist_ok=True)

sys.path.insert(0, str(LEGACY_MODEL_ROOT))
import drugreflector as dr  # noqa: E402


signature = pd.read_csv(SIGNATURE_FILE)
disease_vscore = pd.Series(
    signature["vscore"].astype(float).to_numpy(),
    index=signature["Gene"].astype(str).to_numpy(),
    name="IPF_vs_Control",
)
disease_vscore = disease_vscore[~disease_vscore.index.duplicated(keep="first")].dropna()

checkpoint_paths = [
    LEGACY_MODEL_ROOT / "checkpoints" / f"model_fold_{fold}.pt" for fold in range(3)
]
model = dr.DrugReflector(checkpoint_paths=[str(p) for p in checkpoint_paths])
coverage = model.check_gene_coverage(disease_vscore.index)


def predict_one(vscore: pd.Series, label: str) -> pd.DataFrame:
    predictions = model.predict(vscore, n_top=200)
    rank_col = next(col for col in predictions.columns if col[0] == "rank")
    prob_col = next(col for col in predictions.columns if col[0] == "prob")
    logit_col = next(col for col in predictions.columns if col[0] == "logit")
    result = pd.DataFrame(
        {
            "Drug_ID": predictions.index.astype(str),
            "Rank": predictions[rank_col].astype(int).to_numpy(),
            "Probability": predictions[prob_col].astype(float).to_numpy(),
            "Logit": predictions[logit_col].astype(float).to_numpy(),
        }
    ).sort_values(["Rank", "Probability"], ascending=[True, False])
    result["Query"] = label
    return result


# The network was trained to identify a compound from the transition it
# induces. Therefore the disease signature itself yields phenocopy candidates;
# therapeutic reversal candidates require the sign-inverted disease signature.
phenocopy = predict_one(disease_vscore.rename("IPF_phenocopy"), "Disease_phenocopy")
reversal = predict_one((-disease_vscore).rename("IPF_reversal"), "Disease_reversal")

compound_info = pd.read_csv(ROOT / "00_data" / "LINCS" / "compoundinfo_beta.txt", sep="\t", dtype=str).fillna("")


def collapse_nonempty(values: pd.Series) -> str:
    clean = [x.strip() for x in values.astype(str) if x.strip()]
    return " | ".join(dict.fromkeys(clean))


metadata = (
    compound_info.groupby("pert_id", as_index=False)
    .agg(
        Compound_Name=("cmap_name", collapse_nonempty),
        Target=("target", collapse_nonempty),
        MOA=("moa", collapse_nonempty),
        Canonical_SMILES=("canonical_smiles", collapse_nonempty),
        InChIKey=("inchi_key", collapse_nonempty),
        Aliases=("compound_aliases", collapse_nonempty),
    )
    .rename(columns={"pert_id": "Drug_ID"})
)
for table in (phenocopy, reversal):
    table.reset_index(drop=True, inplace=True)

phenocopy = phenocopy.merge(metadata, how="left", on="Drug_ID")
reversal = reversal.merge(metadata, how="left", on="Drug_ID")

phenocopy.to_csv(OUT / "GSE32537_IPF_DrugReflector_phenocopy_full.csv", index=False)
phenocopy.head(50).to_csv(OUT / "GSE32537_IPF_DrugReflector_phenocopy_top50.csv", index=False)
reversal.to_csv(OUT / "GSE32537_IPF_DrugReflector_reversal_full.csv", index=False)
reversal.head(50).to_csv(OUT / "GSE32537_IPF_DrugReflector_reversal_top50.csv", index=False)
disease_vscore.sort_values(ascending=False).rename("vscore").to_csv(
    OUT / "GSE32537_IPF_input_vscores.csv", header=True
)

found = int(coverage["total_found"])
required = len(model.model.dimensions["var_names"][0])
missing = required - found
run_info = {
    "dataset": "GSE32537",
    "contrast": "IPF/UIP - control",
    "control_samples": 50,
    "ipf_samples": 119,
    "input_genes": int(len(disease_vscore)),
    "model_genes_found": found,
    "model_genes_missing": missing,
    "model_genes_required": required,
    "model_gene_coverage_percent": 100 * found / required,
    "checkpoints": [str(p) for p in checkpoint_paths],
    "n_predictions_requested": 200,
    "method_note": (
        "DrugReflector is a compound-identity classifier. The IPF-minus-control signature "
        "is the phenocopy query; its sign-inverse is the therapeutic reversal query."
    ),
    "interpretation": "Ranked computational hypotheses, not evidence of clinical efficacy.",
}
(OUT / "run_info.json").write_text(json.dumps(run_info, ensure_ascii=False, indent=2), encoding="utf-8")

top = reversal.head(20).sort_values("Probability", ascending=True)
fig, ax = plt.subplots(figsize=(8, 7))
ax.barh(top["Drug_ID"], top["Probability"], color="#2A9D8F")
ax.set_xlabel("DrugReflector probability")
ax.set_ylabel("")
ax.set_title("Top DrugReflector candidates for reversing the GSE32537 IPF signature")
ax.grid(axis="x", alpha=0.25)
fig.tight_layout()
fig.savefig(OUT / "GSE32537_IPF_DrugReflector_top20.png", dpi=300)
fig.savefig(OUT / "GSE32537_IPF_DrugReflector_top20.pdf")
plt.close(fig)

print(json.dumps(run_info, ensure_ascii=False, indent=2))
print("\nTop reversal candidates")
print(reversal.head(20).to_string(index=False))
print("\nTop disease-phenocopy candidates")
print(phenocopy.head(10).to_string(index=False))

