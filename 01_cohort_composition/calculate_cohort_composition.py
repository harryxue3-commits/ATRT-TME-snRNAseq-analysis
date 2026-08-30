#!/usr/bin/env python3
"""Freeze and audit the source data used in manuscript Figure 1.

The script uses the deposited GSE283839 metadata without additional filtering.
The published Harmony coordinates are carried forward without recomputing
an embedding.  ``Final_Annotation_Focus == "TME"`` is the sole rule for the
nonmalignant TME compartment; all other published focus labels are tumor labels.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kruskal


METADATA: Path
COUNT_MATRIX: Path
SOURCE_DATA: Path

EXPECTED_TOTAL = 17_564
EXPECTED_GENES = 36_601
EXPECTED_TME = 516
EXPECTED_TUMOR = 17_048
EXPECTED_TME_COUNTS = {
    "Microglia": 288,
    "OPC": 79,
    "Endothelial": 48,
    "Astrocytes": 47,
    "Neurons": 34,
    "Pericytes": 20,
}

TME_ORDER = ["Microglia", "OPC", "Endothelial", "Astrocytes", "Neurons", "Pericytes"]
SUBTYPE_ORDER = ["ATRT-MYC", "ATRT-SHH", "ATRT-TYR"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def pct(numerator: pd.Series | float, denominator: pd.Series | float):
    return 100.0 * numerator / denominator


def audit_matrix_header(metadata_ids: pd.Series) -> dict[str, object]:
    with gzip.open(COUNT_MATRIX, "rt") as handle:
        matrix_ids = handle.readline().rstrip("\n").split("\t")
        gene_rows = sum(1 for _ in handle)

    converted = [cell_id.replace(".1_", "-1_", 1) for cell_id in matrix_ids]
    metadata_set = set(metadata_ids)
    matrix_set = set(converted)
    duplicate_matrix_ids = len(converted) - len(matrix_set)
    duplicate_metadata_ids = int(metadata_ids.duplicated().sum())
    missing_from_metadata = sorted(matrix_set - metadata_set)
    missing_from_matrix = sorted(metadata_set - matrix_set)

    if len(matrix_ids) != EXPECTED_TOTAL:
        raise ValueError(f"Count-matrix header has {len(matrix_ids):,} cells; expected {EXPECTED_TOTAL:,}.")
    if gene_rows != EXPECTED_GENES:
        raise ValueError(f"Count matrix has {gene_rows:,} gene rows; expected {EXPECTED_GENES:,}.")
    if duplicate_matrix_ids or duplicate_metadata_ids or missing_from_metadata or missing_from_matrix:
        raise ValueError("Count-matrix and metadata identifiers do not form a unique one-to-one match.")

    return {
        "matrix_cell_columns": len(matrix_ids),
        "matrix_gene_rows": gene_rows,
        "matrix_shape_genes_by_nuclei": f"{gene_rows} x {len(matrix_ids)}",
        "metadata_rows": len(metadata_ids),
        "barcode_transform": "replace the first '.1_' with '-1_'",
        "matched_identifiers": len(matrix_set & metadata_set),
        "duplicate_matrix_identifiers_after_transform": duplicate_matrix_ids,
        "duplicate_metadata_identifiers": duplicate_metadata_ids,
        "missing_from_metadata": len(missing_from_metadata),
        "missing_from_matrix": len(missing_from_matrix),
        "one_to_one_match": True,
    }


def write_csv(frame: pd.DataFrame, name: str) -> None:
    frame.to_csv(SOURCE_DATA / name, index=False, float_format="%.9g")


def build() -> None:
    data = pd.read_csv(METADATA, sep="\t")
    required = {
        "Cell", "subtype", "ID", "Final_Annotation", "Final_Annotation_Focus",
        "UMAPHARMONY_1", "UMAPHARMONY_2",
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Metadata is missing required columns: {missing}")
    if len(data) != EXPECTED_TOTAL:
        raise ValueError(f"Metadata has {len(data):,} rows; expected {EXPECTED_TOTAL:,}.")
    if data[sorted(required)].isna().any().any():
        raise ValueError("Required Figure 1 metadata fields contain missing values.")
    if data["Cell"].duplicated().any():
        raise ValueError("Metadata cell identifiers are not unique.")

    data = data.copy()
    data.insert(0, "source_row", np.arange(1, len(data) + 1))
    data["compartment"] = np.where(data["Final_Annotation_Focus"].eq("TME"), "TME", "Tumor")
    data["resolved_immune_class"] = np.select(
        [
            data["Final_Annotation"].eq("Microglia") & data["compartment"].eq("TME"),
            data["compartment"].eq("TME"),
        ],
        ["Resolved innate immune (Microglia)", "Non-immune TME"],
        default="Not applicable (tumor)",
    )

    compartment_counts = data["compartment"].value_counts().to_dict()
    if compartment_counts != {"Tumor": EXPECTED_TUMOR, "TME": EXPECTED_TME}:
        raise ValueError(f"Unexpected compartment counts: {compartment_counts}")
    observed_tme = (
        data.loc[data["compartment"].eq("TME"), "Final_Annotation"]
        .value_counts()
        .to_dict()
    )
    if observed_tme != EXPECTED_TME_COUNTS:
        raise ValueError(f"Unexpected TME annotation counts: {observed_tme}")

    matrix_audit = audit_matrix_header(data["Cell"])

    cell_columns = [
        "source_row", "Cell", "ID", "subtype", "Final_Annotation",
        "Final_Annotation_Focus", "compartment", "resolved_immune_class",
        "UMAPHARMONY_1", "UMAPHARMONY_2",
    ]
    data[cell_columns].to_csv(
        SOURCE_DATA / "Figure_1_cell_level_metadata.tsv.gz",
        sep="\t",
        index=False,
        float_format="%.9g",
        compression={"method": "gzip", "mtime": 0},
    )

    overall = (
        data.groupby("compartment", observed=True)
        .size()
        .rename("nuclei")
        .reset_index()
    )
    overall["percentage_of_all_nuclei"] = pct(overall["nuclei"], len(data))
    overall["all_nuclei_denominator"] = len(data)
    overall["figure_panel"] = "C"
    overall["definition"] = overall["compartment"].map({
        "Tumor": "Final_Annotation_Focus != TME; includes published Unannotated tumor nuclei",
        "TME": "Final_Annotation_Focus == TME",
    })
    overall["compartment_order"] = overall["compartment"].map({"Tumor": 1, "TME": 2})
    overall = overall.sort_values("compartment_order").drop(columns="compartment_order")
    write_csv(overall, "Figure_1_overall_compartments.csv")

    tme = data.loc[data["compartment"].eq("TME")].copy()
    tme_composition = (
        tme.groupby("Final_Annotation", observed=True)
        .size()
        .reindex(TME_ORDER, fill_value=0)
        .rename("nuclei")
        .reset_index()
        .rename(columns={"Final_Annotation": "tme_cell_type"})
    )
    tme_composition["percentage_of_TME"] = pct(tme_composition["nuclei"], len(tme))
    tme_composition["TME_denominator"] = len(tme)
    tme_composition["percentage_of_all_nuclei"] = pct(tme_composition["nuclei"], len(data))
    tme_composition["all_nuclei_denominator"] = len(data)
    tme_composition["figure_panel"] = "D"
    write_csv(tme_composition, "Figure_1_TME_composition.csv")

    immune_audit = pd.DataFrame([
        {
            "lineage": "Microglia / CNS-resident myeloid",
            "published_cell_type_annotation": "Microglia",
            "resolution_status": "Separately resolved published TME label",
            "nuclei": 288,
            "percentage_of_TME": 100 * 288 / EXPECTED_TME,
            "percentage_of_all_nuclei": 100 * 288 / EXPECTED_TOTAL,
            "interpretation": "Resolved innate immune lineage",
        },
        {
            "lineage": "T cells",
            "published_cell_type_annotation": pd.NA,
            "resolution_status": "Not separately resolved by published TME labels",
            "nuclei": pd.NA,
            "percentage_of_TME": pd.NA,
            "percentage_of_all_nuclei": pd.NA,
            "interpretation": "Do not encode as zero or infer biological absence",
        },
        {
            "lineage": "B / plasma cells",
            "published_cell_type_annotation": pd.NA,
            "resolution_status": "Not separately resolved by published TME labels",
            "nuclei": pd.NA,
            "percentage_of_TME": pd.NA,
            "percentage_of_all_nuclei": pd.NA,
            "interpretation": "Do not encode as zero or infer biological absence",
        },
        {
            "lineage": "Natural killer cells",
            "published_cell_type_annotation": pd.NA,
            "resolution_status": "Not separately resolved by published TME labels",
            "nuclei": pd.NA,
            "percentage_of_TME": pd.NA,
            "percentage_of_all_nuclei": pd.NA,
            "interpretation": "Do not encode as zero or infer biological absence",
        },
    ])
    immune_audit["figure_panel"] = "D"
    write_csv(immune_audit, "Figure_1_immune_annotation_audit.csv")

    sample_counts = (
        data.groupby(["ID", "subtype", "compartment"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .rename(columns={"ID": "sample_id", "TME": "TME_nuclei", "Tumor": "tumor_nuclei"})
    )
    for name in ["TME_nuclei", "tumor_nuclei"]:
        if name not in sample_counts:
            sample_counts[name] = 0
    sample_counts["total_nuclei"] = sample_counts["TME_nuclei"] + sample_counts["tumor_nuclei"]
    sample_counts["TME_percentage"] = pct(sample_counts["TME_nuclei"], sample_counts["total_nuclei"])
    sample_counts["tumor_percentage"] = pct(sample_counts["tumor_nuclei"], sample_counts["total_nuclei"])
    sample_counts["subtype_order"] = sample_counts["subtype"].map(dict(zip(SUBTYPE_ORDER, range(3))))
    sample_counts["sample_number"] = sample_counts["sample_id"].str.extract(r"-(\d+)$").astype(int)
    sample_counts = sample_counts.sort_values(["subtype_order", "sample_number"])
    sample_counts["figure_panels"] = "E;F"
    sample_counts = sample_counts.drop(columns=["subtype_order", "sample_number"])
    write_csv(sample_counts, "Figure_1_sample_compartments.csv")

    sample_tme = (
        tme.groupby(["ID", "subtype", "Final_Annotation"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(columns=TME_ORDER, fill_value=0)
        .reset_index()
        .rename(columns={"ID": "sample_id"})
    )
    sample_tme["TME_nuclei"] = sample_tme[TME_ORDER].sum(axis=1)
    sample_tme["figure_panel"] = "F"
    for cell_type in TME_ORDER:
        sample_tme[f"{cell_type}_percentage_of_TME"] = pct(sample_tme[cell_type], sample_tme["TME_nuclei"])
    sample_tme["subtype_order"] = sample_tme["subtype"].map(dict(zip(SUBTYPE_ORDER, range(3))))
    sample_tme["sample_number"] = sample_tme["sample_id"].str.extract(r"-(\d+)$").astype(int)
    sample_tme = sample_tme.sort_values(["subtype_order", "sample_number"]).drop(columns=["subtype_order", "sample_number"])
    write_csv(sample_tme, "Figure_1_sample_TME_composition.csv")

    pooled_subtype = (
        data.groupby(["subtype", "compartment"], observed=True)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .rename(columns={"TME": "TME_nuclei", "Tumor": "tumor_nuclei"})
    )
    pooled_subtype["total_nuclei"] = pooled_subtype["TME_nuclei"] + pooled_subtype["tumor_nuclei"]
    pooled_subtype["pooled_TME_percentage"] = pct(pooled_subtype["TME_nuclei"], pooled_subtype["total_nuclei"])

    subtype_summary = (
        sample_counts.groupby("subtype", observed=True)["TME_percentage"]
        .agg(samples="count", median_TME_percentage="median", mean_TME_percentage="mean", minimum_TME_percentage="min", maximum_TME_percentage="max")
        .reset_index()
        .merge(pooled_subtype, on="subtype", how="left")
    )
    subtype_summary["figure_panel"] = "E"
    subtype_summary["subtype_order"] = subtype_summary["subtype"].map(dict(zip(SUBTYPE_ORDER, range(3))))
    subtype_summary = subtype_summary.sort_values("subtype_order").drop(columns="subtype_order")
    write_csv(subtype_summary, "Figure_1_subtype_summary.csv")

    groups = [
        sample_counts.loc[sample_counts["subtype"].eq(subtype), "TME_percentage"].to_numpy()
        for subtype in SUBTYPE_ORDER
    ]
    full_test = kruskal(*groups)
    sensitivity = sample_counts.loc[~sample_counts["sample_id"].eq("ATRT-MYC-2")].copy()
    sensitivity_groups = [
        sensitivity.loc[sensitivity["subtype"].eq(subtype), "TME_percentage"].to_numpy()
        for subtype in SUBTYPE_ORDER
    ]
    sensitivity_test = kruskal(*sensitivity_groups)
    statistics = pd.DataFrame([
        {
            "analysis": "Primary 12-sample subtype comparison",
            "outcome": "Within-sample TME percentage",
            "test": "Kruskal-Wallis",
            "biological_unit": "Tumor sample",
            "groups": "ATRT-MYC n=3; ATRT-SHH n=5; ATRT-TYR n=4",
            "H": full_test.statistic,
            "df": 2,
            "P_value": full_test.pvalue,
            "multiplicity_correction": "Not applied; one prespecified Figure 1 comparison",
            "result": "No subtype-associated difference detected",
            "figure_panel": "E",
        },
        {
            "analysis": "Sensitivity excluding extreme low-n sample ATRT-MYC-2",
            "outcome": "Within-sample TME percentage",
            "test": "Kruskal-Wallis",
            "biological_unit": "Tumor sample",
            "groups": "ATRT-MYC n=2; ATRT-SHH n=5; ATRT-TYR n=4",
            "H": sensitivity_test.statistic,
            "df": 2,
            "P_value": sensitivity_test.pvalue,
            "multiplicity_correction": "Not applied; prespecified sensitivity analysis",
            "result": "No subtype-associated difference detected",
            "figure_panel": "E",
        },
    ])
    write_csv(statistics, "Figure_1_statistics.csv")

    cell_type_by_sample = (
        data.groupby(["ID", "subtype", "Final_Annotation", "compartment", "resolved_immune_class"], observed=True)
        .size()
        .rename("cell_count")
        .reset_index()
        .rename(columns={"ID": "sample_id", "Final_Annotation": "published_cell_type_annotation"})
        .merge(sample_counts[["sample_id", "total_nuclei", "TME_nuclei"]], on="sample_id", how="left")
    )
    cell_type_by_sample["percentage_of_sample"] = pct(cell_type_by_sample["cell_count"], cell_type_by_sample["total_nuclei"])
    cell_type_by_sample["percentage_of_sample_TME"] = np.where(
        cell_type_by_sample["compartment"].eq("TME"),
        pct(cell_type_by_sample["cell_count"], cell_type_by_sample["TME_nuclei"]),
        np.nan,
    )
    write_csv(cell_type_by_sample, "Supplementary_Table_S1_Cell_Type_Counts_and_Proportions.csv")

    label_map = (
        data.groupby(["Final_Annotation", "Final_Annotation_Focus", "compartment", "resolved_immune_class"], observed=True)
        .size()
        .rename("nuclei")
        .reset_index()
        .rename(columns={"Final_Annotation": "published_cell_type_annotation", "Final_Annotation_Focus": "published_annotation_focus"})
    )
    label_map["mapping_rule"] = np.where(
        label_map["compartment"].eq("TME"),
        "Final_Annotation_Focus == TME",
        "Final_Annotation_Focus != TME",
    )
    label_map["interpretive_note"] = ""
    label_map.loc[label_map["published_cell_type_annotation"].eq("Unannotated"), "interpretive_note"] = (
        "Published tumor nuclei without a specific tumor-state signature"
    )
    label_map.loc[label_map["published_cell_type_annotation"].eq("Immune-like"), "interpretive_note"] = (
        "Published tumor state; not an infiltrating immune-cell label"
    )
    label_map.loc[label_map["published_cell_type_annotation"].eq("OPC"), "interpretive_note"] = (
        "Nonmalignant TME OPC; distinct from the OPC-like tumor state"
    )
    label_map.loc[label_map["published_cell_type_annotation"].eq("OPC-like"), "interpretive_note"] = (
        "Tumor state; distinct from nonmalignant TME OPC"
    )
    write_csv(label_map, "Figure_1_label_map.csv")

    file_inventory = pd.DataFrame([
        {
            "role": "Figure 1 metadata and published Harmony coordinates",
            "file_name": METADATA.name,
            "path": str(METADATA),
            "bytes": METADATA.stat().st_size,
            "sha256": sha256(METADATA),
        },
        {
            "role": "Raw gene-by-nucleus count matrix; header used for cell-ID audit",
            "file_name": COUNT_MATRIX.name,
            "path": str(COUNT_MATRIX),
            "bytes": COUNT_MATRIX.stat().st_size,
            "sha256": sha256(COUNT_MATRIX),
        },
    ])
    file_inventory["accessed_date"] = "2026-07-15"
    write_csv(file_inventory, "Figure_1_source_file_inventory.csv")

    audit = {
        "cohort_freeze_date": "2026-07-15",
        "dataset": "GSE283839 ATRT RNA V3 subseries",
        "additional_filtering_for_figure_1": "None; deposited processed cohort used as supplied",
        "total_nuclei": len(data),
        "samples": int(data["ID"].nunique()),
        "subtype_sample_counts": data[["ID", "subtype"]].drop_duplicates()["subtype"].value_counts().sort_index().to_dict(),
        "published_harmony_coordinates": ["UMAPHARMONY_1", "UMAPHARMONY_2"],
        "coordinate_missing_values": int(data[["UMAPHARMONY_1", "UMAPHARMONY_2"]].isna().sum().sum()),
        "compartment_rule": "Final_Annotation_Focus == 'TME' -> TME; otherwise -> Tumor",
        "tumor_nuclei": EXPECTED_TUMOR,
        "TME_nuclei": EXPECTED_TME,
        "TME_cell_type_counts": EXPECTED_TME_COUNTS,
        "resolved_immune_annotation": (
            "Microglia are the only separately resolved immune lineage in the published TME labels; "
            "adaptive lymphocyte populations were not separately resolved, so their biological absence is not inferred."
        ),
        "matrix_metadata_identifier_audit": matrix_audit,
        "primary_subtype_test": {
            "test": "Kruskal-Wallis",
            "H": float(full_test.statistic),
            "df": 2,
            "P_value": float(full_test.pvalue),
            "biological_unit": "tumor sample",
        },
        "sensitivity_excluding_ATRT_MYC_2": {
            "reason": "102 total nuclei; 89 published TME nuclei",
            "H": float(sensitivity_test.statistic),
            "df": 2,
            "P_value": float(sensitivity_test.pvalue),
        },
    }
    (SOURCE_DATA / "Figure_1_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")

    readme = f"""# Figure 1 frozen source data

Figure 1 uses the deposited processed GSE283839 ATRT RNA V3 cohort without additional filtering: **{len(data):,} nuclei from 12 tumor samples** (ATRT-MYC n=3, ATRT-SHH n=5, ATRT-TYR n=4). The UMAP panels use the supplied `UMAPHARMONY_1/2` coordinates; no embedding was recomputed.

## Classification rule

- `Final_Annotation_Focus == TME` -> nonmalignant TME ({EXPECTED_TME:,} nuclei; {100*EXPECTED_TME/len(data):.3f}%).
- All other published focus labels -> tumor ({EXPECTED_TUMOR:,} nuclei; {100*EXPECTED_TUMOR/len(data):.3f}%).
- `Unannotated` denotes published tumor nuclei without a specific tumor-state signature.
- `Immune-like` is a published tumor state, not an infiltrating immune-cell label.
- `OPC` is a nonmalignant TME label; `OPC-like` is a tumor-state label.

The published TME labels separately resolve microglia but not T, B, or NK populations. The absence of separate adaptive-immune labels is therefore not treated as evidence of biological absence. Figure 1 supports the wording **sparse, microglia-dominated published TME**, not a definitive immune-cold classification.

## Statistical unit

Subtype comparisons use each tumor sample as the biological unit. The primary comparison of within-sample TME percentage used one Kruskal-Wallis test (H={full_test.statistic:.3f}, df=2, P={full_test.pvalue:.3f}); no multiplicity adjustment was required for this single prespecified Figure 1 comparison. Excluding the extreme low-n ATRT-MYC-2 sample gave H={sensitivity_test.statistic:.3f}, P={sensitivity_test.pvalue:.3f}.

## Files

- `Figure_1_cell_level_metadata.tsv.gz`: plotted cell identifiers, published labels, frozen compartment map, and Harmony coordinates.
- `Figure_1_overall_compartments.csv`: source for panel C.
- `Figure_1_TME_composition.csv`: source for panel D.
- `Figure_1_immune_annotation_audit.csv`: explicit resolved-versus-unresolved immune-lineage audit for panel D interpretation.
- `Figure_1_sample_compartments.csv` and `Figure_1_subtype_summary.csv`: source for panel E.
- `Figure_1_sample_TME_composition.csv`: source for panel F.
- `Figure_1_statistics.csv`: primary and sensitivity tests.
- `Figure_1_label_map.csv`: auditable label semantics and collapse.
- `Supplementary_Table_S1_Cell_Type_Counts_and_Proportions.csv`: sample-level long-form source table.
- `Figure_1_source_file_inventory.csv`: input names, sizes, access date, and SHA-256 hashes.
- `Figure_1_audit.json`: machine-readable cohort and identifier audit.
"""
    (SOURCE_DATA / "README_Figure_1_source_data.md").write_text(readme, encoding="utf-8")

    print(f"Wrote frozen Figure 1 source data to {SOURCE_DATA}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True, help="GSE283839 metadata TSV.")
    parser.add_argument("--counts", type=Path, required=True, help="Raw gene-by-nucleus count matrix TSV.GZ.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for cohort-composition outputs.")
    return parser.parse_args()


def main() -> None:
    global METADATA, COUNT_MATRIX, SOURCE_DATA
    args = parse_args()
    METADATA = args.metadata.resolve()
    COUNT_MATRIX = args.counts.resolve()
    SOURCE_DATA = args.output_dir.resolve()
    for path in (METADATA, COUNT_MATRIX):
        if not path.is_file():
            raise FileNotFoundError(path)
    SOURCE_DATA.mkdir(parents=True, exist_ok=True)
    build()


if __name__ == "__main__":
    main()
