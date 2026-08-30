#!/usr/bin/env python3
"""Matched-sample pseudobulk analysis of published GEO astrocytes and malignant cells."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


TUMOR_STATES = {
    "CP-like", "Cilia-like", "Hypoxic", "IPC-like", "Mesenchymal-like",
    "NPC-like", "OPC-like", "RG-like",
}

MODULES = {
    "Cytokine_growth_factor_secretion": [
        "IL6", "LIF", "OSM", "TGFB1", "TGFB2", "TGFB3", "VEGFA", "VEGFB",
        "FGF2", "HGF", "IGF1", "IGF2", "PDGFA", "PDGFB", "PDGFC", "CXCL12",
        "CCL2", "CSF1", "PTN", "MDK", "NRG1", "KITLG", "SPP1", "SERPINE1", "SHH",
    ],
    "ECM_adhesion_contact": [
        "FN1", "TNC", "SPARCL1", "THBS1", "THBS2", "COL1A1", "COL1A2",
        "COL3A1", "COL4A1", "COL4A2", "LAMA2", "LAMA4", "LAMC1", "VCAN",
        "DCN", "HSPG2", "JAG1", "ICAM1", "VCAM1", "GJA1", "CD44", "LGALS1",
        "MMP2", "MMP9",
    ],
    "Lipid_cholesterol_delivery": [
        "APOE", "CLU", "ABCA1", "ABCG1", "LPL", "FASN", "SCD", "SREBF1",
        "SREBF2", "HMGCR", "HMGCS1", "SQLE", "LDLR", "NPC1", "FABP5",
    ],
    "Lactate_glucose_metabolic_support": [
        "SLC2A1", "SLC2A3", "HK2", "PFKFB3", "GAPDH", "PGK1", "ENO1", "LDHA",
        "SLC16A1", "SLC16A3", "GYS1", "PYGB", "GPI",
    ],
    "Glutamate_glutamine_support": [
        "SLC1A2", "SLC1A3", "GLUL", "GLUD1", "GOT1", "GOT2", "SLC38A2",
        "SLC7A5", "SLC1A4",
    ],
    "Antioxidant_glutathione_support": [
        "SLC7A11", "GCLC", "GCLM", "GSS", "GSR", "GPX1", "GPX4", "GSTP1",
        "NQO1", "HMOX1", "SOD1", "SOD2", "PRDX1", "PRDX2", "NFE2L2",
    ],
    "Extracellular_vesicle_release": [
        "RAB27A", "RAB27B", "SMPD3", "PDCD6IP", "TSG101", "CD63", "CD81",
        "SDCBP", "YKT6", "VAMP7",
    ],
    "Mitochondrial_transfer_permissive": [
        "TNFAIP2", "RHOT1", "RHOT2", "TRAK1", "TRAK2", "KIF5B", "MYO10",
        "CDC42", "RAC1", "RHOA", "MFN1", "MFN2", "OPA1", "GJA1", "CD38", "GAP43",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True, help="Raw gene-by-nucleus count matrix TSV.GZ.")
    parser.add_argument("--metadata", type=Path, required=True, help="GSE283839 metadata TSV or TSV.GZ.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for pseudobulk and module outputs.")
    parser.add_argument("--min-astrocytes", type=int, default=5)
    parser.add_argument("--min-tumor-cells", type=int, default=20)
    return parser.parse_args()


def read_metadata(path: Path) -> pd.DataFrame:
    metadata = pd.read_csv(path, sep="\t", compression="infer", dtype=str)
    required = {"Cell", "ID", "Final_Annotation", "Final_Annotation_Focus"}
    missing = sorted(required.difference(metadata.columns))
    if missing:
        raise ValueError(f"Metadata is missing columns: {missing}")
    metadata["count_cell_id"] = metadata["Cell"].str.replace("-1_", ".1_", n=1, regex=False)
    return metadata


def aggregate_counts(count_path: Path, groups: dict[str, str], wanted_genes: set[str]):
    with gzip.open(count_path, "rt") as handle:
        count_cells = handle.readline().rstrip("\n").split("\t")
        group_labels = sorted(set(groups.values()))
        group_index = {label: index for index, label in enumerate(group_labels)}
        membership = np.array([group_index.get(groups.get(cell), -1) for cell in count_cells], dtype=np.int16)
        keep = membership >= 0
        selected_membership = membership[keep]
        cell_counts = np.bincount(selected_membership, minlength=len(group_labels)).astype(int)
        library_sizes = np.zeros(len(group_labels), dtype=np.int64)
        aggregated: dict[str, np.ndarray] = {}
        duplicates = Counter()

        for line_number, line in enumerate(handle, start=2):
            gene, values = line.rstrip("\n").split("\t", 1)
            vector = np.fromstring(values, sep="\t", dtype=np.int64)
            if len(vector) != len(count_cells):
                raise RuntimeError(f"Malformed count row {line_number}: {gene}")
            group_counts = np.bincount(
                selected_membership,
                weights=vector[keep],
                minlength=len(group_labels),
            ).astype(np.int64)
            library_sizes += group_counts
            if gene in wanted_genes:
                duplicates[gene] += 1
                aggregated[gene] = aggregated.get(gene, np.zeros(len(group_labels), dtype=np.int64)) + group_counts

    counts = pd.DataFrame(aggregated, index=group_labels).T.reindex(sorted(wanted_genes), fill_value=0)
    return (
        counts,
        pd.Series(library_sizes, index=group_labels, name="library_size"),
        pd.Series(cell_counts, index=group_labels, name="n_cells"),
        sorted(gene for gene, occurrences in duplicates.items() if occurrences > 1),
    )


def main() -> None:
    args = parse_args()
    counts_path = args.counts.resolve()
    metadata_path = args.metadata.resolve()
    output_dir = args.output_dir.resolve()
    for path in (counts_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    output_dir.mkdir(parents=True, exist_ok=True)

    metadata = read_metadata(metadata_path)
    metadata["compartment"] = np.select(
        [
            metadata["Final_Annotation"].eq("Astrocytes") & metadata["Final_Annotation_Focus"].eq("TME"),
            metadata["Final_Annotation"].isin(TUMOR_STATES),
        ],
        ["Astrocytes", "Tumor"],
        default="Excluded",
    )

    sample_counts = (
        metadata.loc[metadata["compartment"].ne("Excluded")]
        .groupby(["ID", "compartment"])
        .size()
        .unstack(fill_value=0)
    )
    for compartment in ("Astrocytes", "Tumor"):
        if compartment not in sample_counts:
            sample_counts[compartment] = 0
    sample_counts["eligible"] = (
        sample_counts["Astrocytes"].ge(args.min_astrocytes)
        & sample_counts["Tumor"].ge(args.min_tumor_cells)
    )
    eligible_samples = sample_counts.index[sample_counts["eligible"]].tolist()
    if len(eligible_samples) != 3:
        raise RuntimeError(f"Expected three eligible matched samples, found {len(eligible_samples)}: {eligible_samples}")
    sample_counts.reset_index().to_csv(output_dir / "matched_sample_cell_counts.csv", index=False)

    selected_metadata = metadata.loc[
        metadata["ID"].isin(eligible_samples) & metadata["compartment"].ne("Excluded")
    ]
    groups = {
        row.count_cell_id: f"{row.ID}|{row.compartment}"
        for row in selected_metadata.itertuples()
    }
    wanted_genes = {gene for genes in MODULES.values() for gene in genes}
    counts, library_sizes, group_cell_counts, duplicate_genes = aggregate_counts(
        counts_path, groups, wanted_genes
    )

    counts.to_csv(output_dir / "support_gene_pseudobulk_raw_counts.csv")
    log2_cpm = np.log2(counts.div(library_sizes, axis=1) * 1e6 + 1)
    log2_cpm.to_csv(output_dir / "support_gene_pseudobulk_log2_cpm.csv")

    pseudobulk_metadata = pd.DataFrame({
        "group": library_sizes.index,
        "sample": [group.split("|", 1)[0] for group in library_sizes.index],
        "compartment": [group.split("|", 1)[1] for group in library_sizes.index],
        "library_size": library_sizes.values,
        "n_cells": group_cell_counts.values,
    })
    pseudobulk_metadata.to_csv(output_dir / "pseudobulk_metadata.csv", index=False)

    module_rows = []
    for module, genes in MODULES.items():
        scores = log2_cpm.loc[genes].mean(axis=0)
        for group, score in scores.items():
            sample, compartment = group.split("|", 1)
            module_rows.append({"sample": sample, "compartment": compartment, "module": module, "score": score})
    module_scores = pd.DataFrame(module_rows)
    module_scores.to_csv(output_dir / "astrocyte_support_module_scores.csv", index=False)

    module_wide = module_scores.pivot(index=["sample", "module"], columns="compartment", values="score").reset_index()
    module_wide["astrocyte_minus_tumor"] = module_wide["Astrocytes"] - module_wide["Tumor"]
    module_summary = (
        module_wide.groupby("module", as_index=False)["astrocyte_minus_tumor"]
        .agg(median_difference="median", minimum_difference="min", maximum_difference="max")
    )
    module_wide.to_csv(output_dir / "astrocyte_support_module_paired_differences.csv", index=False)
    module_summary.to_csv(output_dir / "astrocyte_support_module_summary.csv", index=False)

    gene_long = log2_cpm.T.reset_index(names="group").melt(id_vars="group", var_name="gene", value_name="log2_cpm")
    gene_long[["sample", "compartment"]] = gene_long["group"].str.split("|", n=1, expand=True)
    gene_wide = gene_long.pivot(index=["sample", "gene"], columns="compartment", values="log2_cpm").reset_index()
    gene_wide["astrocyte_minus_tumor"] = gene_wide["Astrocytes"] - gene_wide["Tumor"]
    gene_wide.to_csv(output_dir / "support_gene_paired_differences.csv", index=False)

    definitions = pd.DataFrame(
        [(module, gene) for module, genes in MODULES.items() for gene in genes],
        columns=["module", "gene"],
    )
    definitions.to_csv(output_dir / "astrocyte_support_module_gene_sets.csv", index=False)

    report = {
        "eligible_samples": eligible_samples,
        "n_eligible_samples": len(eligible_samples),
        "astrocytes_in_matched_analysis": int(sample_counts.loc[eligible_samples, "Astrocytes"].sum()),
        "minimum_astrocytes": args.min_astrocytes,
        "minimum_tumor_cells": args.min_tumor_cells,
        "tumor_states": sorted(TUMOR_STATES),
        "n_modules": len(MODULES),
        "duplicate_gene_rows_combined": duplicate_genes,
        "interpretation_limit": "Scores measure relative transcriptional compatibility and do not measure metabolite flux, secretion, matrix deposition, vesicle release, or organelle transfer.",
    }
    (output_dir / "analysis_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
