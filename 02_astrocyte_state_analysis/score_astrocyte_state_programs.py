"""Conventional rank-based analysis of the 47 GEO-annotated astrocytes.

This is a dependency-light reproduction of the current UCell U-statistic. It
uses raw integer counts, treats the four programs independently, and derives a
descriptive positivity call from expression-matched random signatures. Calls
are not cell-type labels or evidence of biological function.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import rankdata


# Curated, mutually disjoint programs fixed before the present conventional
# rescoring. Functional terms in the names are shorthand for transcriptomic
# evidence and do not imply mutually exclusive cell types or biological effects.
SIGNATURES = {
    "Homeostatic/support": [
        "ALDH1L1", "AQP4", "SLC1A2", "SLC1A3", "GLUL", "GJA1", "S100B",
        "SOX9", "FGFR3", "ATP1A2", "ATP1B2", "KCNJ10", "SLC4A4", "GPC5",
    ],
    "Pan-reactive": [
        "GFAP", "VIM", "CD44", "SERPINA3", "LCN2", "CHI3L1", "OSMR",
        "STAT3", "SOCS3", "JUN", "FOS", "CEBPB", "ICAM1", "VCAM1",
        "TIMP1", "ANXA1",
    ],
    "Inflammatory/IFN/complement-associated": [
        "C3", "SERPING1", "GBP2", "MX1", "IFITM3", "STAT1", "CXCL10",
        "CCL2", "IL1B", "TNF", "NFKBIA", "IRF1", "HLA-DRA", "HLA-DRB1",
        "CD74", "IFIT1", "IFIT3", "ISG15", "CFB", "C1S", "C1R",
    ],
    "ECM/tissue-remodeling-associated": [
        "S100A10", "PTX3", "EMP1", "SPHK1", "CD109", "CLCF1", "TGM1",
        "TM4SF1", "SPARCL1", "CLU", "APOE", "THBS1", "VEGFA",
    ],
}

# Compact, nonredundant marker panel for Figure 2D. These genes are displayed
# but are not given extra weight in the program scores.
DISPLAYED_MARKERS = {
    "Homeostatic/support": ["SLC1A2", "SLC1A3", "AQP4"],
    "Pan-reactive": ["GFAP", "VIM", "CD44"],
    "Inflammatory/IFN/complement-associated": ["STAT1", "IFIT3", "C1S"],
    "ECM/tissue-remodeling-associated": ["SPARCL1", "EMP1", "VEGFA"],
}


def ucell(rank_matrix: np.ndarray, indices: np.ndarray, max_rank: int) -> np.ndarray:
    """Calculate UCell 2.7.6+ normalized Mann-Whitney U scores."""
    n_genes = len(indices)
    min_sum = n_genes * (n_genes + 1) / 2
    max_u = n_genes * max_rank - min_sum
    score = 1 - (rank_matrix[:, indices].sum(axis=1) - min_sum) / max_u
    return np.clip(score, 0, 1)


def bh(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values)
    ordered = values[order]
    adjusted = ordered * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Cell2Sentence-ready TME H5AD containing raw counts in .raw.")
    parser.add_argument("--output-dir", dest="output", type=Path, required=True, help="Directory for astrocyte score outputs.")
    parser.add_argument("--max-rank", type=int, default=1500)
    parser.add_argument("--n-null", type=int, default=10_000)
    parser.add_argument("--pool-size", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20260806)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.input)
    selected = adata.obs["Final_Annotation"].astype(str).eq("Astrocytes")
    astro = adata[selected].copy()
    raw = astro.raw.to_adata()
    counts = raw.X.toarray() if sparse.issparse(raw.X) else np.asarray(raw.X)
    counts = counts.astype(float, copy=False)
    library_size = counts.sum(axis=1)
    log1p_cp10k = np.log1p(
        counts * (10_000 / np.maximum(library_size, 1))[:, None]
    )
    genes = raw.var_names.astype(str).to_numpy()
    gene_index = {gene: index for index, gene in enumerate(genes)}

    missing = {
        program: [gene for gene in signature if gene not in gene_index]
        for program, signature in SIGNATURES.items()
    }
    if any(missing.values()):
        raise RuntimeError(f"Signature genes absent from .raw: {missing}")

    ranks = np.empty_like(counts, dtype=float)
    for cell_index in range(counts.shape[0]):
        ranks[cell_index] = np.minimum(
            rankdata(-counts[cell_index], method="average"), args.max_rank
        )

    tested_genes = {gene for signature in SIGNATURES.values() for gene in signature}
    valid_controls = np.array(
        [
            gene not in tested_genes
            and not gene.startswith(("MT-", "RPL", "RPS"))
            and gene not in {"MALAT1", "NEAT1"}
            for gene in genes
        ]
    )
    valid_indices = np.flatnonzero(valid_controls)
    detection = (counts > 0).mean(axis=0)
    mean_rank = ranks.mean(axis=0)
    detection_scale = np.std(detection[valid_indices]) or 1
    rank_scale = np.std(mean_rank[valid_indices]) or 1
    rng = np.random.default_rng(args.seed)

    records: list[dict] = []
    for program, signature in SIGNATURES.items():
        signature_indices = np.array([gene_index[gene] for gene in signature])
        observed = ucell(ranks, signature_indices, args.max_rank)
        detected = (counts[:, signature_indices] > 0).sum(axis=1)

        matching_pools = []
        for index in signature_indices:
            distance = (
                ((detection[valid_indices] - detection[index]) / detection_scale) ** 2
                + ((mean_rank[valid_indices] - mean_rank[index]) / rank_scale) ** 2
            )
            nearest = valid_indices[
                np.argsort(distance, kind="stable")[: args.pool_size]
            ]
            matching_pools.append(nearest)

        null_scores = np.empty((args.n_null, astro.n_obs), dtype=float)
        for iteration in range(args.n_null):
            chosen: list[int] = []
            used: set[int] = set()
            for pool in matching_pools:
                available = pool[~np.isin(pool, np.fromiter(used, dtype=int))]
                selected_control = int(rng.choice(available if len(available) else pool))
                chosen.append(selected_control)
                used.add(selected_control)
            null_scores[iteration] = ucell(
                ranks, np.asarray(chosen, dtype=int), args.max_rank
            )

        null_p95 = np.quantile(null_scores, 0.95, axis=0, method="higher")
        empirical_p = (
            1 + (null_scores >= observed).sum(axis=0)
        ) / (args.n_null + 1)
        positive = (observed > null_p95) & (detected >= 2)

        for index, cell_id in enumerate(astro.obs_names):
            records.append(
                {
                    "cell_id": str(cell_id),
                    "sample_id": str(astro.obs.iloc[index]["ID"]),
                    "subtype": str(astro.obs.iloc[index]["subtype"]),
                    "program": program,
                    "signature_size": len(signature),
                    "ucell": observed[index],
                    "detected_signature_genes": int(detected[index]),
                    "null_p95": null_p95[index],
                    "empirical_p": empirical_p[index],
                    "positive": bool(positive[index]),
                }
            )

    cells = pd.DataFrame(records)
    cells["q_within_cell"] = cells.groupby("cell_id")["empirical_p"].transform(
        lambda values: bh(values.to_numpy())
    )
    cells["q_global"] = bh(cells["empirical_p"].to_numpy())

    samples = (
        cells.groupby(["sample_id", "subtype", "program"], observed=True)
        .agg(
            n_cells=("cell_id", "nunique"),
            mean_ucell=("ucell", "mean"),
            median_ucell=("ucell", "median"),
            positive_cells=("positive", "sum"),
            positive_fraction=("positive", "mean"),
            mean_detected_genes=("detected_signature_genes", "mean"),
        )
        .reset_index()
    )
    summaries = (
        cells.groupby("program", observed=True)
        .agg(
            signature_size=("signature_size", "first"),
            mean_ucell=("ucell", "mean"),
            median_ucell=("ucell", "median"),
            positive_cells=("positive", "sum"),
            positive_fraction=("positive", "mean"),
            mean_detected_genes=("detected_signature_genes", "mean"),
        )
        .reset_index()
    )

    program_order = list(SIGNATURES)
    positive_wide = cells.pivot(index="cell_id", columns="program", values="positive")
    combinations = []
    for _, row in positive_wide[program_order].iterrows():
        selected_programs = [program for program in program_order if bool(row[program])]
        combinations.append(" + ".join(selected_programs) if selected_programs else "None")
    upset = (
        pd.Series(combinations)
        .value_counts()
        .rename_axis("combination")
        .reset_index(name="n_cells")
    )
    upset["percent"] = 100 * upset["n_cells"] / astro.n_obs

    gene_records = []
    sample_gene_records = []
    for program, signature in SIGNATURES.items():
        for gene in signature:
            index = gene_index[gene]
            detected_mask = counts[:, index] > 0
            gene_records.append(
                {
                    "program": program,
                    "gene": gene,
                    "detected_cells": int(detected_mask.sum()),
                    "detection_percent": 100 * detected_mask.mean(),
                    "mean_raw_count": counts[:, index].mean(),
                    "mean_log1p_cp10k_all_cells": log1p_cp10k[:, index].mean(),
                    "mean_log1p_cp10k_detected_cells": (
                        log1p_cp10k[detected_mask, index].mean()
                        if detected_mask.any()
                        else 0
                    ),
                }
            )
            for sample_id in astro.obs["ID"].astype(str).unique():
                sample_mask = astro.obs["ID"].astype(str).to_numpy() == sample_id
                sample_detected = detected_mask[sample_mask]
                sample_gene_records.append(
                    {
                        "sample_id": sample_id,
                        "subtype": str(astro.obs.loc[sample_mask, "subtype"].iloc[0]),
                        "n_cells": int(sample_mask.sum()),
                        "program": program,
                        "gene": gene,
                        "detected_cells": int(sample_detected.sum()),
                        "detection_percent": 100 * sample_detected.mean(),
                        "mean_log1p_cp10k_all_cells": log1p_cp10k[
                            sample_mask, index
                        ].mean(),
                    }
                )

    cells.to_csv(args.output / "cell_program_scores.csv", index=False)
    samples.to_csv(args.output / "sample_program_summary.csv", index=False)
    summaries.to_csv(args.output / "program_summary.csv", index=False)
    upset.to_csv(args.output / "upset_combinations.csv", index=False)
    pd.DataFrame(gene_records).to_csv(
        args.output / "signature_gene_detection.csv", index=False
    )
    pd.DataFrame(sample_gene_records).to_csv(
        args.output / "sample_marker_expression.csv", index=False
    )

    displayed_markers = pd.DataFrame(
        [
            {
                "program": program,
                "gene": gene,
                "selection_rule": "representative_nonredundant",
            }
            for program, marker_genes in DISPLAYED_MARKERS.items()
            for gene in marker_genes
        ]
    )
    displayed_markers.to_csv(args.output / "displayed_markers.csv", index=False)
    pd.DataFrame(sample_gene_records).merge(
        displayed_markers[["program", "gene"]],
        on=["program", "gene"],
        how="inner",
        validate="many_to_one",
    ).to_csv(args.output / "marker_dotplot_by_sample.csv", index=False)
    (args.output / "analysis_report.json").write_text(
        json.dumps(
            {
                "input": str(args.input),
                "n_cells": int(astro.n_obs),
                "n_tumor_samples": int(astro.obs["ID"].nunique()),
                "tumor_sample_cell_counts": {
                    str(key): int(value)
                    for key, value in astro.obs["ID"].value_counts().sort_index().items()
                },
                "signatures": SIGNATURES,
                "displayed_markers": DISPLAYED_MARKERS,
                "max_rank": args.max_rank,
                "n_null": args.n_null,
                "matching_pool_size": args.pool_size,
                "seed": args.seed,
                "positivity": (
                    "observed UCell > cell-specific 95th percentile of matched null "
                    "and >=2 detected signature genes"
                ),
            },
            indent=2,
        )
    )
    print(summaries.to_string(index=False))
    print(f"\nOutputs: {args.output}")


if __name__ == "__main__":
    main()
