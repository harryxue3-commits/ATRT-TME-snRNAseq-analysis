#!/usr/bin/env python3
"""Build a Cell2Sentence-ready AnnData file from the published GSE283839 TME cells.

Published nonmalignant TME cells have Final_Annotation_Focus == "TME".
The count matrix uses a period in the 10x suffix (BARCODE.1...), whereas the
metadata uses the standard hyphen (BARCODE-1...). Only this one substitution is
allowed, and the script requires a complete, unambiguous one-to-one match.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


BARCODE_DOT_RE = re.compile(r"^([ACGTN]+)\.1(?=_|$)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--counts",
        type=Path,
        required=True,
        help="Gene-by-cell raw-count matrix (TSV.GZ).",
    )
    parser.add_argument(
        "--metadata",
        type=Path,
        required=True,
        help="Published GEO cell metadata (TSV.GZ).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        required=True,
        help="Output AnnData H5AD file.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=None,
        help="Machine-readable build and validation summary.",
    )
    parser.add_argument(
        "--chunk-genes",
        type=int,
        default=500,
        help="Number of gene rows read per chunk from the wide count matrix.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def count_id_to_metadata_id(cell_id: str) -> str:
    converted, replacements = BARCODE_DOT_RE.subn(r"\1-1", cell_id, count=1)
    if replacements != 1:
        raise ValueError(f"Unexpected count-matrix cell ID format: {cell_id!r}")
    return converted


def log10_normalize_total(counts: sparse.csr_matrix) -> tuple[sparse.csr_matrix, float]:
    """Replicate Scanpy normalize_total(target_sum=None) + log1p(base=10)."""
    library_sizes = np.asarray(counts.sum(axis=1)).ravel().astype(np.float64)
    positive = library_sizes[library_sizes > 0]
    if positive.size != counts.shape[0]:
        raise ValueError("At least one selected TME cell has zero total counts.")
    target_sum = float(np.median(positive))
    scale = (target_sum / library_sizes).astype(np.float32)
    normalized = counts.astype(np.float32).multiply(scale[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data) / np.log(10.0)
    normalized = normalized.astype(np.float32)
    normalized.eliminate_zeros()
    return normalized, target_sum


def main() -> None:
    args = parse_args()
    if args.summary is None:
        args.summary = args.output.with_name(args.output.stem + "_summary.json")
    if args.chunk_genes < 1:
        raise ValueError("--chunk-genes must be positive.")
    print("Selection: published GEO nonmalignant TME cells", flush=True)
    for path in (args.counts, args.metadata):
        if not path.is_file():
            raise FileNotFoundError(f"Required input file not found: {path}")

    metadata = pd.read_csv(args.metadata, sep="\t", compression="infer")
    required = {"Cell", "Final_Annotation", "Final_Annotation_Focus"}
    missing_columns = sorted(required.difference(metadata.columns))
    if missing_columns:
        raise ValueError(f"Metadata is missing required columns: {missing_columns}")
    if metadata["Cell"].isna().any() or metadata["Cell"].duplicated().any():
        raise ValueError("Metadata Cell IDs contain missing values or duplicates.")

    metadata["Cell"] = metadata["Cell"].astype(str)
    metadata_by_cell = metadata.set_index("Cell", drop=False)
    tme_mask = metadata["Final_Annotation_Focus"].astype(str).eq("TME")
    tme_metadata = metadata.loc[tme_mask].copy()
    if tme_metadata.empty:
        raise ValueError("No published GEO TME cells were found.")
    if tme_metadata["Cell"].duplicated().any():
        raise ValueError("The published GEO TME subset contains duplicate Cell IDs.")
    print(f"Selected metadata cells: {len(tme_metadata):,}", flush=True)

    count_columns = pd.read_csv(args.counts, sep="\t", compression="gzip", nrows=0).columns.astype(str).tolist()
    converted_count_ids = [count_id_to_metadata_id(cell_id) for cell_id in count_columns]
    if len(set(converted_count_ids)) != len(converted_count_ids):
        raise ValueError("Count cell IDs become duplicated after the period-to-hyphen conversion.")
    count_to_metadata = dict(zip(count_columns, converted_count_ids, strict=True))
    metadata_to_count = {metadata_id: count_id for count_id, metadata_id in count_to_metadata.items()}

    count_id_set = set(converted_count_ids)
    metadata_id_set = set(metadata["Cell"])
    if count_id_set != metadata_id_set:
        only_counts = sorted(count_id_set.difference(metadata_id_set))
        only_metadata = sorted(metadata_id_set.difference(count_id_set))
        raise ValueError(
            "Count/metadata cell IDs do not match completely after the documented conversion. "
            f"Only counts: {len(only_counts)}; only metadata: {len(only_metadata)}."
        )

    tme_ids = tme_metadata["Cell"].tolist()
    missing_tme = sorted(set(tme_ids).difference(metadata_to_count))
    if missing_tme:
        raise ValueError(f"Published GEO TME cells missing from counts: {len(missing_tme)}")
    selected_count_columns = {metadata_to_count[cell_id] for cell_id in tme_ids}

    # Read the extremely wide matrix in gene-row chunks to avoid materializing a
    # dense cells-by-genes array. Pandas recognizes the first, unlabeled field in
    # each data row as the gene index.
    reader = pd.read_csv(
        args.counts,
        sep="\t",
        compression="gzip",
        usecols=lambda column: column in selected_count_columns,
        chunksize=args.chunk_genes,
    )
    sparse_gene_blocks: list[sparse.csr_matrix] = []
    gene_names: list[str] = []
    ordered_count_columns: list[str] | None = None
    for chunk_number, chunk in enumerate(reader, start=1):
        if ordered_count_columns is None:
            ordered_count_columns = chunk.columns.astype(str).tolist()
            if len(ordered_count_columns) != len(tme_ids):
                raise ValueError(
                    f"Expected {len(tme_ids)} selected count columns, "
                    f"read {len(ordered_count_columns)}."
                )
        elif chunk.columns.astype(str).tolist() != ordered_count_columns:
            raise ValueError("Selected count columns changed between gene chunks.")
        if chunk.index.isna().any():
            raise ValueError("The count matrix contains missing gene identifiers.")
        if not all(pd.api.types.is_integer_dtype(dtype) for dtype in chunk.dtypes):
            values64 = chunk.to_numpy(dtype=np.float64, copy=True)
            if not np.all(np.isfinite(values64)) or not np.all(values64 == np.floor(values64)):
                raise ValueError("The raw-count matrix contains non-integer or non-finite values.")
        values = chunk.to_numpy(dtype=np.int32, copy=True)
        if (values < 0).any():
            raise ValueError("The raw-count matrix contains negative values.")
        gene_names.extend(chunk.index.astype(str).tolist())
        sparse_gene_blocks.append(sparse.csr_matrix(values.T, dtype=np.int32))
        if chunk_number % 10 == 0:
            print(f"Read {len(gene_names):,} gene rows", flush=True)

    if ordered_count_columns is None or not sparse_gene_blocks:
        raise ValueError("No expression data were read from the count matrix.")
    raw_counts_all = sparse.hstack(sparse_gene_blocks, format="csr", dtype=np.int32)
    del sparse_gene_blocks
    print(
        f"Sparse raw matrix assembled: {raw_counts_all.shape[0]:,} cells x "
        f"{raw_counts_all.shape[1]:,} genes; {raw_counts_all.nnz:,} nonzero counts",
        flush=True,
    )
    genes_all = pd.Index(gene_names, name="gene_symbol")
    duplicate_gene_rows = int(genes_all.duplicated(keep=False).sum())
    if genes_all.has_duplicates:
        # Aggregate duplicate symbols through a sparse projection so C2S sees
        # each gene token exactly once without densifying the matrix.
        codes, unique_genes = pd.factorize(genes_all, sort=False)
        projection = sparse.csr_matrix(
            (np.ones(len(codes), dtype=np.int32), (np.arange(len(codes)), codes)),
            shape=(len(codes), len(unique_genes)),
        )
        raw_counts_all = (raw_counts_all @ projection).tocsr().astype(np.int32)
        genes_all = pd.Index(unique_genes.astype(str), name="gene_symbol")

    ordered_metadata_ids = [count_to_metadata[cell_id] for cell_id in ordered_count_columns]
    if len(set(ordered_metadata_ids)) != len(ordered_metadata_ids):
        raise ValueError("Selected cell mapping is ambiguous after matrix loading.")
    if set(ordered_metadata_ids) != set(tme_ids):
        raise ValueError("Loaded count columns do not exactly equal the published GEO TME set.")

    obs = metadata_by_cell.loc[ordered_metadata_ids].copy()
    obs.index = pd.Index(ordered_metadata_ids, name="cell_id")
    obs["count_matrix_cell_id"] = ordered_count_columns
    obs["published_tme"] = obs["Final_Annotation_Focus"].astype(str).eq("TME")
    obs["published_cell_type"] = obs["Final_Annotation"].astype(str)
    # Convenient aliases used in the Cell2Sentence tutorials.
    obs["cell_type"] = obs["Final_Annotation"].astype(str)
    obs["batch_condition"] = obs["ID"].astype(str)
    obs["organism"] = "Homo sapiens"
    obs["tissue"] = "ATRT tumor microenvironment"
    obs["assay"] = obs["technology"].astype(str)

    detected_genes_per_cell = np.asarray((raw_counts_all > 0).sum(axis=1)).ravel()
    cell_keep = detected_genes_per_cell >= 200
    cells_removed_min_genes = int((~cell_keep).sum())
    if cells_removed_min_genes:
        raw_counts_all = raw_counts_all[cell_keep].tocsr()
        obs = obs.iloc[np.flatnonzero(cell_keep)].copy()
        ordered_metadata_ids = obs.index.astype(str).tolist()
        detected_genes_per_cell = detected_genes_per_cell[cell_keep]

    n_cells_all = np.asarray((raw_counts_all > 0).sum(axis=0)).ravel().astype(np.int32)
    total_counts_all = np.asarray(raw_counts_all.sum(axis=0)).ravel().astype(np.int64)
    gene_keep = n_cells_all >= 3
    if not gene_keep.any():
        raise ValueError("No genes remain after the Cell2Sentence min_cells=3 filter.")

    counts_filtered = raw_counts_all[:, gene_keep].tocsr().astype(np.int32)
    genes_filtered = genes_all[gene_keep]
    normalized, normalization_target_sum = log10_normalize_total(counts_filtered)
    print(f"C2S matrix prepared: {normalized.shape[0]:,} x {normalized.shape[1]:,}", flush=True)

    raw_library_sizes = np.asarray(raw_counts_all.sum(axis=1)).ravel().astype(np.int64)
    filtered_library_sizes = np.asarray(counts_filtered.sum(axis=1)).ravel().astype(np.int64)
    obs["raw_library_size_recomputed"] = raw_library_sizes
    obs["raw_detected_genes_recomputed"] = detected_genes_per_cell.astype(np.int32)
    obs["c2s_filtered_library_size"] = filtered_library_sizes

    var_all = pd.DataFrame(index=genes_all)
    var_all["gene_name"] = genes_all.astype(str)
    var_all["n_cells_raw"] = n_cells_all
    var_all["total_counts_raw"] = total_counts_all
    var_all["c2s_retained_min_3_cells"] = gene_keep
    var_filtered = var_all.loc[genes_filtered].copy()

    adata = ad.AnnData(X=normalized, obs=obs, var=var_filtered)
    adata.layers["counts"] = counts_filtered
    adata.raw = ad.AnnData(X=raw_counts_all, obs=obs.copy(), var=var_all)
    adata.uns["dataset"] = "GSE283839"
    adata.uns["subset_definition"] = "Final_Annotation_Focus == 'TME'"
    adata.uns["subset_interpretation"] = "Published GEO nonmalignant tumor-microenvironment cells"
    adata.uns["cell_id_transformation"] = (
        "Count matrix only: replace the first period in BARCODE.1[_suffix] with "
        "a hyphen, yielding BARCODE-1[_suffix]; required complete one-to-one matching."
    )
    adata.uns["X_contents"] = "normalize_total(target_sum=median library size) followed by log1p(base=10)"
    adata.uns["counts_layer_contents"] = "Raw integer counts for C2S-retained genes"
    adata.uns["raw_contents"] = "Complete raw integer counts for all genes before min_cells filtering"
    adata.uns["c2s_cell_filter"] = "min_genes=200"
    adata.uns["c2s_gene_filter"] = "min_cells=3"
    adata.uns["normalization_target_sum"] = normalization_target_sum
    adata.uns["source_counts_file"] = str(args.counts.resolve())
    adata.uns["source_metadata_file"] = str(args.metadata.resolve())
    adata.uns["source_counts_sha256"] = sha256(args.counts)
    adata.uns["source_metadata_sha256"] = sha256(args.metadata)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    print(f"Writing {args.output}", flush=True)
    adata.write_h5ad(args.output, compression="gzip")

    # Reopen the final artifact and validate its essential invariants.
    check = ad.read_h5ad(args.output)
    print("Reopened H5AD for validation", flush=True)
    if check.n_obs != len(obs) or check.n_vars != int(gene_keep.sum()):
        raise RuntimeError("Reopened H5AD dimensions differ from the constructed object.")
    if set(check.obs_names) != set(ordered_metadata_ids):
        raise RuntimeError("Reopened H5AD cell IDs differ from the validated TME IDs.")
    if "counts" not in check.layers or check.raw is None:
        raise RuntimeError("Reopened H5AD is missing the raw counts layer or .raw matrix.")
    if check.raw.n_vars != len(genes_all):
        raise RuntimeError("Reopened H5AD .raw matrix does not contain all genes.")
    if not np.allclose(check.X.data, normalized.data, rtol=1e-6, atol=1e-7):
        raise RuntimeError("Reopened normalized expression values failed validation.")

    annotation_counts = (
        obs["Final_Annotation"].astype(str).value_counts().sort_index().astype(int).to_dict()
    )
    sample_counts = obs["ID"].astype(str).value_counts().sort_index().astype(int).to_dict()
    summary = {
        "output_h5ad": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
        "output_bytes": args.output.stat().st_size,
        "source_counts": str(args.counts.resolve()),
        "source_metadata": str(args.metadata.resolve()),
        "metadata_cells_total": int(len(metadata)),
        "selected_cells_requested": int(len(tme_ids)),
        "selected_cells_written": int(check.n_obs),
        "published_tme_cells_written": int(obs["published_tme"].sum()),
        "cells_removed_by_c2s_min_genes_200": cells_removed_min_genes,
        "all_count_cells_matched_to_metadata": len(converted_count_ids),
        "id_match_rate": 1.0,
        "duplicate_gene_rows_aggregated": duplicate_gene_rows,
        "genes_raw_all": int(check.raw.n_vars),
        "genes_c2s_retained_min_3_cells": int(check.n_vars),
        "raw_integer_counts_preserved": True,
        "x_is_c2s_base10_log_normalized": True,
        "normalization_target_sum": normalization_target_sum,
        "published_cell_type_counts": annotation_counts,
        "sample_counts": sample_counts,
        "metadata_columns_preserved": metadata.columns.astype(str).tolist(),
        "derived_obs_columns": [
            "count_matrix_cell_id",
            "published_tme",
            "published_cell_type",
            "cell_type",
            "batch_condition",
            "organism",
            "tissue",
            "assay",
            "raw_library_size_recomputed",
            "raw_detected_genes_recomputed",
            "c2s_filtered_library_size",
        ],
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
