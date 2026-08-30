#!/usr/bin/env python3
"""Generate ranked Cell2Sentence gene sentences from a prepared AnnData file."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import pandas as pd
from cell2sentence import utils


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Cell2Sentence-ready H5AD file.")
    parser.add_argument("--output", type=Path, required=True, help="Output CSV containing cell_name and cell_sentence.")
    parser.add_argument("--top-genes", type=int, default=1000, help="Maximum ranked genes retained per cell.")
    parser.add_argument("--seed", type=int, default=42, help="Seed used by Cell2Sentence to resolve expression ties.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input.resolve()
    output_path = args.output.resolve()
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    if args.top_genes < 1:
        raise ValueError("--top-genes must be positive")

    adata = ad.read_h5ad(input_path)
    if not adata.obs_names.is_unique:
        raise ValueError("Cell identifiers are not unique")
    if not adata.var_names.is_unique:
        raise ValueError("Gene identifiers are not unique")

    vocabulary = utils.generate_vocabulary(adata)
    full_sentences = utils.generate_sentences(
        adata=adata,
        vocab=vocabulary,
        delimiter=" ",
        random_state=args.seed,
    )
    ranked_sentences = [" ".join(sentence.split()[: args.top_genes]) for sentence in full_sentences]

    output = pd.DataFrame({
        "cell_name": adata.obs_names.astype(str),
        "cell_sentence": ranked_sentences,
        "n_genes_in_sentence": [len(sentence.split()) for sentence in ranked_sentences],
    })
    for column in ("Final_Annotation", "Final_Annotation_Focus", "ID", "subtype"):
        if column in adata.obs.columns:
            output[column] = adata.obs[column].astype(str).to_numpy()

    if output["cell_name"].duplicated().any() or len(output) != adata.n_obs:
        raise RuntimeError("Output cell identifiers failed validation")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    print(f"Wrote {len(output):,} ranked cell sentences to {output_path}")


if __name__ == "__main__":
    main()
