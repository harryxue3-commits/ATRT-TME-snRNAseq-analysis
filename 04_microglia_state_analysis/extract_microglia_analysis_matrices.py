#!/usr/bin/env python3
"""Extract the raw-count matrix and gene-set definitions for microglial scoring."""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


MODULE_GENES = {
    "homeostatic": ["P2RY12", "TMEM119", "SALL1", "HEXB", "CX3CR1", "GPR34", "SLC2A5", "BIN1", "CSF1R", "TREM2"],
    "m1_inflammatory": ["IL1B", "TNF", "IL6", "CXCL9", "CXCL10", "CXCL11", "CCL2", "CCL3", "CCL4", "NFKBIA", "STAT1", "IRF1", "NOS2", "CD86"],
    "m2_immunoregulatory": ["MRC1", "CD163", "MSR1", "IL10", "TGFB1", "ARG1", "CCL18", "MAFB", "VSIG4", "FOLR2", "MARCO"],
    "interferon": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "GBP1", "STAT1", "IRF7"],
    "phagocytic_lipid": ["APOE", "APOC1", "LPL", "SPP1", "GPNMB", "TREM2", "CTSD", "CTSB", "LGALS3", "FABP5", "ABCA1"],
    "antigen_presentation": ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "CIITA", "B2M", "TAP1"],
    "hypoxia_stress": ["HIF1A", "VEGFA", "CA9", "LDHA", "BNIP3", "NDRG1", "SLC2A1", "HMOX1"],
    "complement": ["C1QA", "C1QB", "C1QC", "C3", "CFB", "SERPING1"],
    "proliferation": ["MKI67", "TOP2A", "CENPF", "UBE2C", "TYMS", "PCNA"],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--counts", type=Path, required=True, help="Raw gene-by-nucleus count matrix TSV.GZ.")
    parser.add_argument("--metadata", type=Path, required=True, help="GSE283839 metadata TSV or TSV.GZ.")
    parser.add_argument("--output-dir", type=Path, required=True, help="Directory for extracted microglial inputs.")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open("r")


def normalized_cell_id(cell_id: str) -> str:
    return cell_id.replace(".1_", "-1_", 1)


def read_microglia_ids(metadata_path: Path) -> set[str]:
    with open_text(metadata_path) as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        identifiers = {
            row["Cell"]
            for row in rows
            if row["Final_Annotation"] == "Microglia" and row["Final_Annotation_Focus"] == "TME"
        }
    if len(identifiers) != 288:
        raise RuntimeError(f"Expected 288 microglial nuclei, found {len(identifiers)}")
    return identifiers


def extract_microglia_matrix(counts_path: Path, microglia_ids: set[str], output_path: Path) -> None:
    with gzip.open(counts_path, "rt") as source:
        count_ids = source.readline().rstrip("\n").split("\t")
        metadata_ids = [normalized_cell_id(cell_id) for cell_id in count_ids]
        selected = [index for index, cell_id in enumerate(metadata_ids) if cell_id in microglia_ids]
        if len(selected) != len(microglia_ids):
            raise RuntimeError(f"Matched {len(selected)} of {len(microglia_ids)} microglial identifiers")

        with gzip.open(output_path, "wt") as output:
            output.write("gene\t" + "\t".join(metadata_ids[index] for index in selected) + "\n")
            for line_number, line in enumerate(source, start=2):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != len(count_ids) + 1:
                    raise RuntimeError(f"Malformed count-matrix row {line_number}")
                output.write(fields[0] + "\t" + "\t".join(fields[index + 1] for index in selected) + "\n")


def write_gene_sets(output_path: Path) -> None:
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["module", "gene"])
        for module, genes in MODULE_GENES.items():
            writer.writerows((module, gene) for gene in genes)


def main() -> None:
    args = parse_args()
    counts = args.counts.resolve()
    metadata = args.metadata.resolve()
    output_dir = args.output_dir.resolve()
    for path in (counts, metadata):
        if not path.is_file():
            raise FileNotFoundError(path)
    output_dir.mkdir(parents=True, exist_ok=True)
    microglia_ids = read_microglia_ids(metadata)
    extract_microglia_matrix(counts, microglia_ids, output_dir / "microglia_raw_counts.tsv.gz")
    write_gene_sets(output_dir / "microglia_state_program_gene_sets.tsv")
    print(f"Extracted raw counts for {len(microglia_ids)} microglial nuclei into {output_dir}")


if __name__ == "__main__":
    main()
