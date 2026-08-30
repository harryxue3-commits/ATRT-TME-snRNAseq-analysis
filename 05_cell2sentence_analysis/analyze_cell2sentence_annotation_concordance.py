#!/usr/bin/env python3
"""Build Cell2Sentence supplementary analyses, figures, and source tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


OUT: Path
PREDICTIONS: Path
PROGRAM_SCORES: Path
C2S_RANK_SCORES: Path
ANALYSIS_REPORT: Path

GEO_ORDER = [
    "Astrocyte",
    "Microglia/myeloid",
    "OPC/oligodendroglial",
    "Endothelial",
    "Pericyte/vascular",
    "Neuronal",
]
C2S_ORDER = GEO_ORDER + ["Other immune", "Other/ambiguous"]
PROGRAM_ORDER = [
    "Homeostatic/support",
    "Pan-reactive",
    "Inflammatory/IFN/complement-associated",
    "ECM/tissue-remodeling-associated",
]
PROGRAM_SHORT = {
    "Homeostatic/support": "Homeostatic",
    "Pan-reactive": "Reactive",
    "Inflammatory/IFN/complement-associated": "Inflammatory",
    "ECM/tissue-remodeling-associated": "ECM remodeling",
}
GEO_MAP = {
    "Astrocytes": "Astrocyte",
    "Microglia": "Microglia/myeloid",
    "OPC": "OPC/oligodendroglial",
    "Endothelial": "Endothelial",
    "Pericytes": "Pericyte/vascular",
    "Neurons": "Neuronal",
}


def harmonize_c2s(raw_label: str) -> str:
    """Map raw C2S labels to prespecified broad lineages for descriptive comparison."""
    label = str(raw_label).strip().lower()
    if "astrocyte" in label:
        return "Astrocyte"
    if "oligodendro" in label:
        return "OPC/oligodendroglial"
    if label in {
        "microglial cell",
        "mature microglial cell",
        "macrophage",
        "myeloid cell",
        "central nervous system macrophage",
        "monocyte",
        "dendritic cell",
        "langerhans cell",
    }:
        return "Microglia/myeloid"
    if "endothelial" in label:
        return "Endothelial"
    if label in {
        "pericyte",
        "vascular associated smooth muscle cell",
        "vascular leptomeningeal cell",
    }:
        return "Pericyte/vascular"
    if "neuron" in label:
        return "Neuronal"
    if any(
        term in label
        for term in (
            "t cell",
            "b cell",
            "natural killer",
            "innate lymphoid",
            "plasma cell",
            "leukocyte",
            "megakaryocyte",
        )
    ):
        return "Other immune"
    return "Other/ambiguous"


def build_tables() -> dict[str, pd.DataFrame]:
    predictions = pd.read_csv(PREDICTIONS)
    predictions["GEO_label_broad"] = predictions["Final_Annotation"].map(GEO_MAP)
    predictions["C2S_label_broad"] = predictions["predicted_cell_type"].map(harmonize_c2s)
    predictions["broad_label_concordant"] = (
        predictions["GEO_label_broad"] == predictions["C2S_label_broad"]
    )
    predictions["top_50_ranked_genes"] = predictions["cell_sentence"].fillna("").map(
        lambda value: " ".join(str(value).split()[:50])
    )

    all_tme = predictions[
        [
            "cell_name",
            "ID",
            "subtype",
            "Final_Annotation",
            "GEO_label_broad",
            "predicted_cell_type",
            "C2S_label_broad",
            "broad_label_concordant",
            "n_genes_in_sentence",
            "top_50_ranked_genes",
        ]
    ].rename(
        columns={
            "cell_name": "cell_id",
            "ID": "sample_id",
            "Final_Annotation": "GEO_label_original",
            "predicted_cell_type": "C2S_label_original",
        }
    )

    confusion = pd.crosstab(
        predictions["GEO_label_broad"], predictions["C2S_label_broad"]
    ).reindex(index=GEO_ORDER, columns=C2S_ORDER, fill_value=0)
    confusion_long = confusion.stack().rename("cell_count").reset_index()
    confusion_long["row_total"] = confusion_long.groupby("GEO_label_broad")[
        "cell_count"
    ].transform("sum")
    confusion_long["row_fraction"] = (
        confusion_long["cell_count"] / confusion_long["row_total"]
    )

    metrics = []
    for label in GEO_ORDER:
        geo_positive = predictions["GEO_label_broad"].eq(label)
        c2s_positive = predictions["C2S_label_broad"].eq(label)
        tp = int((geo_positive & c2s_positive).sum())
        fp = int((~geo_positive & c2s_positive).sum())
        fn = int((geo_positive & ~c2s_positive).sum())
        precision = tp / (tp + fp) if tp + fp else np.nan
        recall = tp / (tp + fn) if tp + fn else np.nan
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else np.nan
        metrics.append(
            {
                "broad_lineage": label,
                "true_positive_using_GEO_as_reference": tp,
                "false_positive_using_GEO_as_reference": fp,
                "false_negative_using_GEO_as_reference": fn,
                "precision_descriptive": precision,
                "recall_same_broad_lineage": recall,
                "f1_descriptive": f1,
                "GEO_cell_count": int(geo_positive.sum()),
                "C2S_cell_count": int(c2s_positive.sum()),
            }
        )
    metrics = pd.DataFrame(metrics)
    overall = pd.DataFrame(
        [
            {
                "metric": "Broad-lineage agreement",
                "numerator": int(predictions["broad_label_concordant"].sum()),
                "denominator": len(predictions),
                "fraction": predictions["broad_label_concordant"].mean(),
                "definition": "GEO_label_broad equals C2S_label_broad for all 516 published TME nuclei; no cells excluded.",
            }
        ]
    )

    scores = pd.read_csv(PROGRAM_SCORES)
    score_wide = scores.pivot(
        index=["cell_id", "sample_id", "subtype"],
        columns="program",
        values=[
            "ucell",
            "detected_signature_genes",
            "null_p95",
            "empirical_p",
            "positive",
            "q_within_cell",
            "q_global",
        ],
    )
    score_wide.columns = [
        f"{PROGRAM_SHORT[program]}_{measure}" for measure, program in score_wide.columns
    ]
    score_wide = score_wide.reset_index()
    positive_cols = [f"{PROGRAM_SHORT[program]}_positive" for program in PROGRAM_ORDER]
    score_wide["program_count"] = score_wide[positive_cols].sum(axis=1).astype(int)

    def membership(row: pd.Series) -> str:
        selected = [PROGRAM_SHORT[p] for p in PROGRAM_ORDER if bool(row[f"{PROGRAM_SHORT[p]}_positive"])]
        return " + ".join(selected) if selected else "No program"

    score_wide["exact_program_membership"] = score_wide.apply(membership, axis=1)
    astro_predictions = predictions[predictions["Final_Annotation"].eq("Astrocytes")].copy()
    astro = score_wide.merge(
        astro_predictions[
            [
                "cell_name",
                "predicted_cell_type",
                "C2S_label_broad",
                "n_genes_in_sentence",
                "top_50_ranked_genes",
            ]
        ],
        left_on="cell_id",
        right_on="cell_name",
        how="left",
        validate="one_to_one",
    ).drop(columns="cell_name")
    astro = astro.rename(columns={"predicted_cell_type": "C2S_label_original"})
    astro["C2S_astrocyte_like"] = astro["C2S_label_broad"].eq("Astrocyte")

    astro_group = (
        astro.assign(
            conventional_group=np.where(
                astro["program_count"].eq(0), "No program", "At least one program"
            )
        )
        .groupby(["conventional_group", "C2S_label_broad"], observed=True)
        .size()
        .rename("cell_count")
        .reset_index()
    )

    astro_long = scores.merge(
        astro[["cell_id", "C2S_astrocyte_like", "C2S_label_original", "program_count"]],
        on="cell_id",
        how="left",
        validate="many_to_one",
    )
    astro_long["within_program_z"] = astro_long.groupby("program")["ucell"].transform(
        lambda values: (values - values.mean()) / values.std(ddof=0)
        if values.std(ddof=0) > 0
        else 0
    )

    rank_scores = pd.read_csv(C2S_RANK_SCORES).rename(columns={"Cell": "cell_id"})
    harmonization_map = (
        predictions[["predicted_cell_type", "C2S_label_broad"]]
        .drop_duplicates()
        .sort_values(["C2S_label_broad", "predicted_cell_type"])
        .rename(columns={"predicted_cell_type": "C2S_label_original"})
    )
    harmonization_map["mapping_note"] = (
        "Broad lineage used only for descriptive annotation concordance; not a new cell identity."
    )

    report = json.loads(ANALYSIS_REPORT.read_text(encoding="utf-8"))
    signature_rows = []
    for program, genes in report["signatures"].items():
        for gene in genes:
            signature_rows.append({"program": program, "gene_symbol": gene})

    return {
        "all_tme": all_tme,
        "confusion_long": confusion_long,
        "metrics": metrics,
        "overall": overall,
        "astro": astro,
        "astro_group": astro_group,
        "astro_long": astro_long,
        "rank_scores": rank_scores,
        "harmonization_map": harmonization_map,
        "signatures": pd.DataFrame(signature_rows),
    }


def build_figure(tables: dict[str, pd.DataFrame]) -> None:
    confusion_long = tables["confusion_long"]
    metrics = tables["metrics"]
    overall = tables["overall"].iloc[0]
    astro_group = tables["astro_group"]
    astro_long = tables["astro_long"]

    confusion = confusion_long.pivot(
        index="GEO_label_broad", columns="C2S_label_broad", values="cell_count"
    ).reindex(index=GEO_ORDER, columns=C2S_ORDER, fill_value=0)
    row_pct = confusion.div(confusion.sum(axis=1), axis=0) * 100
    annotations = confusion.astype(str) + "\n" + row_pct.round(0).astype(int).astype(str) + "%"

    mpl.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10,
            "axes.titleweight": "bold",
            "axes.edgecolor": "#98A4AC",
            "axes.linewidth": 0.7,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(11.0, 8.0))
    grid = fig.add_gridspec(
        2,
        2,
        left=0.075,
        right=0.965,
        bottom=0.11,
        top=0.955,
        width_ratios=[1.45, 1.0],
        height_ratios=[1.05, 1.0],
        wspace=0.33,
        hspace=0.58,
    )
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 0])
    ax_d = fig.add_subplot(grid[1, 1])

    cmap = sns.light_palette("#356A8A", as_cmap=True)
    sns.heatmap(
        row_pct,
        annot=annotations,
        fmt="",
        cmap=cmap,
        vmin=0,
        vmax=100,
        linewidths=0.7,
        linecolor="white",
        cbar_kws={"label": "Within-GEO-label percentage"},
        annot_kws={"fontsize": 7},
        ax=ax_a,
    )
    ax_a.set_xlabel("")
    ax_a.set_ylabel("Published GEO annotation")
    ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=35, ha="right")
    ax_a.set_yticklabels(ax_a.get_yticklabels(), rotation=0)
    ax_a.set_title("All 516 TME nuclei: counts and row percentages", loc="left")
    ax_a.text(-0.10, 1.08, "A", transform=ax_a.transAxes, fontsize=15, fontweight="bold")

    metrics_plot = metrics.set_index("broad_lineage").reindex(GEO_ORDER).reset_index()
    y = np.arange(len(metrics_plot))
    ax_b.barh(y, metrics_plot["recall_same_broad_lineage"] * 100, color="#3A9D8F")
    ax_b.set_yticks(y, metrics_plot["broad_lineage"])
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 100)
    ax_b.set_xlabel("Assigned to the same broad lineage (%)")
    ax_b.grid(axis="x", color="#E6EBEE", linewidth=0.7)
    ax_b.set_axisbelow(True)
    for index, row in metrics_plot.iterrows():
        ax_b.text(
            min(98, row["recall_same_broad_lineage"] * 100 + 2),
            index,
            f"{int(row['true_positive_using_GEO_as_reference'])}/{int(row['GEO_cell_count'])}",
            va="center",
            fontsize=7.5,
        )
    ax_b.axvline(float(overall["fraction"]) * 100, color="#B44747", linewidth=1.2, linestyle="--")
    ax_b.text(
        float(overall["fraction"]) * 100,
        len(metrics_plot) - 0.15,
        f"Overall {int(overall['numerator'])}/{int(overall['denominator'])} ({float(overall['fraction']):.1%})",
        ha="right",
        va="bottom",
        fontsize=7.5,
        color="#8F3636",
    )
    ax_b.set_title("Broad-lineage concordance by GEO label", loc="left")
    ax_b.text(-0.18, 1.08, "B", transform=ax_b.transAxes, fontsize=15, fontweight="bold")
    ax_b.spines[["top", "right"]].set_visible(False)

    broad_display = ["Astrocyte", "OPC/oligodendroglial", "Other"]
    astro_group = astro_group.copy()
    astro_group["display_group"] = np.where(
        astro_group["C2S_label_broad"].isin(["Astrocyte", "OPC/oligodendroglial"]),
        astro_group["C2S_label_broad"],
        "Other",
    )
    astro_group = (
        astro_group.groupby(["conventional_group", "display_group"])["cell_count"]
        .sum()
        .unstack(fill_value=0)
        .reindex(index=["No program", "At least one program"], columns=broad_display, fill_value=0)
    )
    colors = {
        "Astrocyte": "#356A8A",
        "OPC/oligodendroglial": "#B8902F",
        "Other": "#A6ADB3",
    }
    left = np.zeros(len(astro_group))
    for category in broad_display:
        values = astro_group[category].to_numpy()
        ax_c.barh(np.arange(len(astro_group)), values, left=left, color=colors[category], label=category)
        for row_index, value in enumerate(values):
            if value:
                ax_c.text(left[row_index] + value / 2, row_index, str(int(value)), ha="center", va="center", color="white" if category == "Astrocyte" else "#26323A", fontweight="bold", fontsize=8)
        left += values
    ax_c.set_yticks(np.arange(len(astro_group)), ["Below all four thresholds\n(n = 23)", "At least one program\n(n = 24)"])
    ax_c.invert_yaxis()
    ax_c.set_xlim(0, 24)
    ax_c.set_xlabel("GEO-annotated astrocytes")
    ax_c.set_title("Astrocyte C2S labels by conventional program evidence", loc="left")
    ax_c.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    ax_c.spines[["top", "right", "left"]].set_visible(False)
    ax_c.grid(axis="x", color="#E6EBEE", linewidth=0.7)
    ax_c.set_axisbelow(True)
    ax_c.text(-0.10, 1.08, "C", transform=ax_c.transAxes, fontsize=15, fontweight="bold")

    plot = astro_long.copy()
    plot["program_short"] = plot["program"].map(PROGRAM_SHORT)
    plot["C2S group"] = np.where(plot["C2S_astrocyte_like"], "Astrocyte-like", "Other label")
    sns.boxplot(
        data=plot,
        x="program_short",
        y="within_program_z",
        hue="C2S group",
        order=[PROGRAM_SHORT[p] for p in PROGRAM_ORDER],
        hue_order=["Astrocyte-like", "Other label"],
        palette={"Astrocyte-like": "#356A8A", "Other label": "#B8902F"},
        width=0.62,
        fliersize=0,
        linewidth=0.8,
        ax=ax_d,
    )
    sns.stripplot(
        data=plot,
        x="program_short",
        y="within_program_z",
        hue="C2S group",
        order=[PROGRAM_SHORT[p] for p in PROGRAM_ORDER],
        hue_order=["Astrocyte-like", "Other label"],
        dodge=True,
        palette={"Astrocyte-like": "#356A8A", "Other label": "#B8902F"},
        alpha=0.55,
        size=2.5,
        linewidth=0,
        ax=ax_d,
    )
    handles, labels = ax_d.get_legend_handles_labels()
    ax_d.legend(handles[:2], ["C2S astrocyte-like (n = 30)", "Other C2S label (n = 17)"], frameon=False, fontsize=7, loc="upper center", bbox_to_anchor=(0.5, -0.22))
    ax_d.axhline(0, color="#9AA5AD", linewidth=0.7)
    ax_d.set_xlabel("")
    ax_d.set_ylabel("Within-program standardized score")
    ax_d.set_xticks(
        np.arange(4),
        ["Homeostatic", "Reactive", "Inflammatory", "ECM"],
        rotation=25,
        ha="right",
    )
    ax_d.set_title("Conventional scores by C2S astrocyte prediction", loc="left")
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.grid(axis="y", color="#E6EBEE", linewidth=0.7)
    ax_d.set_axisbelow(True)
    ax_d.text(-0.16, 1.08, "D", transform=ax_d.transAxes, fontsize=15, fontweight="bold")

    fig.savefig(OUT / "Supplementary_Figure_S1_Cell2Sentence.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(OUT / "Supplementary_Figure_S1_Cell2Sentence.pdf", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True, help="Cell2Sentence cell-level predictions and ranked gene sentences CSV.")
    parser.add_argument("--program-scores", type=Path, required=True, help="Conventional astrocyte program-score CSV.")
    parser.add_argument("--rank-scores", type=Path, required=True, help="Cell2Sentence astrocyte rank-score CSV.")
    parser.add_argument("--analysis-report", type=Path, required=True, help="Conventional astrocyte analysis report JSON.")
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    global OUT, PREDICTIONS, PROGRAM_SCORES, C2S_RANK_SCORES, ANALYSIS_REPORT
    args = parse_args()
    PREDICTIONS = args.predictions.resolve()
    PROGRAM_SCORES = args.program_scores.resolve()
    C2S_RANK_SCORES = args.rank_scores.resolve()
    ANALYSIS_REPORT = args.analysis_report.resolve()
    OUT = args.output_dir.resolve()
    for path in (PREDICTIONS, PROGRAM_SCORES, C2S_RANK_SCORES, ANALYSIS_REPORT):
        if not path.is_file():
            raise FileNotFoundError(path)
    OUT.mkdir(parents=True, exist_ok=True)
    tables = build_tables()
    filenames = {
        "all_tme": "all_tme_cell_level.csv",
        "confusion_long": "confusion_counts.csv",
        "metrics": "agreement_metrics.csv",
        "overall": "overall_agreement.csv",
        "astro": "astrocyte_cell_level.csv",
        "astro_group": "astrocyte_c2s_by_program_group.csv",
        "astro_long": "astrocyte_program_scores_long.csv",
        "rank_scores": "c2s_rank_weighted_astrocyte_scores.csv",
        "harmonization_map": "c2s_label_harmonization.csv",
        "signatures": "conventional_signature_definitions.csv",
    }
    for key, filename in filenames.items():
        tables[key].to_csv(OUT / filename, index=False)
    build_figure(tables)
    summary = {
        "n_all_tme": int(len(tables["all_tme"])),
        "broad_agreement_n": int(tables["overall"].iloc[0]["numerator"]),
        "broad_agreement_denominator": int(tables["overall"].iloc[0]["denominator"]),
        "broad_agreement_fraction": float(tables["overall"].iloc[0]["fraction"]),
        "n_geo_astrocytes": int(len(tables["astro"])),
        "n_c2s_astrocyte_like": int(tables["astro"]["C2S_astrocyte_like"].sum()),
        "no_program_n": int(tables["astro"]["program_count"].eq(0).sum()),
        "no_program_c2s_astrocyte_like_n": int(
            (
                tables["astro"]["program_count"].eq(0)
                & tables["astro"]["C2S_astrocyte_like"]
            ).sum()
        ),
        "program_positive_n": int(tables["astro"]["program_count"].gt(0).sum()),
        "program_positive_c2s_astrocyte_like_n": int(
            (
                tables["astro"]["program_count"].gt(0)
                & tables["astro"]["C2S_astrocyte_like"]
            ).sum()
        ),
    }
    (OUT / "analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
