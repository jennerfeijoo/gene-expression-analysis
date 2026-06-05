"""Tests for GPL annotation helpers."""

from __future__ import annotations

import pandas as pd

from src.annotation import (
    annotate_probe_ranking,
    collapse_by_largest_abs_effect,
    collapse_to_gene_level,
    compare_gene_collapse_rules,
    standardize_annotation_columns,
    summarize_gene_probe_consistency,
    summarize_probes_per_gene,
)


def test_standardize_annotation_columns_accepts_alternate_names() -> None:
    annotation = pd.DataFrame(
        {
            "ID": ["probe_1"],
            "Gene Symbol": ["GENE1"],
            "Gene Name": ["Example gene"],
            "Entrez Gene ID": ["101"],
        }
    )

    standardized = standardize_annotation_columns(annotation)

    assert standardized.columns.tolist() == [
        "probe_id",
        "gene_symbol",
        "gene_title",
        "entrez_id",
    ]
    assert standardized.loc[0].to_dict() == {
        "probe_id": "probe_1",
        "gene_symbol": "GENE1",
        "gene_title": "Example gene",
        "entrez_id": "101",
    }


def test_annotate_probe_ranking_preserves_all_ranking_rows() -> None:
    ranking = pd.DataFrame(
        {
            "probe_id": ["probe_1", "probe_2"],
            "adjusted_p_value": [0.01, 0.20],
        }
    )
    annotation = pd.DataFrame(
        {
            "ID": ["probe_1"],
            "Gene Symbol": ["GENE1"],
        }
    )

    annotated = annotate_probe_ranking(ranking, annotation)

    assert len(annotated) == len(ranking)
    assert annotated["probe_id"].tolist() == ranking["probe_id"].tolist()
    assert annotated.loc[0, "gene_symbol"] == "GENE1"
    assert pd.isna(annotated.loc[1, "gene_symbol"])


def test_missing_gene_symbols_are_handled_safely() -> None:
    annotation = pd.DataFrame(
        {
            "Probe Set ID": ["probe_1", "probe_2"],
            "Symbol": ["---", ""],
            "Description": ["First gene", "Second gene"],
        }
    )

    standardized = standardize_annotation_columns(annotation)

    assert standardized["gene_symbol"].isna().all()
    assert standardized["gene_title"].tolist() == ["First gene", "Second gene"]


def test_summarize_probes_per_gene_counts_expanded_symbols() -> None:
    annotated = pd.DataFrame(
        {
            "probe_id": ["probe_1", "probe_2", "probe_3"],
            "gene_symbol": ["GENE1", "GENE1 /// GENE2", pd.NA],
        }
    )

    summary = summarize_probes_per_gene(annotated).set_index("gene_symbol")

    assert summary.loc["GENE1", "probe_count"] == 2
    assert summary.loc["GENE2", "probe_count"] == 1


def test_collapse_to_gene_level_keeps_one_row_per_gene() -> None:
    annotated = pd.DataFrame(
        {
            "probe_id": ["probe_1", "probe_2", "probe_3"],
            "gene_symbol": ["GENE1", "GENE1", "GENE2"],
            "adjusted_p_value": [0.05, 0.01, 0.20],
            "mean_paired_difference": [3.0, 1.0, -2.0],
        }
    )

    collapsed = collapse_to_gene_level(annotated)

    assert collapsed["gene_symbol"].nunique() == 2
    assert len(collapsed) == 2
    assert collapsed.set_index("gene_symbol").loc["GENE1", "probe_count"] == 2


def test_collapse_to_gene_level_selects_smallest_adjusted_p_value() -> None:
    annotated = pd.DataFrame(
        {
            "probe_id": ["probe_less_significant", "probe_representative"],
            "gene_symbol": ["GENE1", "GENE1"],
            "adjusted_p_value": [0.05, 0.001],
            "mean_paired_difference": [5.0, 1.0],
        }
    )

    collapsed = collapse_to_gene_level(annotated)

    assert collapsed.loc[0, "probe_id"] == "probe_representative"
    assert collapsed.loc[0, "adjusted_p_value"] == 0.001


def test_summarize_gene_probe_consistency_flags_mixed_directions() -> None:
    annotated = pd.DataFrame(
        {
            "probe_id": ["probe_1", "probe_2", "probe_3"],
            "gene_symbol": ["GENE1", "GENE1", "GENE2"],
            "mean_paired_difference": [2.0, -1.0, 3.0],
        }
    )

    summary = summarize_gene_probe_consistency(annotated).set_index("gene_symbol")

    assert bool(summary.loc["GENE1", "direction_conflict"])
    assert not bool(summary.loc["GENE1", "same_direction"])
    assert bool(summary.loc["GENE2", "same_direction"])


def test_collapse_by_largest_abs_effect_selects_largest_effect_probe() -> None:
    annotated = pd.DataFrame(
        {
            "probe_id": ["probe_significant", "probe_large_effect"],
            "gene_symbol": ["GENE1", "GENE1"],
            "adjusted_p_value": [0.001, 0.05],
            "mean_paired_difference": [1.0, -4.0],
        }
    )

    collapsed = collapse_by_largest_abs_effect(annotated)

    assert collapsed.loc[0, "probe_id"] == "probe_large_effect"
    assert collapsed.loc[0, "probe_count"] == 2


def test_compare_gene_collapse_rules_reports_top_overlap() -> None:
    table_a = pd.DataFrame({"gene_symbol": ["A", "B", "C", "D"]})
    table_b = pd.DataFrame({"gene_symbol": ["B", "C", "E", "F"]})

    comparison = compare_gene_collapse_rules(table_a, table_b, top_n=3)

    assert comparison["in_both"].sum() == 2
    assert set(comparison.loc[comparison["in_both"], "gene_symbol"]) == {"B", "C"}
