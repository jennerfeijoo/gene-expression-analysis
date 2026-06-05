"""Tests for GPL annotation helpers."""

from __future__ import annotations

import pandas as pd

from src.annotation import annotate_probe_ranking, standardize_annotation_columns


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
