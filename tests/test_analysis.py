"""Tests for exploratory group-comparison helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from src.analysis import (
    benjamini_hochberg,
    build_paired_sample_table,
    derive_pair_ids,
    paired_ttest_by_probe,
    rank_group_comparison,
)


def test_benjamini_hochberg_preserves_shape_and_rank_monotonicity() -> None:
    p_values = np.array([0.04, 0.001, 0.03, 0.20])

    adjusted = benjamini_hochberg(p_values)
    order = np.argsort(p_values)

    assert adjusted.shape == p_values.shape
    assert np.all(np.diff(adjusted[order]) >= 0)
    assert np.all(adjusted >= p_values)
    assert np.all(adjusted <= 1)


def test_rank_group_comparison_returns_expected_columns() -> None:
    expression = pd.DataFrame(
        {
            "ID_REF": ["probe_1", "probe_2"],
            "A1": [8.0, 3.0],
            "A2": [9.0, 3.2],
            "A3": [10.0, 2.8],
            "B1": [1.0, 3.1],
            "B2": [2.0, 2.9],
            "B3": [1.5, 3.0],
        }
    )
    metadata = pd.DataFrame(
        {
            "sample_accession": ["A1", "A2", "A3", "B1", "B2", "B3"],
            "group": ["tumor", "tumor", "tumor", "normal", "normal", "normal"],
        }
    )

    result = rank_group_comparison(
        expression,
        metadata,
        group_column="group",
        group_a="tumor",
        group_b="normal",
    )

    assert result.columns.tolist() == [
        "probe_id",
        "mean_tumor",
        "mean_normal",
        "log2_fold_change",
        "p_value",
        "adjusted_p_value",
    ]
    assert result["probe_id"].tolist()[0] == "probe_1"
    assert result.loc[0, "log2_fold_change"] > 0


def test_build_paired_sample_table_keeps_complete_pairs() -> None:
    metadata = pd.DataFrame(
        {
            "sample_accession": ["N1", "T1", "N2", "T2"],
            "group": ["normal", "tumor", "normal", "tumor"],
            "individual_id": ["P1", "P1", "P2", "P2"],
        }
    )

    metadata = derive_pair_ids(metadata)
    pair_table = build_paired_sample_table(metadata, "group", "pair_id")

    assert pair_table.columns.tolist() == ["pair_id", "normal", "tumor"]
    assert pair_table["pair_id"].tolist() == ["P1", "P2"]
    assert pair_table["tumor"].tolist() == ["T1", "T2"]


def test_build_paired_sample_table_drops_incomplete_pairs() -> None:
    metadata = pd.DataFrame(
        {
            "sample_accession": ["N1", "T1", "N2"],
            "group": ["normal", "tumor", "normal"],
            "pair_id": ["P1", "P1", "P2"],
        }
    )

    pair_table = build_paired_sample_table(metadata, "group", "pair_id")

    assert pair_table["pair_id"].tolist() == ["P1"]


def test_build_paired_sample_table_rejects_duplicate_group_samples() -> None:
    metadata = pd.DataFrame(
        {
            "sample_accession": ["N1", "N1_repeat", "T1"],
            "group": ["normal", "normal", "tumor"],
            "pair_id": ["P1", "P1", "P1"],
        }
    )

    with pytest.raises(ValueError, match="at most one"):
        build_paired_sample_table(metadata, "group", "pair_id")


def test_paired_ttest_by_probe_returns_expected_columns() -> None:
    expression = pd.DataFrame(
        {
            "ID_REF": ["probe_1", "probe_2"],
            "N1": [1.0, 5.0],
            "T1": [3.0, 5.2],
            "N2": [2.0, 4.8],
            "T2": [5.0, 4.7],
            "N3": [1.5, 5.1],
            "T3": [5.5, 5.0],
        }
    )
    pair_table = pd.DataFrame(
        {
            "pair_id": ["P1", "P2", "P3"],
            "normal": ["N1", "N2", "N3"],
            "tumor": ["T1", "T2", "T3"],
        }
    )

    result = paired_ttest_by_probe(
        expression,
        pair_table,
        tumor_column="tumor",
        normal_column="normal",
    )

    assert result.columns.tolist() == [
        "probe_id",
        "mean_tumor",
        "mean_paired_normal",
        "mean_paired_difference",
        "p_value",
        "adjusted_p_value",
    ]
    probe_1 = result.set_index("probe_id").loc["probe_1"]
    assert probe_1["mean_paired_difference"] > 0
