"""Tests for exploratory group-comparison helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.analysis import benjamini_hochberg, rank_group_comparison


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
