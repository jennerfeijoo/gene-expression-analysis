"""Small statistical helpers for exploratory expression analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind


def benjamini_hochberg(p_values: pd.Series | np.ndarray) -> np.ndarray:
    """Adjust p-values with the Benjamini-Hochberg procedure."""
    values = np.asarray(p_values, dtype=float)
    adjusted = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values)

    if not valid.any():
        return adjusted

    valid_values = values[valid]
    order = np.argsort(valid_values)
    ranked = valid_values[order]
    ranks = np.arange(1, len(ranked) + 1)
    ranked_adjusted = ranked * len(ranked) / ranks
    ranked_adjusted = np.minimum.accumulate(ranked_adjusted[::-1])[::-1]
    ranked_adjusted = np.clip(ranked_adjusted, 0.0, 1.0)

    restored = np.empty_like(ranked_adjusted)
    restored[order] = ranked_adjusted
    adjusted[valid] = restored
    return adjusted


def welch_ttest_by_probe(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    group_a: str,
    group_b: str,
) -> pd.DataFrame:
    """Compare two sample groups probe by probe with Welch's t-test."""
    required_metadata = {"sample_accession", group_column}
    missing_metadata = required_metadata - set(metadata.columns)
    if missing_metadata:
        raise ValueError(f"Metadata is missing columns: {sorted(missing_metadata)}")
    if "ID_REF" not in expression.columns:
        raise ValueError("Expression table must contain an ID_REF column.")

    samples_a = metadata.loc[
        metadata[group_column] == group_a,
        "sample_accession",
    ].tolist()
    samples_b = metadata.loc[
        metadata[group_column] == group_b,
        "sample_accession",
    ].tolist()
    if not samples_a or not samples_b:
        raise ValueError("Both comparison groups must contain samples.")

    missing_samples = (set(samples_a) | set(samples_b)) - set(expression.columns)
    if missing_samples:
        raise ValueError(
            f"Expression table is missing {len(missing_samples)} comparison samples."
        )

    values_a = expression[samples_a].to_numpy(dtype=float)
    values_b = expression[samples_b].to_numpy(dtype=float)
    test_result = ttest_ind(
        values_a,
        values_b,
        axis=1,
        equal_var=False,
        nan_policy="omit",
    )

    return pd.DataFrame(
        {
            "probe_id": expression["ID_REF"].astype(str),
            f"mean_{group_a}": np.nanmean(values_a, axis=1),
            f"mean_{group_b}": np.nanmean(values_b, axis=1),
            "log2_fold_change": (
                np.nanmean(values_a, axis=1) - np.nanmean(values_b, axis=1)
            ),
            "p_value": test_result.pvalue,
        }
    )


def rank_group_comparison(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
    group_column: str,
    group_a: str,
    group_b: str,
) -> pd.DataFrame:
    """Run a Welch comparison, adjust p-values, and rank the probes."""
    comparison = welch_ttest_by_probe(
        expression,
        metadata,
        group_column,
        group_a,
        group_b,
    )
    comparison["adjusted_p_value"] = benjamini_hochberg(comparison["p_value"])
    comparison["_absolute_log2_fold_change"] = comparison[
        "log2_fold_change"
    ].abs()
    return (
        comparison.sort_values(
            ["adjusted_p_value", "_absolute_log2_fold_change"],
            ascending=[True, False],
            na_position="last",
        )
        .drop(columns="_absolute_log2_fold_change")
        .reset_index(drop=True)
    )
