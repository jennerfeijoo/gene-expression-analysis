"""Helpers for loading and joining GEO platform annotations."""

from __future__ import annotations

import gzip
from pathlib import Path
import re

import pandas as pd

PLATFORM_TABLE_BEGIN = "!platform_table_begin"
PLATFORM_TABLE_END = "!platform_table_end"

ANNOTATION_ALIASES = {
    "probe_id": ("probe_id", "id", "id_ref", "probe set id", "probeset id"),
    "gene_symbol": ("gene_symbol", "gene symbol", "symbol"),
    "gene_title": (
        "gene_title",
        "gene title",
        "gene_name",
        "gene name",
        "description",
    ),
    "entrez_id": (
        "entrez_id",
        "entrez gene",
        "entrez gene id",
        "entrezgene",
    ),
}


def _normalize_column_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", name.lower()).strip()


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open(mode="rt", encoding="utf-8", errors="replace")


def load_gpl_annotation(path: str | Path) -> pd.DataFrame:
    """Load the annotation table from a GEO platform text or gzip file."""
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"GPL annotation file not found: {file_path}")

    header_line = None
    row_count = None
    with _open_text(file_path) as handle:
        for line_number, line in enumerate(handle):
            if line.strip() == PLATFORM_TABLE_BEGIN:
                header_line = line_number + 1
                continue
            if line.strip() == PLATFORM_TABLE_END and header_line is not None:
                row_count = line_number - header_line - 1
                break

    if header_line is None:
        return pd.read_csv(file_path, sep="\t", compression="infer", dtype=str)

    if row_count is None:
        raise ValueError("GPL platform annotation table is incomplete.")

    return pd.read_csv(
        file_path,
        sep="\t",
        compression="infer",
        skiprows=header_line,
        nrows=row_count,
        dtype=str,
        low_memory=False,
    )


def standardize_annotation_columns(annotation: pd.DataFrame) -> pd.DataFrame:
    """Standardize common GPL annotation columns and missing-value markers."""
    normalized_columns = {
        _normalize_column_name(column): column for column in annotation.columns
    }
    standardized = pd.DataFrame(index=annotation.index)

    for target, aliases in ANNOTATION_ALIASES.items():
        source = next(
            (
                normalized_columns[_normalize_column_name(alias)]
                for alias in aliases
                if _normalize_column_name(alias) in normalized_columns
            ),
            None,
        )
        standardized[target] = annotation[source] if source else pd.NA

    if standardized["probe_id"].isna().all():
        raise ValueError("No probe identifier column was found in the annotation.")

    return standardized.replace(
        {
            "": pd.NA,
            "---": pd.NA,
            "NA": pd.NA,
            "nan": pd.NA,
        }
    )


def annotate_probe_ranking(
    ranking: pd.DataFrame,
    annotation: pd.DataFrame,
) -> pd.DataFrame:
    """Left-join standardized annotation onto a probe-level ranking."""
    if "probe_id" not in ranking.columns:
        raise ValueError("Ranking must contain a probe_id column.")

    standardized = standardize_annotation_columns(annotation)
    standardized = standardized.drop_duplicates("probe_id", keep="first")
    return ranking.merge(
        standardized,
        on="probe_id",
        how="left",
        validate="many_to_one",
    )


def _expand_gene_symbols(annotated_ranking: pd.DataFrame) -> pd.DataFrame:
    if "gene_symbol" not in annotated_ranking.columns:
        raise ValueError("Annotated ranking must contain a gene_symbol column.")

    expanded = annotated_ranking.dropna(subset=["gene_symbol"]).copy()
    expanded["gene_symbol"] = expanded["gene_symbol"].str.split(r"\s*///\s*")
    expanded = expanded.explode("gene_symbol")
    expanded["gene_symbol"] = expanded["gene_symbol"].str.strip()
    return expanded[expanded["gene_symbol"].ne("")]


def summarize_probes_per_gene(
    annotated_ranking: pd.DataFrame,
) -> pd.DataFrame:
    """Count distinct annotated probes assigned to each gene symbol."""
    expanded = _expand_gene_symbols(annotated_ranking)
    return (
        expanded.groupby("gene_symbol", as_index=False)["probe_id"]
        .nunique()
        .rename(columns={"probe_id": "probe_count"})
        .sort_values(["probe_count", "gene_symbol"], ascending=[False, True])
        .reset_index(drop=True)
    )


def collapse_to_gene_level(
    annotated_ranking: pd.DataFrame,
) -> pd.DataFrame:
    """Keep the smallest-adjusted-p probe as a representative for each gene."""
    required_columns = {
        "probe_id",
        "gene_symbol",
        "adjusted_p_value",
        "mean_paired_difference",
    }
    missing_columns = required_columns - set(annotated_ranking.columns)
    if missing_columns:
        raise ValueError(
            f"Annotated ranking is missing columns: {sorted(missing_columns)}"
        )

    expanded = _expand_gene_symbols(annotated_ranking)
    probe_counts = summarize_probes_per_gene(annotated_ranking)
    expanded["_absolute_paired_difference"] = expanded[
        "mean_paired_difference"
    ].abs()
    representatives = (
        expanded.sort_values(
            [
                "gene_symbol",
                "adjusted_p_value",
                "_absolute_paired_difference",
                "probe_id",
            ],
            ascending=[True, True, False, True],
            na_position="last",
        )
        .drop_duplicates("gene_symbol", keep="first")
        .drop(columns="_absolute_paired_difference")
    )
    return (
        representatives.merge(
            probe_counts,
            on="gene_symbol",
            how="left",
            validate="one_to_one",
        )
        .sort_values(
            ["adjusted_p_value", "gene_symbol"],
            ascending=[True, True],
            na_position="last",
        )
        .reset_index(drop=True)
    )


def summarize_gene_probe_consistency(
    annotated_ranking: pd.DataFrame,
) -> pd.DataFrame:
    """Summarize effect-direction agreement among probes for each gene."""
    if "mean_paired_difference" not in annotated_ranking.columns:
        raise ValueError(
            "Annotated ranking must contain a mean_paired_difference column."
        )

    expanded = _expand_gene_symbols(annotated_ranking).drop_duplicates(
        ["gene_symbol", "probe_id"]
    )
    expanded["_direction"] = 0
    expanded.loc[expanded["mean_paired_difference"] > 0, "_direction"] = 1
    expanded.loc[expanded["mean_paired_difference"] < 0, "_direction"] = -1

    summary = (
        expanded.groupby("gene_symbol", as_index=False)
        .agg(
            probe_count=("probe_id", "nunique"),
            positive_probe_count=("_direction", lambda values: (values > 0).sum()),
            negative_probe_count=("_direction", lambda values: (values < 0).sum()),
            zero_probe_count=("_direction", lambda values: (values == 0).sum()),
        )
    )
    summary["direction_conflict"] = (
        summary["positive_probe_count"].gt(0)
        & summary["negative_probe_count"].gt(0)
    )
    summary["same_direction"] = ~summary["direction_conflict"]
    return summary.sort_values(
        ["probe_count", "gene_symbol"],
        ascending=[False, True],
    ).reset_index(drop=True)


def collapse_by_largest_abs_effect(
    annotated_ranking: pd.DataFrame,
) -> pd.DataFrame:
    """Keep the largest-absolute-effect probe as a representative per gene."""
    required_columns = {
        "probe_id",
        "gene_symbol",
        "adjusted_p_value",
        "mean_paired_difference",
    }
    missing_columns = required_columns - set(annotated_ranking.columns)
    if missing_columns:
        raise ValueError(
            f"Annotated ranking is missing columns: {sorted(missing_columns)}"
        )

    expanded = _expand_gene_symbols(annotated_ranking)
    probe_counts = summarize_probes_per_gene(annotated_ranking)
    expanded["_absolute_paired_difference"] = expanded[
        "mean_paired_difference"
    ].abs()
    representatives = (
        expanded.sort_values(
            [
                "gene_symbol",
                "_absolute_paired_difference",
                "adjusted_p_value",
                "probe_id",
            ],
            ascending=[True, False, True, True],
            na_position="last",
        )
        .drop_duplicates("gene_symbol", keep="first")
        .drop(columns="_absolute_paired_difference")
    )
    result = representatives.merge(
        probe_counts,
        on="gene_symbol",
        how="left",
        validate="one_to_one",
    )
    result["_absolute_paired_difference"] = result["mean_paired_difference"].abs()
    return (
        result.sort_values(
            ["_absolute_paired_difference", "adjusted_p_value", "gene_symbol"],
            ascending=[False, True, True],
            na_position="last",
        )
        .drop(columns="_absolute_paired_difference")
        .reset_index(drop=True)
    )


def compare_gene_collapse_rules(
    table_a: pd.DataFrame,
    table_b: pd.DataFrame,
    top_n: int = 50,
) -> pd.DataFrame:
    """Compare membership and rank among the top genes from two rules."""
    if top_n <= 0:
        raise ValueError("top_n must be positive.")
    if "gene_symbol" not in table_a.columns or "gene_symbol" not in table_b.columns:
        raise ValueError("Both tables must contain a gene_symbol column.")

    top_a = table_a["gene_symbol"].drop_duplicates().head(top_n).tolist()
    top_b = table_b["gene_symbol"].drop_duplicates().head(top_n).tolist()
    symbols = list(dict.fromkeys(top_a + top_b))
    rank_a = {symbol: rank for rank, symbol in enumerate(top_a, start=1)}
    rank_b = {symbol: rank for rank, symbol in enumerate(top_b, start=1)}

    comparison = pd.DataFrame({"gene_symbol": symbols})
    comparison["rank_rule_a"] = comparison["gene_symbol"].map(rank_a).astype("Int64")
    comparison["rank_rule_b"] = comparison["gene_symbol"].map(rank_b).astype("Int64")
    comparison["in_rule_a"] = comparison["rank_rule_a"].notna()
    comparison["in_rule_b"] = comparison["rank_rule_b"].notna()
    comparison["in_both"] = comparison["in_rule_a"] & comparison["in_rule_b"]
    return comparison
