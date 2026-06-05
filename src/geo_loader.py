"""Utilities for reading GEO series matrix files."""

from __future__ import annotations

import csv
import gzip
from io import StringIO
from pathlib import Path
import re

import pandas as pd

TABLE_BEGIN = "!series_matrix_table_begin"
TABLE_END = "!series_matrix_table_end"


def read_geo_series_lines(path: str | Path) -> list[str]:
    """Read a plain-text or gzipped GEO series matrix into a list of lines."""
    file_path = Path(path)
    if not file_path.is_file():
        raise FileNotFoundError(f"GEO series matrix not found: {file_path}")

    open_file = gzip.open if file_path.suffix == ".gz" else open
    with open_file(file_path, mode="rt", encoding="utf-8", errors="replace") as handle:
        return handle.readlines()


def extract_geo_metadata(lines: list[str]) -> dict[str, list[str]]:
    """Extract tab-separated values from GEO metadata lines beginning with ``!``."""
    metadata: dict[str, list[str]] = {}

    for line in lines:
        if line.strip() == TABLE_BEGIN:
            break
        if not line.startswith("!"):
            continue

        fields = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"))
        key = fields[0].removeprefix("!")
        metadata.setdefault(key, []).extend(fields[1:])

    return metadata


def extract_characteristics(metadata: dict[str, list[str]]) -> pd.DataFrame:
    """Return sample characteristics as aligned, named columns."""
    sample_count = len(metadata.get("Sample_geo_accession", []))
    values = metadata.get("Sample_characteristics_ch1", [])

    if not values:
        return pd.DataFrame(index=range(sample_count))
    if sample_count == 0 or len(values) % sample_count != 0:
        raise ValueError("Sample characteristics are not aligned with sample accessions.")

    characteristics: dict[str, list[str]] = {}
    for start in range(0, len(values), sample_count):
        block = values[start : start + sample_count]
        parsed = [value.split(":", maxsplit=1) for value in block]
        labels = {parts[0].strip() for parts in parsed if len(parts) == 2}

        if len(labels) == 1:
            label = next(iter(labels))
            column = re.sub(r"\W+", "_", label.strip().lower()).strip("_")
            block_values = [parts[1].strip() for parts in parsed]
        else:
            column = f"characteristic_{start // sample_count + 1}"
            block_values = block

        base_column = column or f"characteristic_{start // sample_count + 1}"
        suffix = 2
        while column in characteristics:
            column = f"{base_column}_{suffix}"
            suffix += 1
        characteristics[column] = block_values

    return pd.DataFrame(characteristics)


def build_sample_metadata_table(lines: list[str]) -> pd.DataFrame:
    """Build one row per sample from GEO series metadata lines."""
    metadata = extract_geo_metadata(lines)
    accessions = metadata.get("Sample_geo_accession", [])
    if not accessions:
        raise ValueError("No sample accessions were found in the GEO metadata.")

    sample_count = len(accessions)

    def sample_values(key: str) -> list[str | None]:
        values = metadata.get(key, [])
        if not values:
            return [None] * sample_count
        if len(values) != sample_count:
            raise ValueError(f"{key} is not aligned with sample accessions.")
        return values

    sample_metadata = pd.DataFrame(
        {
            "sample_accession": accessions,
            "title": sample_values("Sample_title"),
            "source_name": sample_values("Sample_source_name_ch1"),
        }
    )
    return pd.concat(
        [sample_metadata, extract_characteristics(metadata)],
        axis="columns",
    )


def align_expression_with_metadata(
    expression: pd.DataFrame,
    metadata: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Align expression sample columns to metadata row order."""
    if "sample_accession" not in metadata.columns:
        raise ValueError("Metadata must contain a sample_accession column.")
    if metadata["sample_accession"].duplicated().any():
        raise ValueError("Metadata sample accessions must be unique.")

    accessions = metadata["sample_accession"].tolist()
    identifier_columns = ["ID_REF"] if "ID_REF" in expression.columns else []
    expression_samples = [
        column for column in expression.columns if column not in identifier_columns
    ]

    missing_expression = set(accessions) - set(expression_samples)
    missing_metadata = set(expression_samples) - set(accessions)
    if missing_expression or missing_metadata:
        raise ValueError(
            "Expression and metadata sample accessions do not match: "
            f"{len(missing_expression)} missing from expression, "
            f"{len(missing_metadata)} missing from metadata."
        )

    aligned_expression = expression.loc[:, identifier_columns + accessions].copy()
    aligned_metadata = metadata.set_index("sample_accession").loc[accessions].reset_index()
    return aligned_expression, aligned_metadata


def convert_expression_to_numeric(expression: pd.DataFrame) -> pd.DataFrame:
    """Convert expression sample columns to numeric values without mutating input."""
    converted = expression.copy()
    sample_columns = [
        column for column in converted.columns if column != "ID_REF"
    ]
    converted[sample_columns] = converted[sample_columns].apply(
        pd.to_numeric,
        errors="coerce",
    )
    return converted


def summarize_expression_values(expression: pd.DataFrame) -> pd.Series:
    """Summarize numeric expression values across all samples."""
    sample_columns = [
        column for column in expression.columns if column != "ID_REF"
    ]
    values = expression[sample_columns].stack()
    quantiles = values.quantile([0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99])

    return pd.Series(
        {
            "count": values.count(),
            "missing": expression[sample_columns].isna().sum().sum(),
            "minimum": values.min(),
            "q01": quantiles.loc[0.01],
            "q05": quantiles.loc[0.05],
            "q25": quantiles.loc[0.25],
            "median": quantiles.loc[0.5],
            "q75": quantiles.loc[0.75],
            "q95": quantiles.loc[0.95],
            "q99": quantiles.loc[0.99],
            "maximum": values.max(),
        },
        name="expression_value",
    )


def summarize_sample_distributions(expression: pd.DataFrame) -> pd.DataFrame:
    """Summarize the expression distribution within each sample column."""
    sample_columns = [
        column for column in expression.columns if column != "ID_REF"
    ]
    values = expression[sample_columns]
    summary = pd.DataFrame(
        {
            "minimum": values.min(),
            "q25": values.quantile(0.25),
            "median": values.median(),
            "q75": values.quantile(0.75),
            "maximum": values.max(),
            "missing": values.isna().sum(),
        }
    )
    summary["iqr"] = summary["q75"] - summary["q25"]
    summary.index.name = "sample_accession"
    return summary


def load_geo_expression_table(path: str | Path) -> pd.DataFrame:
    """Load the expression table delimited by GEO series matrix markers."""
    lines = read_geo_series_lines(path)

    try:
        start = next(
            index for index, line in enumerate(lines) if line.strip() == TABLE_BEGIN
        )
        end = next(
            index
            for index, line in enumerate(lines[start + 1 :], start=start + 1)
            if line.strip() == TABLE_END
        )
    except StopIteration as exc:
        raise ValueError("GEO expression table markers are missing or incomplete.") from exc

    table_lines = lines[start + 1 : end]
    if not table_lines:
        raise ValueError("GEO expression table is empty.")

    return pd.read_csv(StringIO("".join(table_lines)), sep="\t", quotechar='"')
