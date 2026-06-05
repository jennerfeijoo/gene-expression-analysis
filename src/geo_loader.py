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
