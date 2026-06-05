"""Utilities for reading GEO series matrix files."""

from __future__ import annotations

import csv
import gzip
from io import StringIO
from pathlib import Path

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
