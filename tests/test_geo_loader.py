"""Tests for the GEO series matrix loader."""

from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd
import pytest

from src.geo_loader import (
    build_sample_metadata_table,
    extract_characteristics,
    extract_geo_metadata,
    load_geo_expression_table,
    read_geo_series_lines,
)

GEO_EXAMPLE = """!Series_geo_accession\t"GSE00000"
!Sample_geo_accession\t"GSM1"\t"GSM2"
!Sample_title\t"control"\t"example"
!Sample_source_name_ch1\t"healthy tissue"\t"example tissue"
!Sample_characteristics_ch1\t"group: control"\t"group: example"
!Sample_characteristics_ch1\t"batch: one"\t"batch: two"
!series_matrix_table_begin
"ID_REF"\t"GSM1"\t"GSM2"
"probe_1"\t1.0\t2.0
"probe_2"\t3.0\t4.0
!series_matrix_table_end
"""


def write_gzipped_example(path: Path, content: str = GEO_EXAMPLE) -> None:
    """Write a small gzipped GEO-like fixture."""
    with gzip.open(path, mode="wt", encoding="utf-8") as handle:
        handle.write(content)


def test_read_geo_series_lines_reads_gzip(tmp_path: Path) -> None:
    path = tmp_path / "series_matrix.txt.gz"
    write_gzipped_example(path)

    lines = read_geo_series_lines(path)

    assert lines[0].startswith("!Series_geo_accession")
    assert lines[-1].strip() == "!series_matrix_table_end"


def test_extract_geo_metadata_collects_sample_fields() -> None:
    metadata = extract_geo_metadata(GEO_EXAMPLE.splitlines(keepends=True))

    assert metadata["Series_geo_accession"] == ["GSE00000"]
    assert metadata["Sample_geo_accession"] == ["GSM1", "GSM2"]
    assert metadata["Sample_title"] == ["control", "example"]
    assert "series_matrix_table_begin" not in metadata


def test_build_sample_metadata_table_aligns_sample_fields() -> None:
    sample_metadata = build_sample_metadata_table(
        GEO_EXAMPLE.splitlines(keepends=True)
    )

    expected = pd.DataFrame(
        {
            "sample_accession": ["GSM1", "GSM2"],
            "title": ["control", "example"],
            "source_name": ["healthy tissue", "example tissue"],
            "group": ["control", "example"],
            "batch": ["one", "two"],
        }
    )
    pd.testing.assert_frame_equal(sample_metadata, expected)


def test_extract_characteristics_preserves_repeated_fields() -> None:
    metadata = {
        "Sample_geo_accession": ["GSM1", "GSM2"],
        "Sample_characteristics_ch1": [
            "status: control",
            "status: example",
            "status: untreated",
            "status: treated",
        ],
    }

    characteristics = extract_characteristics(metadata)

    assert characteristics.columns.tolist() == ["status", "status_2"]
    assert characteristics["status"].tolist() == ["control", "example"]
    assert characteristics["status_2"].tolist() == ["untreated", "treated"]


def test_load_geo_expression_table(tmp_path: Path) -> None:
    path = tmp_path / "series_matrix.txt.gz"
    write_gzipped_example(path)

    expression = load_geo_expression_table(path)

    expected = pd.DataFrame(
        {
            "ID_REF": ["probe_1", "probe_2"],
            "GSM1": [1.0, 3.0],
            "GSM2": [2.0, 4.0],
        }
    )
    pd.testing.assert_frame_equal(expression, expected)


def test_load_geo_expression_table_rejects_missing_markers(tmp_path: Path) -> None:
    path = tmp_path / "invalid.txt.gz"
    write_gzipped_example(path, "!Series_geo_accession\t\"GSE00000\"\n")

    with pytest.raises(ValueError, match="markers"):
        load_geo_expression_table(path)
