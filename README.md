# Gene Expression Analysis

This repository is being developed as a reproducible gene-expression analysis
portfolio project. The initial workflow uses GEO accession GSE44076.

The first loading notebook,
[`notebooks/01_load_series_matrix.ipynb`](notebooks/01_load_series_matrix.ipynb),
checks the local series matrix structure, previews available sample metadata,
and loads the expression table without performing downstream biological
analysis.

Raw GEO files are stored under `data/raw/` for local use and are not tracked by
Git.
