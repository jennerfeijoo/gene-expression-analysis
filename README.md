# Gene Expression Analysis

This repository is being developed as a reproducible gene-expression analysis
portfolio project. The initial workflow uses GEO accession GSE44076.

The first loading notebook,
[`notebooks/01_load_series_matrix.ipynb`](notebooks/01_load_series_matrix.ipynb),
checks the local series matrix structure, previews available sample metadata,
and loads the expression table without performing downstream biological
analysis.

The second notebook,
[`notebooks/02_sample_metadata.ipynb`](notebooks/02_sample_metadata.ipynb),
organizes sample metadata and documents source-based group labels.

The third notebook,
[`notebooks/03_first_expression_plots.ipynb`](notebooks/03_first_expression_plots.ipynb),
checks sample alignment and creates initial expression-distribution plots.

The fourth notebook,
[`notebooks/04_check_expression_preprocessing.ipynb`](notebooks/04_check_expression_preprocessing.ipynb),
checks the supplied value scale and sample-level distribution summaries before
PCA or differential analysis.

The fifth notebook,
[`notebooks/05_first_pca.ipynb`](notebooks/05_first_pca.ipynb), adds a first
exploratory PCA of sample-level structure.

The sixth notebook,
[`notebooks/06_first_group_comparison.ipynb`](notebooks/06_first_group_comparison.ipynb),
adds an exploratory tumor versus paired-normal probe comparison.

The seventh notebook,
[`notebooks/07_check_paired_design.ipynb`](notebooks/07_check_paired_design.ipynb),
checks the matched design and runs a paired comparison when pairing is valid.
The earlier unpaired ranking remains exploratory.

The eighth notebook,
[`notebooks/08_probe_annotation.ipynb`](notebooks/08_probe_annotation.ipynb),
checks GPL13667 probe-to-gene annotation before interpretation.

Raw GEO files and the generated `data/processed/sample_metadata.csv` table are
kept locally and are not tracked by Git. Generated figures under
`reports/figures/` are tracked.
