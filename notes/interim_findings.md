# Interim findings

This note records the current state of the GSE44076 analysis before any
biological interpretation or pathway analysis.

## Dataset and sample groups

The project uses the processed series matrix for GEO accession GSE44076 on
platform GPL13667 (Affymetrix Human Genome U219 Array). The metadata support
three groups:

- Healthy mucosa: 50 samples
- Paired normal mucosa: 98 samples
- Tumor: 98 samples

The tumor and paired-normal samples form 98 complete patient pairs.

## Work completed

The ten notebooks currently cover:

1. loading the GEO series matrix;
2. organizing sample metadata and group labels;
3. checking expression distributions;
4. reviewing preprocessing and sample-level summaries;
5. exploratory PCA;
6. an initial unpaired group comparison;
7. recovery and use of the paired design;
8. GPL13667 probe annotation;
9. a first gene-level summary; and
10. a sensitivity check of the probe-to-gene collapse.

The expression matrix contains 49,386 probes across 246 aligned samples.
Its range and distributions appear consistent with values that have already
been transformed to a log2-like scale. This is a descriptive assessment of
the supplied matrix, not a reconstruction of the original preprocessing.
Sample medians and IQRs were reviewed, and flagged samples were retained.

## Exploratory structure and comparison

PCA showed group-associated structure, with overlap between groups. It is
useful as a view of broad sample-level variation, but it does not establish
biological separation or explain the causes of that variation.

The paired tumor versus paired-normal comparison is preferred over the
earlier unpaired analysis because it uses the matched samples from the same
98 individuals. This accounts for pairing in the test and avoids treating
matched observations as independent. The unpaired ranking remains only an
exploratory reference.

## Annotation and gene-level summary

GPL13667 annotation mapped 48,784 of 49,386 ranked probes to gene symbols;
602 probes had no gene symbol. Mapping is incomplete and can be many-to-one,
and some probe records contain multiple gene symbols.

The first gene-level table keeps, for each gene symbol, the probe with the
smallest probe-level adjusted p-value. It also records how many probes mapped
to that gene. Adjusted p-values were not recalculated after this collapse.

The sensitivity check compared that rule with selecting the probe having the
largest absolute paired mean difference. Only 13 of the top 50 genes
overlapped. Probe directions also disagreed for some multi-probe genes. The
top gene-level ordering therefore depends materially on annotation and
collapse choices.

## What can be said so far

- The local files support a reproducible 98-pair probe-level comparison.
- Broad expression structure is associated with the recorded sample groups,
  while substantial overlap remains.
- Probe annotation enables exploratory gene-symbol summaries.
- Gene-level rankings are sensitive to the representative-probe rule.

## What should not be claimed yet

- The ranked probes or genes are not established biomarkers.
- The results do not demonstrate mechanism, causality, diagnostic value, or
  clinical utility.
- A single representative probe does not provide a definitive gene-level
  differential-expression result.
- Pathway interpretation is premature while annotation and collapse choices
  remain under review.

## Next steps

- Review multi-symbol probe mappings and genes with conflicting probe
  directions.
- Decide whether a different gene-level aggregation method is justified.
- Reconsider multiple-testing handling at the gene-summary level.
- Document any remaining platform and preprocessing assumptions.
- Begin pathway analysis only after the gene-level input definition is
  methodologically defensible.
