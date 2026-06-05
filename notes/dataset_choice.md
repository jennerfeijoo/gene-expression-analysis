# Dataset choice

For the first version of this project, I will start with a public GEO gene expression dataset instead of TCGA.

The reason is practical: I want to focus first on the basic workflow of expression analysis before working with larger cancer genomics resources.

## Selected starting dataset

Candidate dataset: GSE44076

Initial reason for choosing it:

- It is related to colon cancer gene expression.
- It has been used in public bioinformatics analyses involving colon cancer and healthy/normal comparisons.
- It should be more manageable for a first expression-analysis repository than TCGA.
- It fits the project goal: practice metadata review, expression matrix handling, exploratory analysis, group comparison, and cautious biomarker-style ranking.

## What I need to check next

Before using the dataset, I still need to verify:

- sample groups
- number of samples
- platform used
- whether processed expression data are available
- whether gene identifiers are already mapped
- whether the metadata are clean enough for a first analysis
- whether the dataset supports a simple comparison such as tumor vs normal

## Current decision

Use GSE44076 as the first dataset to investigate.

If the metadata or processed files are too complicated, I will switch to another small GEO dataset rather than forcing the analysis.