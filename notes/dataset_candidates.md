# Dataset candidates

I am using this note to choose a small public gene expression dataset for the first version of the project.

The goal is not to find a perfect dataset. The goal is to practice a clean workflow:
metadata review, expression matrix loading, exploratory plots, simple group comparison, and cautious interpretation.

## Candidate source

### NCBI GEO

GEO is a public repository for functional genomics data, including gene expression datasets from microarray and sequencing experiments.

For this project, I want a dataset that has:

- clear sample groups
- manageable sample size
- human disease relevance
- downloadable expression matrix or supplementary table
- usable metadata
- a publication or clear GEO description

## Candidate datasets to check

### Option 1: breast cancer gene expression dataset

Reason to check:
- familiar biomedical context
- useful for practicing tumor vs normal or subtype comparisons
- relevant to biomarker-style analysis

Possible issue:
- datasets can be large or have complicated metadata

### Option 2: colon cancer gene expression dataset

Reason to check:
- common benchmark area in gene expression analysis
- often has clearer case/control comparisons

Possible issue:
- older microarray platforms may require extra annotation work

### Option 3: small inflammation or immune-response expression dataset

Reason to check:
- may be easier to interpret biologically
- useful for pathway-style thinking later

Possible issue:
- may be less directly connected to personalized medicine than cancer datasets

## Current decision

I will first look for a small cancer-related GEO dataset with clear groups and accessible processed expression data.

I will avoid starting with TCGA in this repository because it adds more complexity than I need for the first gene expression project.