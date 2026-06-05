# GEO metadata check

Dataset checked: GSE44076  
Source: NCBI Gene Expression Omnibus

## What the GEO page says

GSE44076 contains gene expression data from healthy, adjacent normal, and tumor colon cells.

Basic information:

- Organism: Homo sapiens
- Experiment type: expression profiling by array
- Platform: GPL13667, Affymetrix Human Genome U219 Array
- Number of samples: 246
- Public since: March 14, 2014
- Last GEO update checked: August 31, 2023

## Study design notes

The dataset includes:

- healthy colon mucosa samples from donors without colonic lesions
- paired adjacent normal mucosa samples from patients
- paired tumor samples from patients with colorectal adenocarcinoma

The GEO summary describes 98 individuals with paired adjacent normal mucosa and tumor samples, plus 50 healthy colon mucosa samples.

## Why this dataset is usable for this project

This looks suitable for a first gene expression analysis because:

- it has clear biological groups
- it has a manageable number of samples
- it uses human colon tissue
- it allows exploratory comparison of healthy, adjacent normal, and tumor samples
- it gives room to practice metadata handling before doing any biomarker-style analysis

## Things to be careful about

- Tumor vs adjacent normal is not the same comparison as tumor vs healthy donor mucosa.
- Adjacent normal tissue may already reflect field effects or patient-specific context.
- Samples are not all independent if paired tumor/normal samples from the same patients are used.
- The platform is microarray, not RNA-seq.
- Any candidate genes from this project would be exploratory only.

## Current decision

Proceed with GSE44076 for the first pass.

The next step is to find the easiest processed data file and inspect the sample metadata before doing any analysis.