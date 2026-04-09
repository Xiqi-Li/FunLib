# WoodmanLab Function Library

Curated R function library for common Woodman Lab analysis workflows. The repository currently exposes a set of numbered script modules that can be sourced independently or combined in an analysis session.

## Contents

- `00_LIBRARY_INDEX.R`: master index, provenance notes, duplicate-function review, and usage guidance
- `01_genomic_data.R`: genomic ranges, coordinates, and copy-number utilities
- `02_mutation_analysis.R`: mutation filtering, mutation summaries, and oncoplot helpers
- `03_expression_analysis.R`: differential expression and related expression-analysis utilities
- `04_clustering.R`: clustering workflows and helper functions
- `05_gsea_pathway.R`: gene set enrichment and pathway analysis helpers
- `06_survival_analysis.R`: Kaplan-Meier and Cox survival utilities
- `07_immune_deconvolution.R`: immune deconvolution wrappers
- `08_spatial_immune_profiling.R`: spatial / mIF immune profiling functions
- `09_oncokb_annotation.R`: OncoKB annotation helpers
- `10_statistics_utilities.R`: general statistics, data handling, and survival utilities
- `11_rtilandscape_utils.R`: RTI landscape, fusion, and annotation helpers

## Usage

Source only the modules you need:

```r
source("01_genomic_data.R")
source("02_mutation_analysis.R")
source("09_oncokb_annotation.R")
```

Use `00_LIBRARY_INDEX.R` to identify the recommended implementation when similar functions exist in more than one historical source.

## Notes

- Several functions depend on external packages that are not bundled in this repository.
- Some functions expect access to external systems or tokens, for example OncoKB (`ONCOKB_API_TOKEN`) or Foundry credentials. Tokens are not stored in this repository and should be provided through the environment or function arguments.
- This repository contains code only. Input datasets, credentials, and generated analysis outputs should remain outside version control.
