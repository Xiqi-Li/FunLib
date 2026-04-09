# FunLib

`FunLib` is an R package that collects common Woodman Lab analysis utilities for genomic, mutation, expression, clustering, pathway, survival, immune, and spatial workflows.

## Installation

Install from the local repository root:

```r
install.packages(".", repos = NULL, type = "source")
```

Or with `devtools`:

```r
devtools::install(".")
```

## Contents

- `00_LIBRARY_INDEX.R`: master index, provenance notes, duplicate-function review, and usage guidance
- `R/01_genomic_data.R`: genomic ranges, coordinates, and copy-number utilities
- `R/02_mutation_analysis.R`: mutation filtering, mutation summaries, and oncoplot helpers
- `R/03_expression_analysis.R`: differential expression and related expression-analysis utilities
- `R/04_clustering.R`: clustering workflows and helper functions
- `R/05_gsea_pathway.R`: gene set enrichment and pathway analysis helpers
- `R/06_survival_analysis.R`: Kaplan-Meier and Cox survival utilities
- `R/07_immune_deconvolution.R`: immune deconvolution wrappers
- `R/08_spatial_immune_profiling.R`: spatial / mIF immune profiling functions
- `R/09_oncokb_annotation.R`: OncoKB annotation helpers
- `R/10_statistics_utilities.R`: general statistics, data handling, and survival utilities
- `R/11_rtilandscape_utils.R`: RTI landscape, fusion, and annotation helpers

## Usage

Load the package and call exported functions directly:

```r
library(FunLib)

get_oncokb_version(token = Sys.getenv("ONCOKB_API_TOKEN"))
```

If you want to inspect provenance and module contents, use `00_LIBRARY_INDEX.R`.

## Notes

- This is now a package-style repository; the functional code lives under `R/`.
- Several functions depend on external packages that are not bundled in this repository.
- Some functions expect access to external systems or tokens, for example OncoKB (`ONCOKB_API_TOKEN`) or Foundry credentials. Tokens are not stored in this repository and should be provided through the environment or function arguments.
- This repository contains code only. Input datasets, credentials, and generated analysis outputs should remain outside version control.
