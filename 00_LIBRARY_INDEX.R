# =============================================================================
# WOODMANLAB R FUNCTION LIBRARY — MASTER INDEX
# Generated: 2026-04-08
# =============================================================================
#
# This directory contains 11 categorized function library files consolidating
# ~330 R/Rmd files found across the WoodmanLab project tree.
#
# HOW TO USE:
#   library(FunLib)
#   # Functions are loaded from the package namespace.
#   # Raw module files live under the package R/ directory.
#
# =============================================================================
# FILE INVENTORY
# =============================================================================
#
#  File                           | Functions | Source lines | Description
#  -------------------------------|-----------|--------------|-----------------------------
#  01_genomic_data.R              |    12     |    657       | GRanges, CNV plots, coordinate utils
#  02_mutation_analysis.R         |    13     |   1919       | SNV/indel filtering, oncoplots
#  03_expression_analysis.R       |    10     |   1279       | DGE, unsupervised, deconvolution
#  04_clustering.R                |    20     |    616       | 13-algorithm clustering, consensus
#  05_gsea_pathway.R              |     7     |    344       | fGSEA, pathway plots, GO enrichment
#  06_survival_analysis.R         |     4     |    266       | KM plots, Cox regression
#  07_immune_deconvolution.R      |     5     |    399       | CIBERSORT, DeconRNASeq, wrapper
#  08_spatial_immune_profiling.R  |    35     |   2170       | mIF/Vectra stats + visualization
#  09_oncokb_annotation.R         |    13     |    513       | OncoKB API queries + MAF annotation
#  10_statistics_utilities.R      |    19     |    773       | Data utils, Cox, signatures, stats (+ compareFeatures)
#  11_rtilandscape_utils.R        |    16     |    547       | RTI-specific: fusions, oncoprint, drug annot
#  -------------------------------|-----------|--------------|-----------------------------
#  TOTAL                          |   154     |   9483       |
#
# =============================================================================
# DUPLICATE FUNCTIONS — ANALYSIS AND RECOMMENDATIONS
# =============================================================================
#
# The following functions appear in multiple source files across the project.
# Each entry lists all locations, key differences, and the recommended version.
#
# ─────────────────────────────────────────────────────────────────────────────
# 1. filter_mutects  (SNV mutation filtering)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mutations.R  [2026-03-02] ← CANONICAL
#     B) Pembro/funcsInPembro.R                                              [2025-07-02]
#     C) IBC/DNAfuncs.R                                                      [2024-03-18]
#     D) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#   DIFFERENCES:
#     A (newest): Adds treat_noisy_samples / noisy_samples parameters to
#                 handle noisy samples explicitly. Otherwise API-identical to B.
#     B: Identical API to A, missing treat_noisy_samples; last tuned for Pembro study.
#     C: Older — missing treat_noisy_samples. Otherwise same filtering logic.
#     D: Oldest version; some filter defaults differ.
#   RECOMMENDATION: Use A (02_mutation_analysis.R). Contains the most parameters.
#
# ─────────────────────────────────────────────────────────────────────────────
# 2. filter_pindels  (Indel filtering)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS: Same 4 as filter_mutects (same provenance and age progression).
#   DIFFERENCES: Same pattern as filter_mutects — A adds treat_noisy_samples.
#   RECOMMENDATION: Use A (02_mutation_analysis.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 3. dge_limma  (Differential gene expression via limma)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/mRNA.R       [2026-03-02] ← CANONICAL
#     B) Pembro/funcsInPembro.R                                              [2025-07-02]
#     C) GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R                  [2025-08-11]
#     D) woodman_lab.XLi23/expr/R/dge_limma.R                                [2024-07-11]
#     E) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#   DIFFERENCES:
#     A: Adds robust_lm and weight parameters (robust linear model + observational weights).
#     B: Same as A (Pembro copy).
#     C: Missing robust_lm / weight params; slightly older.
#     D/E: Older implementations; parameter names differ (sparse_filter_cutoff vs
#          sample_frequency_threshold).
#   RECOMMENDATION: Use A (03_expression_analysis.R). Most complete parameter set.
#
# ─────────────────────────────────────────────────────────────────────────────
# 4. gsea  (Gene Set Enrichment Analysis)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) Pembro/funcsInPembro.R                                              [2025-07-02] ← CANONICAL
#     B) GBM/funcsInGBM.R                                                    [2024-09-06]
#     C) GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R                  [2025-08-11]
#     D) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#   DIFFERENCES:
#     A/B/C: Functionally identical signatures. A is most recently touched.
#     D: Older; some internal differences in output handling.
#   RECOMMENDATION: Use A (05_gsea_pathway.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 5. gsea-adjacent: pathway_gene_table, fgseatableplot, sankey_heatmap
# ─────────────────────────────────────────────────────────────────────────────
#   These three functions travel together and appear in the same 4 locations
#   as gsea above. Identical function bodies across all copies.
#   RECOMMENDATION: Use versions in 05_gsea_pathway.R (from Pembro/funcsInPembro.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 6. check_low_expression
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) Pembro/funcsInPembro.R                                              [2025-07-02] ← CANONICAL
#     B) GBM/funcsInGBM.R                                                    [2024-09-06]
#     C) Pembro/Pembro_codebase/funcsInPembro_codebase.R                    [2026-01-23]
#   DIFFERENCES: All three are byte-for-byte identical.
#   RECOMMENDATION: Use A (03_expression_analysis.R). No meaningful difference.
#
# ─────────────────────────────────────────────────────────────────────────────
# 7. run_dge_gsea_pipeline
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) Pembro/funcsInPembro.R                                              [2025-07-02] ← CANONICAL
#     B) Pembro/Pembro_codebase/funcsInPembro_codebase.R                    [2026-01-23]
#   DIFFERENCES: B adds gsea_rank_p_thresh parameter filtering logic.
#     However A has more complete visualization branches.
#   RECOMMENDATION: Use A (03_expression_analysis.R); check B if you need gsea_rank_p_thresh.
#
# ─────────────────────────────────────────────────────────────────────────────
# 8. immu_deconvolution
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) Pembro/funcsInPembro.R                                              [2025-07-02] ← CANONICAL
#     B) GBM/funcsInGBM.R                                                    [2024-09-06]
#     C) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#   DIFFERENCES: A/B identical. C is an older version.
#   RECOMMENDATION: Use A (07_immune_deconvolution.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 9. df2granges / df2grangelist / bingranges / gr2linear
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS: These 4 utility functions appear in 7+ files:
#     A) woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/functions.R  [2026-03-02] ← CANONICAL
#     B) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#     C) GBM/funcsInGBM.R                                                    [2024-09-06]
#     D) Pembro/funcsInPembro.R                                              [2025-07-02]
#     E) IBC/DNAfuncs.R                                                      [2024-03-18]
#     F) GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R                  [2025-08-11]
#   DIFFERENCES: All essentially identical. A has the most complete roxygen2 docs.
#   RECOMMENDATION: Use A (01_genomic_data.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 10. ggarms  (Chromosomal arm CNV visualization)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/myscripts.R                                       [2023-06-29] ← used in 01
#     B) GBM/funcsInGBM.R                                                    [2024-09-06]
#     C) IBC/DNAfuncs.R                                                      [2024-03-18]
#     D) Pembro/funcsInPembro.R                                              [2025-07-02]
#     E) GBM/GBM_Weathers.Rmd / GBM_figures.Rmd (multiple Rmd versions)
#   DIFFERENCES: B/C/D are essentially identical to A. Rmd versions are embedded
#     copies for self-contained notebooks.
#   RECOMMENDATION: Use version in 01_genomic_data.R.
#
# ─────────────────────────────────────────────────────────────────────────────
# 11. plotKM  (Kaplan-Meier survival plot)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R                  [2025-08-11] ← CANONICAL
#     B) GBM/funcsInGBM.R                                                    [2024-09-06]
#     C) GBM/XL_figures/results/TCGA_Glass_EIS_IDH1wt_KMcurves_BAM.R       [~2023]
#   DIFFERENCES:
#     A/B: Nearly identical; A adds timeYLab parameter.
#     C: Older simplified version — fewer parameters, no statOnly mode.
#   RECOMMENDATION: Use A (06_survival_analysis.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 12. estimate_bestNumberofClusters / map_clusters / *Cluster wrappers / multiCluster
#     consensusCluster / future_consensusCluster
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) GBM/XL_figures/vDotJun0823/GBMnat/R/funcsGBMnat.R                  [2025-08-11] ← CANONICAL
#     B) GBM/XL_figures/vDotJun0823/GBMnat_copy/funcsGBMnat.R               [2024-10-14]
#     C) woodman_lab.XLi23/expr/R/clustering.R                               [2024-07-11]
#     D) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#   DIFFERENCES:
#     A vs B: Files are IDENTICAL in content despite different dates (B is an
#             earlier copy; A was touched more recently).
#     A vs C: C lacks dge_limma and gsea functions that are mixed into funcsGBMnat.R.
#             Clustering logic itself is the same.
#     A vs D: D is older; some cluster wrappers have slightly different internals.
#   RECOMMENDATION: Use A (04_clustering.R).
#
# ─────────────────────────────────────────────────────────────────────────────
# 13. vectraStatsFun.R  (Spatial immune profiling statistics)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/mIFvectra/R/vectraStatsFun.R                      [2024-07-02] ← CANONICAL (general)
#     B) Updated_Panel_10/vectraStatsFun.R                                    [2024-08-05] ← Use for Panel 10 work
#     C) woodman_lab.XLi23/spatialImf/R/vectraStatsFun.R                     [2023-06-13]
#   DIFFERENCES:
#     A: defineCell uses na.omit(). Includes 3-panel cell definitions (Panel 10, 1, 2).
#        23 essential cell definitions. Most complete feature set.
#     B: defineCell omits na.omit(). Panel 10 only. 36 essential cell definitions
#        — more comprehensive for Panel 10-specific work. Newest (2024-08-05).
#     C: Oldest. defineCell has no na.omit(). No panel/cell definitions block.
#        Core stats functions identical to A.
#   RECOMMENDATION: Use A (08_spatial_immune_profiling.R) for general use;
#     use B if working specifically with Panel 10 and needing all 36 cell def variants.
#
# ─────────────────────────────────────────────────────────────────────────────
# 14. myscripts.R  (Large monolithic utility script)
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/myscripts.R                                       [2023-06-29]
#     B) woodman_lab.XLi23/HaifengScript/myscripts.R                        [2023-07-18]
#   DIFFERENCES: ~99% identical. B has minor documentation cleanup, some roxygen2
#     comments rearranged. No functional differences.
#   NOTE: The functions from myscripts.R have been superseded by the modular
#     HaifengPackages.Mar2026 package (functions.R + mRNA.R + mutations.R).
#   RECOMMENDATION: Prefer the HaifengPackages.Mar2026 modular files. myscripts.R
#     is useful as a legacy reference only.
#
# ─────────────────────────────────────────────────────────────────────────────
# 15. GBMnat/funcsGBMnat.R vs GBMnat_copy/funcsGBMnat.R
# ─────────────────────────────────────────────────────────────────────────────
#   These files are BYTE-FOR-BYTE IDENTICAL in function content.
#   GBMnat_copy/ was created as a working backup.
#   RECOMMENDATION: Use GBMnat/R/funcsGBMnat.R (modification date 2025-08-11).
#     GBMnat_copy/ can be deleted safely.
#
# ─────────────────────────────────────────────────────────────────────────────
# 16. HaifengPackages.Mar2026 vs bulkExpr/utiltools
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/HaifengPackages.Mar2026/utiltools/R/  [2026-03-02] ← CANONICAL
#     B) bulkExpr/utiltools/R/                                    [2025-02-17]
#   DIFFERENCES: The three files (mRNA.R, functions.R, mutations.R) are IDENTICAL
#     in content — bulkExpr has an older timestamp but same byte count.
#   RECOMMENDATION: Use A. bulkExpr/ appears to be a deployment copy.
#
# ─────────────────────────────────────────────────────────────────────────────
# 17. oncokb annotation: woodman_lab.XLi23/oncokbAnn vs RTI_landscape/R/RTI_funcs.R
# ─────────────────────────────────────────────────────────────────────────────
#   LOCATIONS:
#     A) woodman_lab.XLi23/oncokbAnn/R/oncokb.R                             [2026-01-24] ← CANONICAL
#     B) RTI_landscape/draftFeb/RTIlndscp/R/RTI_funcs.R                     [2026-04-07]
#   DIFFERENCES:
#     A: Cleaner API, uses httr2, proper token management, deduplication logic.
#        Organized as a proper R package with roxygen2 docs.
#     B: Older query_oncokb() with hardcoded API token in source (security risk!).
#        Also contains oncokb_annotate_maf_bypos but with slightly different internals.
#   RECOMMENDATION: Use A (09_oncokb_annotation.R). Remove hardcoded token from B.
#
# =============================================================================
# FUNCTION QUICK-REFERENCE
# =============================================================================
#
# GENOMIC DATA (01):
#   isar_transform, df2granges, df2grangelist, bingranges, gr2linear,
#   column2namedVector, align_df, chisq_dist,
#   ggcopynumber, ggscores, ggarms, gggenome
#
# MUTATION ANALYSIS (02):
#   get_combined_tcga_mutation_table, filter_mutects, filter_pindels,
#   convert_mutations, complex_oncoplot, complex_oncoplot_o,
#   annotate_genomic_mutations_with_panel_info,
#   annotate_mdl_mutations_with_timepoints,
#   annotate_mdl_mutations_with_mutationdb, get_min_multihit_distance,
#   harmonize_mutations, oncogenicpathway_tableplot, cancerhallmark_heatmap
#
# EXPRESSION ANALYSIS (03):
#   check_low_expression, prepare_unsupervised_data, dge_limma, dge_edgeR,
#   dge_DESeq, run_dge_gsea_pipeline, plot_volcano_with_annotations,
#   unsupervised_analysis, consensus_immunedeconvolute, geneset_activity
#
# CLUSTERING (04):
#   prepare_unsupervised_data, estimate_bestNumberofClusters, map_clusters,
#   kmeansCluster, pamCluster, hclustCluster, fuzzyCluster, mclustCluster,
#   apclustCluster, hdbscanCluster, mclCluster, speccCluster, kkmeansCluster,
#   skmeansCluster, nmfCluster, somCluster,
#   multiCluster, consensusCluster, future_consensusCluster
#
# GSEA & PATHWAY (05):
#   gsea, pathway_gene_table, fgseatableplot, sankey_heatmap,
#   read_pathways, goEnrich
#
# SURVIVAL (06):
#   plotKM, getKM_medium, uniCoxPh, multiCoxPh
#
# IMMUNE DECONVOLUTION (07):
#   DeconRNASeq_, CoreAlg, doPerm, CIBERSORT, immu_deconvolution
#
# SPATIAL IMMUNE PROFILING (08):
#   defineCell, reassignNA, getCtsPcDens, compute_all_nearest_distance,
#   find_nearest_distance_dist, unique_phenotypes, distance_matrix,
#   select_rows, getMeanNearestNeighbourDistance, count_within_many,
#   count_within, getCountWithin, performTtestsAllRows,
#   performPairedTtestsAllRows, performTtestsAllClassesOneVsRest,
#   performTtestsAllClassesEachPair, getAssociation, getCor, robustscale,
#   getSigPair, bi_ripleys_k_mod, getMeanCountWithin,
#   PieDonut, plotPieDonutBySlides, plot_immunoflo, trimMif4plots,
#   plotImmuneHighlights
#
# ONCOKB ANNOTATION (09):
#   `%||%`, .oncokb_get_token, .oncokb_request, get_oncokb_version,
#   extract_diagnostic_implications, extract_treatments,
#   query_oncokb_bypos, query_oncokb_CN, oncokb_annotate_maf,
#   oncokb_annotate_maf_bypos, query_oncokb_fusion,
#   oncokb_annotate_fusion, expand_oncokb_res_tibbles
#
# STATISTICS & UTILITIES (10):
#   is.empty.data.frame, getmode, removeoutlier,
#   set_column_as_rownames, discrete_quantile,
#   filter_assay_info, format_infotable, extract_paired_sample_info,
#   read_info, read_foundry_datasets, write_foundry_datasets,
#   uni_cox, cal_riskscore, survivalSignatures, predictGroup,
#   compare_two_groups, safe_compare_two_groups, add_p_adjust,
#   title_case_smart, compareFeatures
#
# RTI LANDSCAPE UTILS (11):
#   auto_convert_types, normalize_spaces, capitalize_words_with_underscores,
#   collapse_columns_dt,
#   cluster_only_within_group,
#   fit_normal_threshold, run_category_enrichment, recompute_ranks,
#   normalize_fus, getStrandInfo, score_fusions,
#   oncoprint_column_order, make_column_order, make_row_order,
#   map_drugs, pick_specific_pathway
#
# =============================================================================
