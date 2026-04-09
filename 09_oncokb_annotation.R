# =============================================================================
# ONCOKB API ANNOTATION FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-08
# =============================================================================
#
# CONTENTS:
#   %||%                        - Null-coalescing operator
#   .oncokb_get_token           - Internal: retrieve API token from env
#   .oncokb_request             - Internal: authenticated HTTP request
#   get_oncokb_version          - Get current OncoKB database version
#   extract_diagnostic_implications - Parse diagnostic/prognostic implications
#   extract_treatments          - Parse treatment information from response
#   query_oncokb_bypos          - Query by protein change or genomic position
#   query_oncokb_CN             - Query copy number alterations
#   oncokb_annotate_maf         - Annotate MAF by protein change (parallel)
#   oncokb_annotate_maf_bypos   - Annotate MAF by genomic coordinates (parallel)
#   query_oncokb_fusion         - Query fusion/structural variant events
#   oncokb_annotate_fusion      - Annotate fusion dataframe (parallel)
#   expand_oncokb_res_tibbles   - Expand nested OncoKB result tibble columns
#
# SOURCE PROVENANCE:
#   woodman_lab.XLi23/oncokbAnn/R/oncokb.R (2026-01-24) — canonical, most recent
#
# USAGE:
#   Sys.setenv(ONCOKB_API_TOKEN = "your_token_here")
#   result <- query_oncokb_bypos(gene="BRAF", protein_change="V600E")
#   maf_ann <- oncokb_annotate_maf_bypos(maf, cancer_types="BRCA")
#
# DEPENDENCIES: httr2, tibble, purrr, dplyr, rlang, furrr (for parallel), stringr
# =============================================================================

# --- INTERNAL HELPERS --------------------------------------------------------

#' Null-coalescing helper
#' @keywords internal
# SOURCE: oncokb.R (2026-01-24)
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Retrieve OncoKB API token from argument or ONCOKB_API_TOKEN env var
#' @keywords internal
# SOURCE: oncokb.R (2026-01-24)
.oncokb_get_token <- function(token = NULL) {
  if (!is.null(token)) return(token)
  tok <- Sys.getenv("ONCOKB_API_TOKEN", unset = NA_character_)
  if (is.na(tok) || tok == "") stop(
    "OncoKB API token not provided. Set ONCOKB_API_TOKEN or pass `token=` explicitly.",
    call. = FALSE
  )
  tok
}

#' Perform an authenticated HTTP request to OncoKB API
#' @keywords internal
# SOURCE: oncokb.R (2026-01-24)
.oncokb_request <- function(url, query = list(), token) {
  httr2::request(url) |>
    httr2::req_headers(
      Accept = "application/json",
      Authorization = sprintf("Bearer %s", token)
    ) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_perform() |>
    httr2::resp_body_json()
}

# --- VERSION -----------------------------------------------------------------

#' Get current OncoKB data version
#'
#' @param token OncoKB API token. If NULL, reads from ONCOKB_API_TOKEN env var.
#' @return Character scalar version string (e.g. "v4.14")
#' @export
# SOURCE: oncokb.R (2026-01-24)
get_oncokb_version <- function(token = NULL) {
  tok <- .oncokb_get_token(token)
  res <- .oncokb_request(
    url   = "https://www.oncokb.org/api/v1/info",
    token = tok
  )
  res$dataVersion$version %||% NA_character_
}

# --- PARSERS -----------------------------------------------------------------

#' Extract diagnostic or prognostic implications from OncoKB response
#'
#' @param oncokb_response Parsed JSON response from OncoKB
#' @param Implication Field name: "diagnosticImplications" or "prognosticImplications"
#' @return One-row tibble with tumorType, alterations, level, pmids
#' @export
# SOURCE: oncokb.R (2026-01-24)
extract_diagnostic_implications <- function(oncokb_response, Implication = "diagnosticImplications") {
  di <- oncokb_response[[Implication]]
  if (is.null(di) || length(di) == 0) {
    return(tibble::tibble(
      tumorType = NA_character_,
      alterations = NA_character_,
      level = NA_character_,
      pmids = NA_character_
    ))
  }

  tibble::tibble(
    tumorType   = paste(na.omit(purrr::map_chr(di, ~ .x$tumorType$name %||% NA_character_)), collapse = ";"),
    alterations = paste(na.omit(purrr::map_chr(di, ~ paste(.x$alterations %||% NA, collapse = ","))), collapse = ";"),
    level       = paste(na.omit(purrr::map_chr(di, ~ .x$levelOfEvidence %||% NA_character_)), collapse = ";"),
    pmids       = paste(na.omit(purrr::map_chr(di, function(x) {
      if (!is.null(x$pmids)) {
        paste(x$pmids, collapse = ",")
      } else if (!is.null(x$citations$pmids)) {
        paste(x$citations$pmids, collapse = ",")
      } else NA_character_
    })), collapse = ";")
  )
}

#' Extract treatment information from OncoKB response
#'
#' @param oncokb_response Parsed JSON response from OncoKB
#' @return List with elements \code{onerow} (collapsed tibble) and \code{detailed} (per-alteration tibble)
#' @export
# SOURCE: oncokb.R (2026-01-24)
extract_treatments <- function(oncokb_response) {
  tr <- oncokb_response$treatments
  if (is.null(tr) || length(tr) == 0) {
    return(list(
      onerow = tibble::tibble(
        alterations = NA_character_,
        level = NA_character_,
        fdaLevel = NA_character_,
        tumorType = NA_character_,
        drugs = NA_character_
      ),
      detailed = NULL
    ))
  }

  detailed <- purrr::map_dfr(tr, function(x) {
    tibble::tibble(
      alterations = paste(x$alterations %||% NA, collapse = "+"),
      drugs       = paste(unique(purrr::map_chr(x$drugs %||% list(), "drugName")), collapse = "+"),
      level       = x$level %||% NA_character_,
      fdaLevel    = x$fdaLevel %||% NA_character_,
      tumorType   = x$levelAssociatedCancerType$name %||% NA_character_
    )
  })

  onerow <- detailed |>
    dplyr::group_by(.data$alterations) |>
    dplyr::summarise(
      level     = paste(unique(na.omit(.data$level)), collapse = ","),
      fdaLevel  = paste(unique(na.omit(.data$fdaLevel)), collapse = ","),
      tumorType = paste(unique(na.omit(.data$tumorType)), collapse = ","),
      drugs     = paste(unique(na.omit(.data$drugs)), collapse = ","),
      .groups = "drop"
    ) |>
    dplyr::summarise(
      alterations = paste(.data$alterations, collapse = ";"),
      level       = paste(.data$level, collapse = ";"),
      fdaLevel    = paste(.data$fdaLevel, collapse = ";"),
      tumorType   = paste(.data$tumorType, collapse = ";"),
      drugs       = paste(.data$drugs, collapse = ";"),
      .groups = "drop"
    )

  list(onerow = onerow, detailed = detailed)
}

# --- QUERY FUNCTIONS ---------------------------------------------------------

#' Query OncoKB by protein change or genomic position
#'
#' @param gene Hugo gene symbol
#' @param protein_change Protein alteration string (e.g. "V600E")
#' @param variant_type Variant consequence (e.g. "Missense_Mutation")
#' @param chrom Chromosome (used if protein_change is NULL)
#' @param start Start position
#' @param end End position
#' @param ref Reference allele
#' @param alt Alternate allele
#' @param cancer_type OncoTree tumor type code (optional)
#' @param token OncoKB API token
#' @param referenceGenome "GRCh37" or "GRCh38"
#' @param oncokb_version OncoKB data version string
#' @return Tibble with oncogenic classification, treatment, diagnostic, and prognostic fields
#' @export
# SOURCE: oncokb.R (2026-01-24)
query_oncokb_bypos <- function(
    gene = NULL,
    protein_change = NULL,
    variant_type = NULL,
    chrom = NULL,
    start = NULL,
    end = NULL,
    ref = NULL,
    alt = NULL,
    cancer_type = NA,
    token = NULL,
    referenceGenome = "GRCh37",
    oncokb_version = get_oncokb_version()
) {
  tok <- .oncokb_get_token(token)
  if (is.null(cancer_type) || is.na(cancer_type)) cancer_type <- NULL

  if (!is.null(protein_change) && protein_change != '') {
    url <- "https://www.oncokb.org/api/v1/annotate/mutations/byProteinChange"
    query <- list(
      hugoSymbol = gene,
      alteration = protein_change,
      consequence = variant_type,
      tumorType = cancer_type
    )
  } else if (!is.null(chrom) && !is.null(start) && !is.null(end) &&
             !is.null(ref) && !is.null(alt)) {
    url <- "https://www.oncokb.org/api/v1/annotate/mutations/byGenomicChange"
    genomicLocation <- paste(chrom, start, end, ref, alt, sep = ",")
    query <- list(
      genomicLocation = genomicLocation,
      referenceGenome = referenceGenome,
      tumorType = cancer_type
    )
  } else return(tibble::tibble(oncogenic = ''))

  res <- .oncokb_request(url, query, tok)

  tibble::tibble(
    hugoSymbol = res$query$hugoSymbol %||% NA,
    alteration = res$query$alteration %||% NA,
    consequence = res$query$consequence %||% NA,
    geneExist = res$geneExist %||% NA,
    variantExist = res$variantExist %||% NA,
    alleleExist = res$alleleExist %||% NA,
    oncogenic = as.character(res$oncogenic),
    mutationEffect = res$mutationEffect$knownEffect %||% NA,
    hotspot = res$hotspot %||% NA,
    exon = ifelse(is.null(res$exon), NA, paste(res$exon, collapse = ';')),
    highestDiagnosticImplicationLevel = res$highestDiagnosticImplicationLevel %||% NA,
    highestPrognosticImplicationLevel = res$highestPrognosticImplicationLevel %||% NA,
    diagnosticImplications = extract_diagnostic_implications(res),
    prognosticImplications = extract_diagnostic_implications(res, Implication = "prognosticImplications"),
    highestSensitiveLevel = res$highestSensitiveLevel %||% NA,
    highestResistanceLevel = res$highestResistanceLevel %||% NA,
    treatment = extract_treatments(res)$onerow,
    oncokb_version = oncokb_version
  )
}

#' Query OncoKB for copy number alterations
#'
#' @param gene Hugo symbol
#' @param copy_number_type "AMPLIFICATION" or "DELETION"
#' @param cancer_type OncoTree code (optional)
#' @param token OncoKB API token
#' @param referenceGenome Genome build
#' @param oncokb_version OncoKB version string
#' @return Tibble with oncogenicity, treatment, diagnostic/prognostic fields
#' @export
# SOURCE: oncokb.R (2026-01-24)
query_oncokb_CN <- function(
    gene,
    copy_number_type = c("AMPLIFICATION", "DELETION"),
    cancer_type = NULL,
    token = NULL,
    referenceGenome = "GRCh37",
    oncokb_version = get_oncokb_version()
) {
  copy_number_type <- match.arg(toupper(copy_number_type))
  tok <- .oncokb_get_token(token)
  url <- "https://www.oncokb.org/api/v1/annotate/copyNumberAlterations"
  query <- list(
    hugoSymbol = gene,
    copyNameAlterationType = copy_number_type,
    referenceGenome = referenceGenome,
    tumorType = cancer_type
  )
  res <- .oncokb_request(url, query, tok)

  tibble::tibble(
    hugoSymbol = res$query$hugoSymbol %||% gene,
    oncogenic = as.character(res$oncogenic %||% NA),
    mutationEffect = res$mutationEffect$knownEffect %||% NA,
    geneExist = res$geneExist %||% NA,
    highestSensitiveLevel = res$highestSensitiveLevel %||% NA,
    highestResistanceLevel = res$highestResistanceLevel %||% NA,
    highestDiagnosticImplicationLevel = res$highestDiagnosticImplicationLevel %||% NA,
    highestPrognosticImplicationLevel = res$highestPrognosticImplicationLevel %||% NA,
    treatment = extract_treatments(res)$onerow,
    diagnosticImplications = extract_diagnostic_implications(res),
    prognosticImplications = extract_diagnostic_implications(res, Implication = "prognosticImplications"),
    oncokb_version = oncokb_version
  )
}

#' Annotate MAF with OncoKB using protein change
#'
#' Requires Hugo_Symbol and Protein_Change columns. Uses consequence_map for
#' Variant_Classification to consequence conversion.
#'
#' @param maf data.frame in MAF format
#' @param cancer_types Optional cancer type vector (or NULL to use maf$cancer_type)
#' @param parallelize Logical, use furrr::future_pmap_dfr (default TRUE)
#' @return Input MAF with OncoKB annotation columns appended
#' @export
# SOURCE: oncokb.R (2026-01-24)
oncokb_annotate_maf <- function(maf, cancer_types = NULL, parallelize = TRUE) {
  if (is.null(cancer_types) & !'cancer_type' %in% names(maf)) {
    message('No cancer type(s) specified, defaults to NA')
    maf$cancer_type <- NA
  } else if (is.character(cancer_types)) {
    maf <- tibble::add_column(maf, cancer_type = cancer_types)
  }

  oncokb_cols <- dplyr::mutate(maf,
                               gene = Hugo_Symbol,
                               protein_change = stringr::str_replace(Protein_Change, 'p.', ''),
                               variant_type = consequence_map[Variant_Classification]) %>%
    dplyr::mutate(
      start = stringr::str_extract(protein_change, '^[0-9]+'),
      end   = stringr::str_extract(protein_change, '[0-9]+[a-zA-Z]+$') %>% stringr::str_extract("[0-9]+")
    ) %>%
    dplyr::select(gene, protein_change, variant_type, start, end, cancer_type)

  if (parallelize) {
    oncokb_cols <- furrr::future_pmap_dfr(oncokb_cols, query_oncokb_bypos)
  } else {
    oncokb_cols <- purrr::pmap_dfr(oncokb_cols, query_oncokb_bypos)
  }

  dplyr::bind_cols(maf, oncokb_cols)
}

#' Annotate MAF with OncoKB using genomic coordinates
#'
#' Uses Chromosome, Start_Position, End_Position, Reference_Allele, Tumor_Seq_Allele2.
#' Deduplicates queries to avoid redundant API calls.
#'
#' @param maf data.frame in MAF format
#' @param cancer_types Optional cancer type vector (or NULL to use maf$cancer_type)
#' @param parallelize Logical, use furrr (default FALSE — safer for large MAFs)
#' @return Input MAF with OncoKB annotation columns appended
#' @export
# SOURCE: oncokb.R (2026-01-24)
oncokb_annotate_maf_bypos <- function(maf, cancer_types = NULL, parallelize = F) {
  if (is.null(cancer_types) & !'cancer_type' %in% names(maf)) {
    message('No cancer type(s) specified, defaults to NA')
    maf$cancer_type <- NA_character_
  } else if (is.character(cancer_types)) {
    maf <- tibble::add_column(maf, cancer_type = cancer_types)
  }

  queries <- maf %>%
    dplyr::mutate(
      chrom = Chromosome,
      start = Start_Position,
      end   = End_Position,
      ref   = Reference_Allele,
      alt   = Tumor_Seq_Allele2,
      cancer_type = cancer_type
    ) %>%
    dplyr::select(chrom, start, end, ref, alt, cancer_type)

  queries$key <- paste(queries$chrom, queries$start, queries$end,
                       queries$ref, queries$alt, queries$cancer_type, sep = "_")

  unique_queries <- queries %>% dplyr::distinct(key, .keep_all = TRUE)

  if (parallelize) {
    results <- furrr::future_pmap_dfr(unique_queries %>% dplyr::select(-key), query_oncokb_bypos)
  } else {
    results <- purrr::pmap_dfr(unique_queries %>% dplyr::select(-key), query_oncokb_bypos)
  }

  results$key <- unique_queries$key

  annotated <- dplyr::left_join(queries, results, by = "key")

  dplyr::bind_cols(maf, annotated %>% dplyr::select(-key, -chrom, -start, -end, -ref, -alt, -cancer_type))
}

#' Query OncoKB for fusion / structural variant events
#'
#' @param geneA Hugo symbol of 5' gene
#' @param geneB Hugo symbol of 3' gene
#' @param structuralVariantType Default "FUSION"
#' @param isFunctionalFusion Logical
#' @param cancer_type OncoTree code (optional)
#' @param token OncoKB API token
#' @param referenceGenome Genome build
#' @param oncokb_version OncoKB version string
#' @return Tibble with oncogenicity, treatment, and implication fields
#' @export
# SOURCE: oncokb.R (2026-01-24)
query_oncokb_fusion <- function(
    geneA, geneB, structuralVariantType = "FUSION", isFunctionalFusion = FALSE,
    cancer_type = NULL, token = NULL, referenceGenome = "GRCh37",
    oncokb_version = get_oncokb_version()
) {
  tok <- .oncokb_get_token(token)
  url <- "https://www.oncokb.org/api/v1/annotate/structuralVariants"
  query <- list(
    hugoSymbolA = geneA,
    hugoSymbolB = geneB,
    structuralVariantType = toupper(structuralVariantType),
    isFunctionalFusion = tolower(isFunctionalFusion),
    referenceGenome = referenceGenome,
    tumorType = cancer_type
  )
  res <- .oncokb_request(url, query, tok)

  tibble::tibble(
    hugoSymbol = res$query$hugoSymbol %||% NA,
    alteration = res$query$alteration %||% NA,
    geneExist = res$geneExist %||% NA,
    variantExist = res$variantExist %||% NA,
    alleleExist = res$alleleExist %||% NA,
    oncogenic = as.character(res$oncogenic %||% NA),
    hotspot = res$hotspot %||% NA,
    mutationEffect = res$mutationEffect$knownEffect %||% NA,
    highestSensitiveLevel = res$highestSensitiveLevel %||% NA,
    highestResistanceLevel = res$highestResistanceLevel %||% NA,
    highestDiagnosticImplicationLevel = res$highestDiagnosticImplicationLevel %||% NA,
    highestPrognosticImplicationLevel = res$highestPrognosticImplicationLevel %||% NA,
    treatment = extract_treatments(res)$onerow,
    diagnosticImplications = extract_diagnostic_implications(res),
    prognosticImplications = extract_diagnostic_implications(res, Implication = "prognosticImplications"),
    oncokb_version = oncokb_version
  )
}

#' Annotate a fusion dataframe with OncoKB
#'
#' Deduplicates queries and supports optional parallel processing.
#'
#' @param fusion_df Dataframe with fusion calls
#' @param geneA_col Column name for 5' gene (default "gene1")
#' @param geneB_col Column name for 3' gene (default "gene2")
#' @param functional_col Column name for functional status (default "functional")
#' @param cancer_type_col Optional column name for cancer type
#' @param token OncoKB API token
#' @param parallelize Logical (default TRUE)
#' @param referenceGenome Genome build
#' @param oncokb_version OncoKB version string
#' @return fusion_df with OncoKB annotation columns appended
#' @export
# SOURCE: oncokb.R (2026-01-24)
oncokb_annotate_fusion <- function(
    fusion_df, geneA_col = "gene1", geneB_col = "gene2",
    functional_col = "functional", cancer_type_col = NULL,
    token = NULL, parallelize = TRUE, referenceGenome = "GRCh37",
    oncokb_version = get_oncokb_version()
) {
  queries <- fusion_df %>%
    dplyr::mutate(
      geneA = .data[[geneA_col]],
      geneB = .data[[geneB_col]],
      isFunctionalFusion = if (!is.null(functional_col)) .data[[functional_col]] else FALSE,
      cancer_type = if (!is.null(cancer_type_col) && cancer_type_col %in% names(.)) .data[[cancer_type_col]] else NA_character_
    ) %>%
    dplyr::select(geneA, geneB, isFunctionalFusion, cancer_type)

  queries$key <- paste(queries$geneA, queries$geneB, queries$isFunctionalFusion, queries$cancer_type, sep = "_")
  unique_queries <- queries %>% dplyr::distinct(key, .keep_all = TRUE)

  query_fun <- function(geneA, geneB, isFunctionalFusion, cancer_type) {
    query_oncokb_fusion(
      geneA = geneA, geneB = geneB,
      isFunctionalFusion = isFunctionalFusion,
      cancer_type = cancer_type, token = token,
      referenceGenome = referenceGenome,
      oncokb_version = oncokb_version
    )
  }

  if (parallelize) {
    results <- furrr::future_pmap_dfr(unique_queries %>% dplyr::select(-key), query_fun)
  } else {
    results <- purrr::pmap_dfr(unique_queries %>% dplyr::select(-key), query_fun)
  }

  results$key <- unique_queries$key
  annotated <- dplyr::left_join(queries, results, by = "key")

  dplyr::bind_cols(fusion_df, annotated %>% dplyr::select(-key, -geneA, -geneB, -isFunctionalFusion, -cancer_type))
}

#' Expand nested OncoKB result tibble columns
#'
#' Flattens treatment, diagnosticImplications, and prognosticImplications list-columns
#' into flat named columns prefixed with "treatment.".
#'
#' @param tmp Tibble from query_oncokb_bypos / query_oncokb_CN / query_oncokb_fusion
#' @return Flat data.frame with expanded columns
#' @export
# SOURCE: oncokb.R (2026-01-24)
expand_oncokb_res_tibbles=function(tmp){
  tmp.treatment=tmp$treatment
  tmp.treatment$.groups=NULL
  colnames(tmp.treatment)=paste0("treatment.",colnames(tmp.treatment))
  tmp.diagnosticImplications=tmp$diagnosticImplications
  tmp.diagnosticImplications$pmids=NULL
  colnames(tmp.diagnosticImplications)=paste0(colnames(tmp.diagnosticImplications))
  tmp.prognosticImplications=tmp$prognosticImplications
  tmp.prognosticImplications$pmids=NULL
  colnames(tmp.prognosticImplications)=paste0(colnames(tmp.prognosticImplications))
  tmp=tmp%>%
    dplyr::select(-c(treatment,diagnosticImplications,prognosticImplications))%>%
    cbind.data.frame(.,
                     diagnosticImplications=tmp.diagnosticImplications,
                     prognosticImplications=tmp.prognosticImplications,
                     tmp.treatment
    )
}
