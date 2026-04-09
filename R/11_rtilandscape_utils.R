# =============================================================================
# RTI LANDSCAPE UTILITY FUNCTIONS
# Consolidated library: WoodmanLab
# Generated: 2026-04-09
# =============================================================================
#
# CONTENTS — Data wrangling:
#   auto_convert_types              - Auto-detect and convert df column types
#   normalize_spaces                - Normalize whitespace/underscore in strings
#   capitalize_words_with_underscores - Smart capitalize with underscore delimiters
#   collapse_columns_dt             - Collapse repeated columns by key (data.table)
#
# CONTENTS — Clustering:
#   cluster_only_within_group       - Hierarchical clustering within group strata
#
# CONTENTS — Statistical helpers:
#   fit_normal_threshold            - Fit normal distribution and return threshold
#   run_category_enrichment         - Fisher/chi-sq enrichment of variant x category
#   recompute_ranks                 - Recompute prediction ranks per patient/model
#
# CONTENTS — Fusion analysis:
#   normalize_fus                   - Normalize fusion gene pair (sorted, "--" sep)
#   getStrandInfo                   - Get strand info for gene symbols via AnnotationDbi
#   score_fusions                   - Score fusion calls by frame/domain/read support
#
# CONTENTS — Oncoprint / column ordering:
#   oncoprint_column_order          - Score-based column ordering for oncoprints
#   make_column_order               - Group-aware column ordering (wraps above)
#   make_row_order                  - Row ordering by drug class + mutation burden
#
# CONTENTS — Drug annotation:
#   map_drugs                       - Map drug names to canonical synonyms
#   pick_specific_pathway           - Pick most specific pathway from a vector
#
# NOTE: prep_drug_heatmap() from RTI_funcs.R is NOT included here because it
#   depends on undeclared helpers (refactor_annotation_df, get_anno, refactor_track)
#   that are defined elsewhere in the RTI analysis context.
#
# SOURCE PROVENANCE:
#   RTI_landscape/draftFeb/RTIlndscp/R/RTI_funcs.R  (2026-04-07)
#
# =============================================================================


# --- SECTION 1: DATA WRANGLING -----------------------------------------------

# SOURCE: RTI_funcs.R (2026-04-07)
# Auto-detect column types (numeric > date > logical > factor > character)
auto_convert_types <- function(df) {
  df[] <- lapply(df, function(col) {
    # Try numeric
    num_col <- suppressWarnings(as.numeric(col))
    if (!any(is.na(num_col)) || sum(is.na(num_col)) < length(num_col) / 2) {
      return(num_col)
    }

    # Try Date
    date_col <- suppressWarnings(as.Date(col, format = "%Y-%m-%d"))
    if (!all(is.na(date_col))) {
      return(date_col)
    }

    # Try logical
    if (all(tolower(col) %in% c("true", "false", "t", "f", "yes", "no", "1", "0"))) {
      return(as.logical(tolower(col) %in% c("true", "t", "yes", "1")))
    }

    # Default to factor if few unique values
    if (length(unique(col)) < length(col) / 2) {
      return(as.factor(col))
    }

    # Otherwise, keep as character
    return(as.character(col))
  })
  return(df)
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Normalize whitespace: replace non-breaking spaces, collapse multiple spaces,
# trim, then convert spaces to underscores.
normalize_spaces <- function(x) {
  x <- gsub("\u00A0", " ", x)         # Replace non-breaking spaces
  x <- gsub("_", " ", x)              # Replace underscores with spaces
  x <- gsub("\\s+", " ", x)           # Collapse multiple spaces
  x <- trimws(x)                      # Trim leading/trailing spaces
  x <- gsub(" ", "_", x)
  return(x)
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Capitalize words separated by underscores. Preserves all-caps acronyms (e.g. EGFR),
# hyphenated terms, and parentheses.
capitalize_words_with_underscores <- function(input_vector) {
  input_vector <- normalize_spaces(input_vector)

  sapply(input_vector, function(input_string) {
    if (is.na(input_string)) return(NA_character_)

    parts <- unlist(strsplit(input_string, "_"))

    transformed_parts <- sapply(parts, function(part) {
      tokens <- unlist(strsplit(part, "(?=[() ])|(?<=[() ])", perl = TRUE))

      tokens <- sapply(tokens, function(tok) {
        core <- gsub("[^A-Za-z]", "", tok)

        if (nchar(core) > 1 && core == toupper(core)) {
          tok  # keep acronyms (MEK, EGFR, etc.)
        } else if (grepl("^[() ]+$", tok)) {
          tok  # keep delimiters as-is
        } else if (grepl("-", tok)) {
          tok  # keep hyphenated terms as-is
        } else {
          stringr::str_to_title(tok)
        }
      }, USE.NAMES = FALSE)

      paste0(tokens, collapse = "")
    }, USE.NAMES = FALSE)

    paste(transformed_parts, collapse = "_")
  }, USE.NAMES = FALSE)
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Collapse specified columns within a data.table by key columns, joining unique
# non-NA values with `sep`. Accepts column names or a regex pattern.
collapse_columns_dt <- function(data, collapse_cols, sep = ";") {
  dt <- data.table::as.data.table(data)
  all_names <- names(dt)

  # Resolve collapse columns (accept regex or vector of names)
  if (length(collapse_cols) == 1 && !collapse_cols %in% all_names) {
    cols_to_collapse <- grep(collapse_cols, all_names, value = TRUE, ignore.case = TRUE)
  } else {
    cols_to_collapse <- collapse_cols
  }
  if (length(cols_to_collapse) == 0) {
    stop("No columns matched 'collapse_cols'.")
  }
  key_cols <- setdiff(all_names, cols_to_collapse)

  result_dt <- dt[, lapply(.SD, function(x) {
    vals <- x[!is.na(x) & x != ""]
    vals <- vals[!duplicated(vals)]
    if (length(vals) == 0) {
      NA_character_
    } else {
      paste(vals, collapse = sep)
    }
  }),
  .SDcols = cols_to_collapse,
  by = key_cols]

  # Clean spaces around separators
  for (col in cols_to_collapse) {
    result_dt[[col]] <- gsub(paste0("\\s*", sep, "\\s*"), sep, result_dt[[col]])
  }
  result_dt[]
}


# --- SECTION 2: CLUSTERING ---------------------------------------------------

# SOURCE: RTI_funcs.R (2026-04-07)
# Hierarchical clustering within group strata. Returns dendrogram(s) ordered
# within each group. Supports euclidean/jaccard/cor/bicor distance and all
# standard hclust linkage methods. Set merge_dendrogram=TRUE to merge group
# dendrograms into a single tree.
cluster_only_within_group = function(
    data, group_info, merge_dendrogram = FALSE,
    simMethods = c("euclidean", "maximum", "manhattan", "canberra", "minkowski",
                   "binary", "jaccard", "cor", "bicor"),
    clusterMethods = c("complete", "ward.D", "ward.D2", "single", "average",
                       "mcquitty", "median", "centroid")) {

  require(dendextend)
  if (missing(data))
    stop("data should be a matrix or data frame.")
  if (missing(simMethods)) {
    simMethods = match.arg(simMethods)
  } else {
    if (!simMethods %in% c("euclidean", "maximum", "manhattan", "canberra", "minkowski",
                           "binary", "jaccard", "cor", "bicor")) {
      stop("simMethods should be one of 'euclidean', 'maximum', 'manhattan', 'canberra', 'minkowski', 'binary', 'jaccard', 'cor' or 'bicor'.")
    } else {
      if (simMethods == "cor") {
        require(WGCNA)
        sim_mx = function(group_data) {
          cor_matrix <- 1 - WGCNA::corFast(t(group_data), use = "pairwise.complete.obs")
          as.dist(cor_matrix)
        }
      } else if (simMethods == "bicor") {
        require(WGCNA)
        sim_mx = function(group_data) {
          cor_matrix = 1 - bicor(group_data, use = "pairwise.complete.obs", nThreads = 4)
          as.dist(cor_matrix)
        }
      } else {
        sim_mx = function(group_data) dist(t(group_data), method = simMethods)
      }
    }
  }
  if (missing(clusterMethods)) {
    clusterMethods = match.arg(clusterMethods)
  }

  if (!is.factor(group_info)) {
    group_info = factor(group_info)
    group_levels = levels(group_info)
  } else {
    group_levels = levels(group_info)
  }

  group_dend = list()
  group_datas = list()
  order_list = list()

  for (i in 1:length(group_levels)) {
    group_data = data[, group_info == group_levels[i], drop = FALSE]

    gene_variances <- apply(group_data, 1, var, na.rm = TRUE)
    non_zero_var_genes <- which(gene_variances > 0)

    if (length(non_zero_var_genes) < 2) {
      warning(paste("Group '", group_levels[i], "' has no genes with variance > 0. Skipping clustering.", sep = ""))
      next
    }

    group_data_filtered <- group_data[non_zero_var_genes, ]
    group_datas[[i]] = group_data_filtered
    group_data = group_data_filtered

    if (ncol(group_data) > 1) {
      group_dist = sim_mx(group_data)
      group_hc = hclust(group_dist, method = clusterMethods)
      group_dend[[i]] = as.dendrogram(group_hc)
      order_list[[i]] = which(group_info == group_levels[i])[order.dendrogram(group_dend[[i]])]
      order.dendrogram(group_dend[[i]]) = order_list[[i]]
    } else {
      group_dend[[i]] = structure(which(group_info == group_levels[i]),
                                  class = "dendrogram", leaf = TRUE, height = 0,
                                  label = 1, members = 1)
    }
    attr(group_dend[[i]], ".class_label") = group_levels[i]
  }

  if (merge_dendrogram) {
    parent_dist = sim_mx(sapply(order_list, function(x) rowMeans(data[, x, drop = FALSE])))
    dend_p = as.dendrogram(hclust(parent_dist, method = clusterMethods))
    dend_m = merge_dendrogram(dend_p, group_dend)
    order.dendrogram(dend_m) = unlist(order_list[order.dendrogram(dend_p)])
  } else {
    dend_m = group_dend
  }

  return(list(dend = dend_m, order_list = order_list))
}


# --- SECTION 3: STATISTICAL HELPERS ------------------------------------------

# SOURCE: RTI_funcs.R (2026-04-07)
# Fit a normal distribution to x and return the threshold at a given quantile.
# direction="left" returns the lower tail threshold; "right" returns upper tail.
# Optionally pass flag_values to get a logical flag vector relative to threshold.
fit_normal_threshold <- function(x, quantile_cutoff = 0.05, direction = "left", flag_values = NULL) {
  x <- x[!is.na(x)]

  mu <- mean(x)
  sigma <- sd(x)

  if (direction == "left") {
    threshold <- qnorm(quantile_cutoff, mean = mu, sd = sigma)
  } else if (direction == "right") {
    threshold <- qnorm(1 - quantile_cutoff, mean = mu, sd = sigma)
  } else {
    stop("direction must be 'left' or 'right'")
  }

  flags <- NULL
  if (!is.null(flag_values)) {
    if (direction == "left") {
      flags <- flag_values < threshold
    } else {
      flags <- flag_values > threshold
    }
  }

  return(list(
    mean = mu,
    sd = sigma,
    threshold = threshold,
    flags = flags
  ))
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Fisher exact + chi-squared enrichment of variant x category combinations.
# Returns per-combination contingency results with BH-adjusted Fisher p.
run_category_enrichment <- function(
  variant_df,
  sample_info_df,
  category_col,
  variant_col,
  sample_id_col,
  patient_id_col
) {
  total_patients <- sample_info_df %>%
    dplyr::pull(.data[[patient_id_col]]) %>%
    unique() %>%
    length()

  variant_category_combos <- variant_df %>%
    dplyr::group_by(.data[[variant_col]], .data[[category_col]]) %>%
    dplyr::summarise(
      n = dplyr::n_distinct(.data[[sample_id_col]]),
      .groups = "drop"
    ) %>%
    dplyr::distinct()

  results <- vector("list", nrow(variant_category_combos))

  for (i in seq_len(nrow(variant_category_combos))) {
    variant_val  <- variant_category_combos[[variant_col]][i]
    category_val <- variant_category_combos[[category_col]][i]

    category_total <- sample_info_df %>%
      dplyr::filter(.data[[category_col]] == category_val) %>%
      dplyr::pull(.data[[patient_id_col]]) %>%
      unique() %>%
      length()

    variant_total <- variant_df %>%
      dplyr::filter(.data[[variant_col]] == variant_val) %>%
      dplyr::pull(.data[[sample_id_col]]) %>%
      unique() %>%
      length()

    a <- variant_category_combos$n[i]
    b <- category_total - a
    c <- variant_total - a
    d <- total_patients - (a + b + c)

    if (min(c(a, b, c, d)) < 0) next

    contingency <- matrix(c(a, b, c, d), nrow = 2)
    fisher_res <- tryCatch(fisher.test(contingency), error = function(e) NULL)
    chisq_p    <- tryCatch(chisq.test(contingency)$p.value, error = function(e) NA)

    results[[i]] <- tibble::tibble(
      variant      = variant_val,
      category     = category_val,
      var_in_cat_n = a,
      var_total_n  = variant_total,
      cat_total_n  = category_total,
      fisher_p     = if (!is.null(fisher_res)) fisher_res$p.value else NA,
      fisher_or    = if (!is.null(fisher_res)) fisher_res$estimate else NA,
      chisq_p      = chisq_p
    )
  }

  dplyr::bind_rows(results) %>%
    dplyr::mutate(fisher_p_adj = p.adjust(fisher_p, method = "BH"))
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Recompute prediction ranks within Patient_ID x model groups (ties = "min").
recompute_ranks <- function(df) {
  df %>%
    dplyr::group_by(Patient_ID, model) %>%
    dplyr::mutate(rank = rank(pred, ties.method = "min")) %>%
    dplyr::ungroup()
}


# --- SECTION 4: FUSION ANALYSIS ----------------------------------------------

# SOURCE: RTI_funcs.R (2026-04-07)
# Normalize a fusion pair to a canonical form: sort gene symbols alphabetically
# and join with "--". Used to deduplicate A--B vs B--A.
normalize_fus <- function(gene1, gene2) {
  sort(c(gene1, gene2)) %>% paste(collapse = "--")
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Get strand orientation for a vector of gene symbols using AnnotationDbi.
# Returns a named character vector (symbol -> strand "+"/"-").
getStrandInfo <- function(genes, gene_info) {
  map_ids <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = genes,
    keytype = "SYMBOL",
    columns = c("ENTREZID")
  )
  map_ids <- map_ids[map_ids$ENTREZID %in% gene_info$gene_id, ]
  map_ids <- map_ids[!duplicated(map_ids$SYMBOL), ]

  gene_info_subset <- gene_info[map_ids$ENTREZID]
  strand_info <- setNames(as.character(strand(gene_info_subset)), map_ids$SYMBOL)
  return(strand_info)
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Score fusion calls using a weighted scheme: frame (InFrame/UTR3/OutFrame),
# domain retention, split-read support, multi-sample evidence, cancer DB hits,
# TCGA hotspot genes, oncogene annotation, and multiple-mapping penalty.
# Returns the input data frame with added scoring columns and a composite `score`.
score_fusions <- function(fusion_df,
                          w_frame_in = 3,
                          w_frame_utr3 = 1,
                          w_domain = 2,
                          nread_support1 = 7,
                          nread_support2 = 15,
                          w_split1 = 1,
                          w_split2 = 2,
                          w_nsample = 1,
                          w_db = 2,
                          w_tcga_hot = 1,
                          w_oncogene = 1,
                          penalty_multimap = -1,
                          penalize_outframe = 0) {

  df <- fusion_df

  df$has_domain <- ifelse(is.na(df$KeptDomains1), FALSE,
                          ifelse(df$KeptDomains1 != "." & df$KeptDomains1 != "", TRUE, FALSE)) |
                   ifelse(is.na(df$KeptDomains2), FALSE,
                          ifelse(df$KeptDomains2 != "." & df$KeptDomains2 != "", TRUE, FALSE))

  df$frame_points <- 0
  df$frame_points[df$Frame == "InFrame"]  <- w_frame_in
  df$frame_points[df$Frame == "UTR3-CD"]  <- w_frame_utr3
  df$frame_points[df$Frame == "OutFrame"] <- penalize_outframe

  df$domain_points <- ifelse(df$has_domain, w_domain, 0)

  df$split_points <- 0
  df$split_points[df$split_cnt >= nread_support1] <- w_split1
  df$split_points[df$split_cnt >= nread_support2] <- w_split2

  df$nsample_points <- ifelse(df$n_sample > 1, w_nsample, 0)

  df$db_points <- ifelse(!is.na(df$cancer_db_hits) & df$cancer_db_hits != "", w_db, 0)

  df$tcga_points <- (ifelse(df$TCGA_Gene1_isHot == "Yes", w_tcga_hot, 0) +
                     ifelse(df$TCGA_Gene2_isHot == "Yes", w_tcga_hot, 0))

  df$oncogene_points <- 0
  df$oncogene_points[grepl("Oncogene", df$Genes1_Fun, ignore.case = TRUE)] <-
    df$oncogene_points[grepl("Oncogene", df$Genes1_Fun, ignore.case = TRUE)] + w_oncogene
  df$oncogene_points[grepl("Oncogene", df$Genes2_Fun, ignore.case = TRUE)] <-
    df$oncogene_points[grepl("Oncogene", df$Genes2_Fun, ignore.case = TRUE)] + w_oncogene

  df$multimap_points <- 0
  df$multimap_points[df$multiple_mapping > 1 |
                     df$Genes1_multiple_mapping == 1 |
                     df$Genes2_multiple_mapping == 1] <- penalty_multimap

  df$score <- rowSums(df[, c("frame_points", "domain_points", "split_points", "nsample_points",
                              "db_points", "tcga_points", "oncogene_points", "multimap_points")],
                      na.rm = TRUE)
  return(df)
}


# --- SECTION 5: ONCOPRINT / COLUMN ORDERING ----------------------------------

# SOURCE: RTI_funcs.R (2026-04-07)
# Score-based column ordering for oncoprints. Scores columns by the presence
# of TRUE values weighted by row position (higher rows = higher weight).
oncoprint_column_order = function(mat, row_order) {
  scoreCol = function(x) {
    score = 0
    for (i in 1:length(x)) {
      if (x[i]) {
        score = score + 2^(length(x) - i * 1 / x[i])
      }
    }
    score
  }
  scores = apply(mat[row_order, , drop = FALSE], 2, scoreCol)
  order(scores, decreasing = TRUE)
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Group-aware column ordering: orders columns within each group defined by
# order_by_col, respecting factor levels. Within each group, calls
# oncoprint_column_order(). Returns a vector of column names in display order.
make_column_order <- function(count_matrix, row_order, patient_meta,
                              order_by_col = "Diagnosis_II", id_col = "Patient_ID") {
  ord <- c()

  if (is.factor(patient_meta[[order_by_col]])) {
    groups <- levels(patient_meta[[order_by_col]])
  } else {
    groups <- unique(patient_meta[[order_by_col]])
  }

  for (grp in groups) {
    patients_in_grp <- patient_meta[[id_col]][patient_meta[[order_by_col]] == grp]
    mat_sub <- count_matrix[, patients_in_grp, drop = FALSE]
    idx <- oncoprint_column_order(mat_sub, row_order)
    ord <- c(ord, patients_in_grp[idx])
  }
  ord
}

# SOURCE: RTI_funcs.R (2026-04-07)
# Row ordering for a drug/alteration matrix: order by drug class first, then
# by descending total sample count, then by descending number of mutated samples.
make_row_order = function(count_matrix, drug_classs) {
  n_mut = rowSums(count_matrix > 0)
  order(drug_classs, -rowSums(count_matrix), -n_mut)
}


# --- SECTION 6: DRUG ANNOTATION ----------------------------------------------

# SOURCE: RTI_funcs.R (2026-04-07)
# Map a vector of drug names to canonical synonyms using a synonym map (data.table).
# synmap must have columns: alias_clean (lowercase alias) and canonical_name.
map_drugs <- function(drug_vec, synmap) {
  query <- data.table::data.table(drug_query = drug_vec)
  query[, drug_query_clean := tolower(drug_query)]
  out <- merge(query, synmap,
               by.x = "drug_query_clean",
               by.y = "alias_clean",
               all.x = TRUE)
  return(out)
}

# SOURCE: RTI_funcs.R (2026-04-07)
# From a vector of pathway strings, return the single most specific pathway:
# removes "Other"/"Unclassified" if alternatives exist, then picks the longest.
pick_specific_pathway <- function(paths) {
  cleaned <- paths[!grepl("Other|Unclassified", paths)]

  if (length(cleaned) == 0) cleaned <- paths

  if (length(cleaned) > 1) {
    cleaned <- cleaned[which.max(nchar(cleaned))]
  }

  return(cleaned[1])
}
