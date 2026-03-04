#' Run differential methylation analysis with methylKit
#'
#' High-level differential methylation wrapper that enforces CoMMA's manuscript
#' defaults while keeping all key thresholds configurable.
#'
#' The function expects CoMMA's canonical site-table schema (see
#' [validate_site_table()]) and performs:
#' 
#' 1. Optional removal of mutated sites before model fitting.
#' 2. Coverage filtering (`n_total >= coverage_min`, default 10).
#' 3. Coverage normalization by median coverage.
#' 4. Differential methylation with methylKit logistic regression and
#'    SLIM multiple-testing correction.
#' 5. Differential-site flagging at default thresholds `q < 0.05` and
#'    `abs(percent_diff) > 10`.
#'
#' @param site_table Canonical CoMMA site table with columns at minimum:
#'   `seqname`, `pos`, `strand`, `n_mod`, `n_total`, `sample_id`, and `group`.
#'   `mod_base` and `motif` are supported; missing values are defaulted to
#'   `"6mA"` and `"GATC"` respectively.
#' @param mod_base Character scalar specifying which modification type to
#'   analyze. Default `"6mA"` preserves historical behavior. Set to `NULL` to
#'   skip modification-type filtering.
#' @param motif Character scalar specifying which motif to analyze. Default
#'   `"GATC"` preserves historical behavior. Set to `NULL` to skip motif
#'   filtering.
#' @param remove_mutated_sites Logical; whether to remove mutated sites before
#'   fitting the model. Default `TRUE` to preserve manuscript behavior.
#' @param mutated_sites Optional specification of mutated loci used when
#'   `remove_mutated_sites = TRUE`. Supported forms:
#'   * `NULL` (default): if `site_table` has logical column `mutated` or
#'     `site_mutated`, loci are inferred from rows marked `TRUE` and all rows
#'     for those loci are removed (locus identity uses `seqname` + `pos`, and
#'     includes `strand` when present).
#'   * Numeric/integer vector: interpreted as genomic positions (`pos`) to
#'     remove.
#'   * Data frame containing `seqname` and `pos` columns (optional `strand`).
#' @param coverage_min Minimum total coverage required to keep a site.
#'   Default is `10`.
#' @param qvalue_threshold Q-value significance threshold. Default is `0.05`.
#' @param percent_diff_threshold Absolute percent methylation-difference
#'   threshold. Default is `10`.
#' @param reference_group Character scalar naming the group that should be
#'   assigned `treatment = 0` (the reference/denominator in methylKit). The
#'   sign of `percent_diff` in the result is `non-reference minus reference`.
#'   Default `NULL` preserves the historical behavior of assigning treatment
#'   by alphabetical order of group names.
#' @param destrand Logical; passed to [methylKit::unite()]. Default `FALSE`.
#' @param context_label Character scalar metadata label stored in methylKit
#'   objects via `methylRaw(..., context = context_label)`. Default `"6mA_GATC"`
#'   reflects the canonical bacterial modification/motif context. Pass `"CpG"`
#'   explicitly to restore the previous default if needed.
#'
#' @return A named list with stable components:
#' * `result_table`: tidy differential methylation table containing
#'   `seqname`, `pos`, `strand`, `pvalue`, `qvalue`, `percent_diff`,
#'   `significant`.
#' * `plot_data`: list with plotting-ready components:
#'   * `volcano_data`: `result_table` plus `neg_log10_qvalue`.
#'   * `beta_by_sample`: per-site beta values with sample/group metadata.
#' * `methylkit_objects`: list containing intermediate methylKit objects
#'   (`methyl_raw`, `normalized`, `united`, `diff`).
#'
#' @export
run_differential_methylation <- function(
  site_table,
  mod_base = "6mA",
  motif = "GATC",
  reference_group = NULL,
  remove_mutated_sites = TRUE,
  mutated_sites = NULL,
  coverage_min = 10,
  qvalue_threshold = 0.05,
  percent_diff_threshold = 10,
  destrand = FALSE,
  context_label = "6mA_GATC"
) {
  if (!is.null(mod_base) && (!is.character(mod_base) || length(mod_base) != 1 || is.na(mod_base) || mod_base == "")) {
    stop("`mod_base` must be NULL or a single non-empty character value.", call. = FALSE)
  }
  if (!is.null(motif) && (!is.character(motif) || length(motif) != 1 || is.na(motif) || motif == "")) {
    stop("`motif` must be NULL or a single non-empty character value.", call. = FALSE)
  }
  if (!is.null(reference_group) && (!is.character(reference_group) || length(reference_group) != 1 || is.na(reference_group) || reference_group == "")) {
    stop("`reference_group` must be NULL or a single non-empty character value.", call. = FALSE)
  }
  if (!is.logical(remove_mutated_sites) || length(remove_mutated_sites) != 1 || is.na(remove_mutated_sites)) {
    stop("`remove_mutated_sites` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  if (!is.numeric(coverage_min) || length(coverage_min) != 1 || is.na(coverage_min) || coverage_min < 0) {
    stop("`coverage_min` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(qvalue_threshold) || length(qvalue_threshold) != 1 || is.na(qvalue_threshold) || qvalue_threshold < 0 || qvalue_threshold > 1) {
    stop("`qvalue_threshold` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(percent_diff_threshold) || length(percent_diff_threshold) != 1 || is.na(percent_diff_threshold) || percent_diff_threshold < 0) {
    stop("`percent_diff_threshold` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.character(context_label) || length(context_label) != 1 || is.na(context_label) || context_label == "") {
    stop("`context_label` must be a single non-empty character value.", call. = FALSE)
  }

  canonical <- validate_site_table(site_table)

  if (!is.null(mod_base)) {
    canonical <- canonical[canonical$mod_base == mod_base, , drop = FALSE]
  }
  if (!is.null(motif) && "motif" %in% names(canonical)) {
    canonical <- canonical[canonical$motif == motif, , drop = FALSE]
  }

  if (nrow(canonical) == 0) {
    stop("No rows remain after `mod_base`/`motif` filtering. Relax filters or inspect inputs.", call. = FALSE)
  }

  if (remove_mutated_sites) {
    canonical <- .remove_mutated_sites(canonical, mutated_sites)
  }

  filtered <- canonical[canonical$n_total >= coverage_min, , drop = FALSE]
  if (nrow(filtered) == 0) {
    stop("No rows remain after coverage filtering. Lower `coverage_min` or provide deeper data.", call. = FALSE)
  }

  # Require at least two groups for differential testing.
  group_by_sample <- unique(filtered[c("sample_id", "group")])
  if (nrow(group_by_sample) < 2 || length(unique(group_by_sample$group)) < 2) {
    stop("Differential methylation requires at least two groups.", call. = FALSE)
  }

  if (!is.null(reference_group) && !reference_group %in% group_by_sample$group) {
    stop(sprintf(
      "`reference_group` '%s' not found in observed groups: %s.",
      reference_group, paste(sort(unique(group_by_sample$group)), collapse = ", ")
    ), call. = FALSE)
  }

  mk_raw <- .site_table_to_methylkit(filtered, context_label = context_label, reference_group = reference_group)
  mk_norm <- .normalize_methylkit_coverage(mk_raw)
  mk_united <- unite(mk_norm, destrand = destrand)
  mk_diff <- .fit_methylkit_diff(mk_united)

  diff_df <- .methylkit_diff_to_df(mk_diff)
  tidy <- .build_tidy_result(diff_df, qvalue_threshold, percent_diff_threshold)

  beta_by_sample <- filtered[, c("seqname", "pos", "strand", "mod_base", "sample_id", "group", "beta"), drop = FALSE]
  volcano_data <- tidy
  volcano_data$neg_log10_qvalue <- -log10(pmax(volcano_data$qvalue, .Machine$double.xmin))

  list(
    result_table = tidy,
    plot_data = list(
      volcano_data = volcano_data,
      beta_by_sample = beta_by_sample
    ),
    methylkit_objects = list(
      methyl_raw = mk_raw,
      normalized = mk_norm,
      united = mk_united,
      diff = mk_diff
    )
  )
}

.remove_mutated_sites <- function(site_table, mutated_sites) {
  if (is.null(mutated_sites)) {
    mutated_col <- intersect(c("site_mutated", "mutated"), names(site_table))
    if (length(mutated_col) == 1 && is.logical(site_table[[mutated_col]])) {
      mutated_rows <- site_table[[mutated_col]] %in% TRUE
      if (!any(mutated_rows)) {
        return(site_table)
      }

      has_strand <- "strand" %in% names(site_table)
      key_site <- if (has_strand) {
        paste(site_table$seqname, site_table$pos, site_table$strand, sep = "::")
      } else {
        paste(site_table$seqname, site_table$pos, sep = "::")
      }
      key_mut <- unique(key_site[mutated_rows])
      return(site_table[!(key_site %in% key_mut), , drop = FALSE])
    }
    return(site_table)
  }

  if (is.numeric(mutated_sites) || is.integer(mutated_sites)) {
    return(site_table[!(site_table$pos %in% mutated_sites), , drop = FALSE])
  }

  if (is.data.frame(mutated_sites)) {
    if (!all(c("seqname", "pos") %in% names(mutated_sites))) {
      stop("`mutated_sites` data.frame must include `seqname` and `pos`.", call. = FALSE)
    }
    has_strand <- "strand" %in% names(mutated_sites)
    key_site <- if (has_strand) {
      paste(site_table$seqname, site_table$pos, site_table$strand, sep = "::")
    } else {
      paste(site_table$seqname, site_table$pos, sep = "::")
    }
    key_mut <- if (has_strand) {
      paste(mutated_sites$seqname, mutated_sites$pos, mutated_sites$strand, sep = "::")
    } else {
      paste(mutated_sites$seqname, mutated_sites$pos, sep = "::")
    }
    return(site_table[!(key_site %in% key_mut), , drop = FALSE])
  }

  stop("`mutated_sites` must be NULL, numeric positions, or a data.frame of loci.", call. = FALSE)
}

.site_table_to_methylkit <- function(site_table, context_label = "6mA_GATC", reference_group = NULL) {
  sample_map <- unique(site_table[c("sample_id", "group")])
  sample_ids <- as.character(sample_map$sample_id)
  if (!is.null(reference_group)) {
    group_levels <- c(reference_group, setdiff(as.character(sample_map$group), reference_group))
    treatment <- as.integer(factor(sample_map$group, levels = group_levels)) - 1L
  } else {
    treatment <- as.integer(as.factor(sample_map$group)) - 1L
  }

  methyl_list <- vector("list", length(sample_ids))
  for (i in seq_along(sample_ids)) {
    sid <- sample_ids[[i]]
    sample_df <- site_table[site_table$sample_id == sid, , drop = FALSE]
    sample_df <- sample_df[order(sample_df$seqname, sample_df$pos), , drop = FALSE]

    sample_df$numTs <- sample_df$n_total - sample_df$n_mod
    sample_df$chrBase <- paste(sample_df$seqname, sample_df$pos, sep = ".")
    sample_df$coverage <- sample_df$n_total
    sample_df$freqC <- (sample_df$n_mod / pmax(sample_df$n_total, 1)) * 100
    sample_df$freqT <- 100 - sample_df$freqC

    mk_df <- sample_df[, c("chrBase", "seqname", "pos", "strand", "coverage", "n_mod", "numTs", "freqC", "freqT")]
    names(mk_df) <- c("chrBase", "chr", "start", "strand", "coverage", "numCs", "numTs", "freqC", "freqT")

    methyl_list[[i]] <- methylRaw(
      mk_df,
      sample.id = sid,
      assembly = "unknown",
      context = context_label,
      resolution = "base"
    )
  }

  methylRawList(methyl_list, treatment = treatment)
}

.normalize_methylkit_coverage <- function(mk_raw) {
  normalizeCoverage(mk_raw, method = "median")
}

.fit_methylkit_diff <- function(mk_united) {
  calculateDiffMeth(
    mk_united,
    overdispersion = "MN",
    test = "Chisq",
    mc.cores = 1,
    correction = "SLIM"
  )
}

.methylkit_diff_to_df <- function(mk_diff) {
  if ("getData" %in% getNamespaceExports("methylKit")) {
    as.data.frame(getData(mk_diff))
  } else {
    as.data.frame(mk_diff)
  }
}

.build_tidy_result <- function(diff_df, qvalue_threshold, percent_diff_threshold) {
  required <- c("chr", "start", "strand", "pvalue", "qvalue", "meth.diff")
  missing <- setdiff(required, names(diff_df))
  if (length(missing) > 0) {
    stop(
      sprintf("Differential results missing expected columns: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  out <- diff_df[, c("chr", "start", "strand", "pvalue", "qvalue", "meth.diff"), drop = FALSE]
  names(out) <- c("seqname", "pos", "strand", "pvalue", "qvalue", "percent_diff")
  out$significant <- !is.na(out$qvalue) & !is.na(out$percent_diff) &
    out$qvalue < qvalue_threshold & abs(out$percent_diff) > percent_diff_threshold
  out
}
