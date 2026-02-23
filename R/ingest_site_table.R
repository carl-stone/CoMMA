#' Validate CoMMA canonical methylation site table
#'
#' Validate and standardize a data frame into CoMMA's canonical internal schema.
#' Required columns are common metadata columns
#' (`seqname`, `pos`, `strand`, `sample_id`, `group`) plus exactly one
#' measurement mode:
#' * count mode: `n_mod` + `n_total`
#' * fraction mode: `beta` + `coverage`
#'
#' In fraction mode, approximate counts are backfilled as
#' `n_mod = round(beta * coverage)` and `n_total = coverage` to keep a stable
#' canonical schema for downstream methods requiring integer counts.
#' `mod_base` is supported and defaults to `"6mA"` when absent. `motif` is
#' optional and defaults to `"GATC"` when absent.
#'
#' The output always includes canonical columns:
#' `seqname`, `pos`, `strand`, `mod_base`, `motif`, `n_mod`, `n_total`,
#' `sample_id`, `group`, and `beta`.
#'
#' @param site_table A data frame containing per-site methylation calls.
#' @param default_mod_base Default modification label used when `mod_base`
#'   column is missing. Defaults to `"6mA"`.
#' @param default_motif Default motif label used when `motif` column is
#'   missing. Defaults to `"GATC"`.
#'
#' @return A validated data frame with canonical columns and derived/backfilled
#'   measurement fields.
#' @export
#'
#' @examples
#' site_table <- data.frame(
#'   seqname = c("chr1", "chr1"),
#'   pos = c(100L, 200L),
#'   strand = c("+", "-"),
#'   mod_base = c("6mA", "6mA"),
#'   n_mod = c(15, 8),
#'   n_total = c(20, 10),
#'   sample_id = c("S1", "S1"),
#'   group = c("WT", "WT")
#' )
#' validate_site_table(site_table)
validate_site_table <- function(
  site_table,
  default_mod_base = "6mA",
  default_motif = "GATC"
) {
  if (!is.data.frame(site_table)) {
    stop("`site_table` must be a data.frame.", call. = FALSE)
  }

  if (!is.character(default_mod_base) || length(default_mod_base) != 1 || is.na(default_mod_base) || default_mod_base == "") {
    stop("`default_mod_base` must be a single non-empty character value.", call. = FALSE)
  }
  if (!is.character(default_motif) || length(default_motif) != 1 || is.na(default_motif) || default_motif == "") {
    stop("`default_motif` must be a single non-empty character value.", call. = FALSE)
  }

  common_required_cols <- c("seqname", "pos", "strand", "sample_id", "group")

  missing_common_cols <- setdiff(common_required_cols, names(site_table))
  if (length(missing_common_cols) > 0) {
    stop(
      paste0(
        "Missing required columns: ",
        paste(missing_common_cols, collapse = ", "),
        ". Common required columns are: ",
        paste(common_required_cols, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  has_n_mod <- "n_mod" %in% names(site_table)
  has_n_total <- "n_total" %in% names(site_table)
  has_beta <- "beta" %in% names(site_table)
  has_coverage <- "coverage" %in% names(site_table)

  has_count_mode <- has_n_mod && has_n_total
  has_fraction_mode <- has_beta && has_coverage

  partial_count_mode <- xor(has_n_mod, has_n_total)
  partial_fraction_mode <- xor(has_beta, has_coverage)

  if (partial_count_mode || partial_fraction_mode || (has_count_mode && has_fraction_mode) || (!has_count_mode && !has_fraction_mode)) {
    stop(
      paste0(
        "Provide exactly one complete measurement mode: ",
        "count mode (`n_mod` + `n_total`) or fraction mode (`beta` + `coverage`)."
      ),
      call. = FALSE
    )
  }

  out <- site_table

  if (!"mod_base" %in% names(out)) {
    out$mod_base <- rep(default_mod_base, nrow(out))
  }
  if (!"motif" %in% names(out)) {
    out$motif <- rep(default_motif, nrow(out))
  }

  out$seqname <- as.character(out$seqname)
  out$strand <- as.character(out$strand)
  out$mod_base <- as.character(out$mod_base)
  out$motif <- as.character(out$motif)
  out$sample_id <- as.character(out$sample_id)
  out$group <- as.character(out$group)

  if (any(is.na(out$seqname) | out$seqname == "")) {
    stop("`seqname` must not contain missing or empty values.", call. = FALSE)
  }
  if (any(is.na(out$mod_base) | out$mod_base == "")) {
    stop("`mod_base` must not contain missing or empty values.", call. = FALSE)
  }
  if (any(is.na(out$sample_id) | out$sample_id == "")) {
    stop("`sample_id` must not contain missing or empty values.", call. = FALSE)
  }
  if (any(is.na(out$group) | out$group == "")) {
    stop("`group` must not contain missing or empty values.", call. = FALSE)
  }

  if (!is.numeric(out$pos)) {
    stop("`pos` must be numeric/integer genomic positions.", call. = FALSE)
  }
  if (any(!is.na(out$pos) & (out$pos <= 0 | out$pos != floor(out$pos)))) {
    stop("`pos` must contain positive integer positions.", call. = FALSE)
  }

  valid_strands <- c("+", "-")
  if (any(is.na(out$strand) | !(out$strand %in% valid_strands))) {
    stop("`strand` must contain only '+' or '-' values.", call. = FALSE)
  }

  if (has_count_mode) {
    if (!is.numeric(out$n_mod)) {
      stop("`n_mod` must be numeric/integer counts.", call. = FALSE)
    }
    if (!is.numeric(out$n_total)) {
      stop("`n_total` must be numeric/integer counts.", call. = FALSE)
    }

    if (any(!is.na(out$n_mod) & (out$n_mod < 0 | out$n_mod != floor(out$n_mod)))) {
      stop("`n_mod` must contain non-negative integer counts.", call. = FALSE)
    }
    if (any(!is.na(out$n_total) & (out$n_total <= 0 | out$n_total != floor(out$n_total)))) {
      stop("`n_total` must contain positive integer counts.", call. = FALSE)
    }

    if (any(!is.na(out$n_mod) & !is.na(out$n_total) & out$n_mod > out$n_total)) {
      stop("`n_mod` must not be greater than `n_total`.", call. = FALSE)
    }

    missing_counts <- is.na(out$n_mod) | is.na(out$n_total)
    if (any(missing_counts)) {
      warning(
        sprintf(
          "Detected %d rows with missing `n_mod` or `n_total`; setting `beta` to NA for those rows.",
          sum(missing_counts)
        ),
        call. = FALSE
      )
    }

    out$beta <- out$n_mod / out$n_total
    out$beta[missing_counts] <- NA_real_
  } else {
    if (!is.numeric(out$beta)) {
      stop("`beta` must be numeric fractions in [0, 1].", call. = FALSE)
    }
    if (any(!is.na(out$beta) & (out$beta < 0 | out$beta > 1))) {
      stop("`beta` must contain values in [0, 1].", call. = FALSE)
    }

    if (!is.numeric(out$coverage)) {
      stop("`coverage` must be numeric/integer counts.", call. = FALSE)
    }
    if (any(!is.na(out$coverage) & (out$coverage <= 0 | out$coverage != floor(out$coverage)))) {
      stop("`coverage` must contain positive integer counts.", call. = FALSE)
    }

    out$n_total <- out$coverage
    out$n_mod <- as.integer(round(out$beta * out$coverage))
  }

  canonical_cols <- c(
    "seqname", "pos", "strand", "mod_base", "motif",
    "n_mod", "n_total", "sample_id", "group", "beta"
  )
  out <- out[, unique(c(canonical_cols, names(out))), drop = FALSE]

  out
}

#' Convert modkit summarized bedMethyl-like output to canonical site schema
#'
#' Convert modkit summarized output (bedMethyl-like table) into CoMMA's
#' canonical site schema and validate it with [validate_site_table()].
#'
#' Column mapping defaults:
#' * `chrom` -> `seqname`
#' * `start` (0-based) -> `pos` (1-based, `start + 1`)
#' * `strand` -> `strand`
#' * `modified_primary_base` (or `mod_base` if already present) -> `mod_base`
#' * `n_mod` -> `n_mod`
#' * `valid_coverage` (or `n_total`) -> `n_total`
#'
#' @param modkit_df A data frame from modkit summarized output.
#' @param sample_id Sample identifier to assign to every row.
#' @param group Group/condition label to assign to every row.
#'
#' @return A validated canonical site table with derived `beta`.
#' @export
#'
#' @examples
#' modkit_df <- data.frame(
#'   chrom = c("chr1", "chr1"),
#'   start = c(99L, 199L),
#'   strand = c("+", "-"),
#'   modified_primary_base = c("6mA", "6mA"),
#'   n_mod = c(15, 8),
#'   valid_coverage = c(20, 10)
#' )
#' convert_modkit_bedmethyl(modkit_df, sample_id = "S1", group = "WT")
convert_modkit_bedmethyl <- function(modkit_df, sample_id, group) {
  if (!is.data.frame(modkit_df)) {
    stop("`modkit_df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(sample_id) || length(sample_id) != 1 || is.na(sample_id) || sample_id == "") {
    stop("`sample_id` must be a single non-empty character value.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1 || is.na(group) || group == "") {
    stop("`group` must be a single non-empty character value.", call. = FALSE)
  }

  seqname_col <- if ("chrom" %in% names(modkit_df)) "chrom" else "seqname"
  pos_col <- if ("start" %in% names(modkit_df)) "start" else "pos"
  mod_base_col <- if ("modified_primary_base" %in% names(modkit_df)) {
    "modified_primary_base"
  } else {
    "mod_base"
  }
  n_total_col <- if ("valid_coverage" %in% names(modkit_df)) {
    "valid_coverage"
  } else {
    "n_total"
  }

  required_source <- c(seqname_col, pos_col, "strand", mod_base_col, "n_mod", n_total_col)
  missing_source <- setdiff(required_source, names(modkit_df))
  if (length(missing_source) > 0) {
    stop(
      paste0(
        "modkit input is missing required columns: ",
        paste(missing_source, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  out <- data.frame(
    seqname = modkit_df[[seqname_col]],
    pos = modkit_df[[pos_col]],
    strand = modkit_df[["strand"]],
    mod_base = modkit_df[[mod_base_col]],
    n_mod = modkit_df[["n_mod"]],
    n_total = modkit_df[[n_total_col]],
    sample_id = rep(sample_id, nrow(modkit_df)),
    group = rep(group, nrow(modkit_df)),
    stringsAsFactors = FALSE
  )

  if (identical(pos_col, "start")) {
    out$pos <- out$pos + 1
  }

  validate_site_table(out)
}

#' Convert Megalodon-like per-site call table to canonical site schema
#'
#' Convert a Megalodon-like TSV/table into CoMMA's canonical site schema,
#' then validate with [validate_site_table()]. This converter is designed to be
#' pragmatic across minor column-name differences.
#'
#' Default column matching (first existing name is used):
#' * sequence: `chromosome` or `contig` or `seqname`
#' * position: `position` or `pos`
#' * strand: `strand`
#' * modified base: `mod_base` or `modified_base`
#' * modified count: `n_mod` or `mod_count`
#' * total count: `n_total` or `coverage` or `read_depth`
#' * motif (optional): `motif`
#'
#' @param megalodon_df A data frame with Megalodon-like per-site rows.
#' @param sample_id Sample identifier to assign to every row.
#' @param group Group/condition label to assign to every row.
#' @param position_is_zero_based Logical; if `TRUE`, `pos` is incremented by 1
#'   during conversion. Defaults to `FALSE`.
#'
#' @return A validated canonical site table with derived `beta`.
#' @export
#'
#' @examples
#' megalodon_df <- data.frame(
#'   chromosome = c("chr1", "chr1"),
#'   position = c(99L, 199L),
#'   strand = c("+", "-"),
#'   modified_base = c("6mA", "6mA"),
#'   mod_count = c(8L, 6L),
#'   read_depth = c(10L, 12L)
#' )
#' convert_megalodon_tsv(megalodon_df, sample_id = "S1", group = "WT", position_is_zero_based = TRUE)
convert_megalodon_tsv <- function(megalodon_df, sample_id, group, position_is_zero_based = FALSE) {
  if (!is.data.frame(megalodon_df)) {
    stop("`megalodon_df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(sample_id) || length(sample_id) != 1 || is.na(sample_id) || sample_id == "") {
    stop("`sample_id` must be a single non-empty character value.", call. = FALSE)
  }
  if (!is.character(group) || length(group) != 1 || is.na(group) || group == "") {
    stop("`group` must be a single non-empty character value.", call. = FALSE)
  }
  if (!is.logical(position_is_zero_based) || length(position_is_zero_based) != 1 || is.na(position_is_zero_based)) {
    stop("`position_is_zero_based` must be TRUE or FALSE.", call. = FALSE)
  }

  pick_col <- function(candidates, data, label) {
    found <- candidates[candidates %in% names(data)]
    if (length(found) == 0) {
      stop(
        sprintf("Megalodon-like input is missing a %s column. Accepted names: %s.", label, paste(candidates, collapse = ", ")),
        call. = FALSE
      )
    }
    found[[1]]
  }

  seqname_col <- pick_col(c("chromosome", "contig", "seqname"), megalodon_df, "sequence")
  pos_col <- pick_col(c("position", "pos"), megalodon_df, "position")
  mod_base_col <- pick_col(c("mod_base", "modified_base"), megalodon_df, "modified base")
  n_mod_col <- pick_col(c("n_mod", "mod_count"), megalodon_df, "modified-count")
  n_total_col <- pick_col(c("n_total", "coverage", "read_depth"), megalodon_df, "total-count")

  if (!"strand" %in% names(megalodon_df)) {
    stop("Megalodon-like input is missing required column: strand.", call. = FALSE)
  }

  out <- data.frame(
    seqname = megalodon_df[[seqname_col]],
    pos = megalodon_df[[pos_col]],
    strand = megalodon_df[["strand"]],
    mod_base = megalodon_df[[mod_base_col]],
    n_mod = megalodon_df[[n_mod_col]],
    n_total = megalodon_df[[n_total_col]],
    motif = if ("motif" %in% names(megalodon_df)) megalodon_df[["motif"]] else NA_character_,
    sample_id = rep(sample_id, nrow(megalodon_df)),
    group = rep(group, nrow(megalodon_df)),
    stringsAsFactors = FALSE
  )

  if (isTRUE(position_is_zero_based)) {
    out$pos <- out$pos + 1
  }

  validate_site_table(out)
}

#' Convert generic long beta+coverage methylation table to canonical site schema
#'
#' Convert generic long-format methylation data with beta and coverage columns
#' into CoMMA's canonical site schema. This is useful when callers already
#' report methylation fractions instead of modified/unmodified counts.
#'
#' Required source columns are `seqname`, `pos`, `strand`, `beta`, `coverage`,
#' `sample_id`, and `group`. Optional columns are `mod_base` and `motif`.
#'
#' @param long_df A generic long-format methylation table.
#' @param position_is_zero_based Logical; if `TRUE`, `pos` is incremented by 1
#'   during conversion. Defaults to `FALSE`.
#' @param default_mod_base Default `mod_base` when absent. Defaults to `"6mA"`.
#' @param default_motif Default `motif` when absent. Defaults to `"GATC"`.
#'
#' @return A validated canonical site table with backfilled count columns.
#' @export
#'
#' @examples
#' long_df <- data.frame(
#'   seqname = c("chr1", "chr1"),
#'   pos = c(100L, 200L),
#'   strand = c("+", "-"),
#'   beta = c(0.7, 0.1),
#'   coverage = c(10L, 20L),
#'   sample_id = c("S1", "S2"),
#'   group = c("WT", "mut"),
#'   mod_base = c("5mC", "5mC")
#' )
#' convert_beta_coverage_long(long_df)
convert_beta_coverage_long <- function(
  long_df,
  position_is_zero_based = FALSE,
  default_mod_base = "6mA",
  default_motif = "GATC"
) {
  if (!is.data.frame(long_df)) {
    stop("`long_df` must be a data.frame.", call. = FALSE)
  }
  if (!is.logical(position_is_zero_based) || length(position_is_zero_based) != 1 || is.na(position_is_zero_based)) {
    stop("`position_is_zero_based` must be TRUE or FALSE.", call. = FALSE)
  }

  required_source <- c("seqname", "pos", "strand", "beta", "coverage", "sample_id", "group")
  missing_source <- setdiff(required_source, names(long_df))
  if (length(missing_source) > 0) {
    stop(
      paste0(
        "generic beta+coverage input is missing required columns: ",
        paste(missing_source, collapse = ", "),
        "."
      ),
      call. = FALSE
    )
  }

  out <- data.frame(
    seqname = long_df[["seqname"]],
    pos = long_df[["pos"]],
    strand = long_df[["strand"]],
    beta = long_df[["beta"]],
    coverage = long_df[["coverage"]],
    mod_base = if ("mod_base" %in% names(long_df)) long_df[["mod_base"]] else rep(default_mod_base, nrow(long_df)),
    motif = if ("motif" %in% names(long_df)) long_df[["motif"]] else rep(default_motif, nrow(long_df)),
    sample_id = long_df[["sample_id"]],
    group = long_df[["group"]],
    stringsAsFactors = FALSE
  )

  if (isTRUE(position_is_zero_based)) {
    out$pos <- out$pos + 1
  }

  validate_site_table(out, default_mod_base = default_mod_base, default_motif = default_motif)
}
