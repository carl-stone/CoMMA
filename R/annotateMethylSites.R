#' Annotate methylation sites by overlap with genomic features
#'
#' For each methylation site, identify all rows in `meta_df` whose `Left..Right`
#' interval contains that site coordinate. Feature hits are written into output
#' columns named by feature `Type` values.
#'
#' This legacy helper now uses CoMMA's canonical annotation workflow internally:
#' [validate_site_table()] + [normalize_annotation_table()] +
#' [annotate_sites_with_features()]. Output schema and markers are preserved for
#' backward compatibility.
#'
#' ## Input schema
#' * `methyl_df`: data.frame containing a numeric coordinate column specified by
#'   `location`.
#' * `meta_df`: data.frame containing required columns:
#'   `Type` (feature class), `Site` (feature identifier), `Left`, and `Right`.
#'
#' ## Return schema
#' Returns `methyl_df` with additional annotation columns:
#' * one column per matched feature `Type`, populated with `Site` labels
#' * stable marker column `No_Feature = "1"` for rows with no overlaps
#'
#' @param methyl_df Data frame of methylation sites.
#' @param meta_df Data frame of genomic features. Required columns:
#'   `Type`, `Site`, `Left`, `Right`.
#' @param location Character scalar naming a column in `methyl_df` containing
#'   genomic coordinates.
#'
#' @return A data.frame with feature-annotation columns appended.
#' @export
#'
#' @name annotateMethylSites
#'
#' @examples
#' methyl_df <- data.frame(Position = c(100L, 300L))
#' meta_df <- data.frame(
#'   Type = c("Gene", "Promoter"),
#'   Site = c("geneA", "proA"),
#'   Left = c(90L, 95L),
#'   Right = c(200L, 110L)
#' )
#' annotateMethylSites(methyl_df, meta_df, location = "Position")
annotateMethylSites <- function(methyl_df, meta_df, location) {
  required_meta_cols <- c("Type", "Site", "Left", "Right")
  missing_meta_cols <- setdiff(required_meta_cols, names(meta_df))

  if (!is.character(location) || length(location) != 1L || is.na(location)) {
    stop("`location` must be a single, non-missing column name.", call. = FALSE)
  }
  if (!location %in% names(methyl_df)) {
    stop(
      sprintf("`methyl_df` is missing required location column `%s`.", location),
      call. = FALSE
    )
  }
  if (length(missing_meta_cols) > 0L) {
    stop(
      sprintf(
        "`meta_df` is missing required columns: %s.",
        paste(missing_meta_cols, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  site_table <- data.frame(
    seqname = "legacy_seqname",
    pos = as.integer(methyl_df[[location]]),
    strand = "+",
    beta = 0,
    coverage = 1L,
    sample_id = "legacy_sample",
    group = "legacy_group",
    legacy_row_id = seq_len(nrow(methyl_df)),
    stringsAsFactors = FALSE
  )
  site_table <- validate_site_table(site_table)

  feature_table <- normalize_annotation_table(
    annotation_df = data.frame(meta_df, seqname = "legacy_seqname", stringsAsFactors = FALSE),
    feature_type_col = "Type",
    feature_id_col = "Site",
    seqname_col = "seqname",
    start_col = "Left",
    end_col = "Right",
    strand_col = NULL
  )

  annotated <- annotate_sites_with_features(
    site_table = site_table,
    feature_table = feature_table,
    site_pos_col = "pos",
    site_seqname_col = "seqname",
    include_unannotated = TRUE
  )

  out <- methyl_df
  for (row_i in seq_len(nrow(methyl_df))) {
    row_hits <- annotated[annotated$legacy_row_id == row_i, , drop = FALSE]
    if (nrow(row_hits) == 0 || all(is.na(row_hits$feature_type))) {
      out[row_i, "No_Feature"] <- "1"
      next
    }

    for (hit_i in seq_len(nrow(row_hits))) {
      feature_type <- row_hits$feature_type[hit_i]
      feature_id <- row_hits$feature_id[hit_i]
      if (is.na(feature_type) || is.na(feature_id)) {
        next
      }
      out[row_i, feature_type] <- feature_id
    }
  }

  out
}
