#' Annotate methylation sites by overlap with genomic features
#'
#' For each methylation site, identify all rows in `meta_df` whose `Left..Right`
#' interval contains that site coordinate. Feature hits are written into output
#' columns named by feature `Type` values.
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

  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[meta_df$Left <= position &
                     meta_df$Right >= position, ]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'No_Feature'] <- '1'
      next
    }
    sites_at_position <- meta_df[meta_df$Left <= position &
                                   meta_df$Right >= position, ]
    for (i in 1:nrow(sites_at_position)) {
      methyl_df[methyl_df[[location]] == position,
                toString(sites_at_position[i, 'Type'])] <-
        sites_at_position[i, 'Site']
    }
  }
  return(methyl_df)
}
