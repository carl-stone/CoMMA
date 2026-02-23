#' Annotate methylation sites relative to nearby transcription start sites
#'
#' Matches methylation sites to transcription units in `meta_df` and computes
#' relative position from TSS on the sense (`RelPos_+`) and antisense
#' (`RelPos_-`) strands within a user-specified window.
#'
#' ## Input schema
#' * `methyl_df`: data.frame containing a genomic coordinate column specified by
#'   `location`.
#' * `meta_df`: data.frame with columns `Type`, `Strand`, `Left`, and `Right`.
#'   Only rows where `Type == "Transcription-Units"` are used.
#'
#' ## Return schema
#' * If `long = TRUE` (default): returns long-format rows with all original
#'   `methyl_df` columns plus stable columns `TSS_strand` and `RelPos`, one
#'   row per site-by-TSS overlap. Sites with no overlap are omitted; if no sites
#'   overlap at all, a zero-row data.frame is returned with these columns.
#' * If `long = FALSE`: returns the original rows plus zero or more
#'   `RelPos_+<n>` / `RelPos_-<n>` columns and optional `NoTSS = "X"` marker
#'   for rows without nearby TSS hits.
#'
#' @param methyl_df Data frame of methylation sites.
#' @param meta_df Data frame of genomic features. Required columns:
#'   `Type`, `Strand`, `Left`, `Right`.
#' @param location Character scalar naming a column in `methyl_df` containing
#'   integer positions.
#' @param size Integer scalar window around TSS coordinates to search.
#' @param long Logical; return long-format output when `TRUE`.
#'
#' @return A data.frame with TSS-relative methylation annotations.
#' @export
#'
#' @examples
#' methyl_df <- data.frame(Position = c(100L, 250L))
#' meta_df <- data.frame(
#'   Type = c("Transcription-Units", "Transcription-Units"),
#'   Strand = c("+", "-"),
#'   Left = c(90L, 200L),
#'   Right = c(130L, 260L)
#' )
#' annotateTSS(methyl_df, meta_df, location = "Position", size = 25)
annotateTSS <-
  function(methyl_df, meta_df, location, size, long = TRUE) {
    legacy_soft_deprecate("annotateTSS", "annotate_sites_with_features")
    required_meta_cols <- c("Type", "Strand", "Left", "Right")
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
    if (!is.numeric(size) || length(size) != 1L || is.na(size) || size < 0) {
      stop("`size` must be a single non-negative number.", call. = FALSE)
    }
    if (!is.logical(long) || length(long) != 1L || is.na(long)) {
      stop("`long` must be a single TRUE/FALSE value.", call. = FALSE)
    }

    meta_df <- dplyr::filter(meta_df, Type == "Transcription-Units")
    for (position in methyl_df[[location]]) {
      if (nrow(meta_df[(
        meta_df$Strand == '+' &
        meta_df$Left - size <= position &
        meta_df$Left + size >= position
      ) |
      (
        meta_df$Strand == '-' &
        meta_df$Right - size <= position &
        meta_df$Right + size >= position
      ), ]) == 0) {
        methyl_df[methyl_df[[location]] == position, 'NoTSS'] <- 'X'
        next
      }
      SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
                                       meta_df$Left - size <= position &
                                       meta_df$Left + size >= position, ]
      AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
                                           meta_df$Right - size <= position &
                                           meta_df$Right + size >= position, ]
      for (i in 1:nrow(SenseTU_at_position)) {
        if (nrow(SenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_+', i)] <-
          position - SenseTU_at_position[i, 'Left']
      }
      for (i in 1:nrow(AntisenseTU_at_position)) {
        if (nrow(AntisenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_-', i)] <-
          AntisenseTU_at_position[i, 'Right'] - position
      }
    }
    if (long) {
      relpos_cols <- grep("^RelPos_", names(methyl_df), value = TRUE)
      if (length(relpos_cols) == 0L) {
        methyl_df <- methyl_df[0, , drop = FALSE]
        methyl_df$TSS_strand <- character()
        methyl_df$RelPos <- numeric()
        return(methyl_df)
      }
      methyl_df <- tidyr::pivot_longer(
        methyl_df,
          cols = dplyr::all_of(relpos_cols),
          names_to = 'TSS_strand',
          names_pattern = 'RelPos_(.)[0-9]*',
          values_to = 'RelPos'
        )
      distinct_keys <- paste(methyl_df[[location]], methyl_df$RelPos, sep = "\r")
      methyl_df <- methyl_df[!duplicated(distinct_keys), , drop = FALSE]
      methyl_df <- dplyr::filter(methyl_df, !is.na(RelPos))
      return(methyl_df)
    } else {
      return(methyl_df)
    }
  }
