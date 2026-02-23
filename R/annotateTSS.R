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
#' * If `long = TRUE` (default): returns long-format rows with columns
#'   `TSS_strand` and `RelPos`, one row per site-by-TSS overlap.
#' * If `long = FALSE`: returns the original rows plus zero or more
#'   `RelPos_+<n>` / `RelPos_-<n>` columns and optional `NoTSS` marker.
#'
#' @param methyl_df Data frame of methylation sites.
#' @param meta_df Data frame of genomic features.
#' @param location Column name in `methyl_df` containing integer positions.
#' @param size Integer window around TSS coordinates to search.
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
    meta_df <- meta_df %>%
      dplyr::filter(Type == 'Transcription-Units')
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
      methyl_df <-
        methyl_df %>% tidyr::pivot_longer(
          cols = dplyr::starts_with('RelPos'),
          names_to = 'TSS_strand',
          names_pattern = 'RelPos_(.)[0-9]*',
          values_to = 'RelPos'
        )
      methyl_df <-
        methyl_df %>% dplyr::distinct(Position, RelPos, .keep_all = TRUE)
      methyl_df <- methyl_df %>% dplyr::filter(!is.na(RelPos))
      return(methyl_df)
    } else {
      return(methyl_df)
    }
  }
