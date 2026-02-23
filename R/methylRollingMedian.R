#' Compute circular sliding-window methylation median across a genome
#'
#' Calculates rolling medians of per-position methylation values on a circular
#' bacterial chromosome. In `method = "exact"`, every genomic position from
#' `1..genome_size` is represented. In `method = "fast"`, medians are computed
#' only for input sites.
#'
#' ## Input schema
#' `df` must contain:
#' * `position_col`: positive integer genomic coordinates
#' * `methyl_col`: numeric methylation values
#'
#' ## Return schema
#' * `method = "exact"`: columns `position`, `methyl`, `med_methyl`
#' * `method = "fast"`: columns `position`, `mean_methyl` (historical name;
#'   values are medians)
#'
#' @param df Data frame with methylation positions and values.
#' @param position_col Column name containing genomic positions.
#' @param methyl_col Column name containing methylation values.
#' @param w_size Sliding window size.
#' @param genome_size Genome size for circular wrap-around.
#' @param method Either `"exact"` (default) or `"fast"`.
#' @param verbose Logical; when `TRUE`, prints loop progress in fast mode.
#'
#' @return A data.frame of sliding-window methylation summary values.
#' @export
#'
#' @examples
#' toy <- data.frame(pos = c(1L, 3L, 5L), beta = c(0.1, 0.8, 0.3))
#' methylRollingMedian(toy, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6)
methylRollingMedian <- function(df,
                                position_col,
                                methyl_col,
                                w_size,
                                genome_size = 4641652,
                                method = "exact",
                                verbose = FALSE) {
  df <- .validate_rolling_input(df, position_col, methyl_col, w_size, genome_size)
  if (!is.character(method) || length(method) != 1L || is.na(method) || !(method %in% c("exact", "fast"))) {
    stop("`method` must be one of: 'exact', 'fast'.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
  }

  if (method == "exact") {
    # Create empty df of every genomic position + extra for the end to wrap around
    all_pos <- data.frame(position = seq(1, genome_size + w_size, 1),
                          methyl = NA)
    all_pos[all_pos$position %in% df$position, ]$methyl <- df$methyl
    # Take the first w_size of sites from the beginning of the chromosome and add them to the end
    all_pos[all_pos$position > genome_size, "methyl"] <- all_pos[all_pos$position <= w_size, "methyl"]
    # Calculate median sliding window
    all_pos$med_methyl <- zoo::rollapply(all_pos$methyl, w_size + 1,
                                         median, na.rm = TRUE, partial = TRUE,
                                         align = "center")
    out_df <- all_pos[1:genome_size, ]
  }

  if (method == "fast") {
    nsites <- nrow(df)
    beginning_sites <- df |>
      dplyr::filter(position <= w_size) |>
      dplyr::mutate(position = position + genome_size)
    df <- dplyr::bind_rows(df, beginning_sites) |>
      dplyr::arrange(position)
    df <- data.matrix(df)
    out_df <- tibble::tibble(position = as.numeric(rep(NA, nsites)),
                             mean_methyl = as.double(rep(NA, nsites)))
    out_df <- data.matrix(out_df)
    for (i in 1:nsites) {
      if (verbose && i %% 1000 == 0) {
        print(i)
      }
      out_df[[i, 'position']] <- df[[i, 'position']]
      out_df[[i, 'mean_methyl']] <- median(df[df[, 1] >= df[[i, 1]] & df[, 1] < (df[[i, 1]] + w_size), 2])
    }
    out_df <- tibble::as_tibble(out_df)
  }
  return(out_df)
}

.validate_rolling_input <- function(df, position_col, methyl_col, w_size, genome_size) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.", call. = FALSE)
  }
  if (!is.character(position_col) || length(position_col) != 1L || is.na(position_col) || !nzchar(position_col)) {
    stop("`position_col` must be a single, non-empty column name.", call. = FALSE)
  }
  if (!is.character(methyl_col) || length(methyl_col) != 1L || is.na(methyl_col) || !nzchar(methyl_col)) {
    stop("`methyl_col` must be a single, non-empty column name.", call. = FALSE)
  }
  missing_cols <- setdiff(c(position_col, methyl_col), names(df))
  if (length(missing_cols) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  position <- df[[position_col]]
  methyl <- df[[methyl_col]]
  if (!is.numeric(position) && !is.integer(position)) {
    stop("`position_col` must contain numeric/integer genomic positions.", call. = FALSE)
  }
  if (any(is.na(position)) || any(!is.finite(position)) || any(position <= 0) || any(position %% 1 != 0)) {
    stop("`position_col` must contain positive, finite integer-like values.", call. = FALSE)
  }
  if (!is.numeric(methyl) && !is.integer(methyl)) {
    stop("`methyl_col` must contain numeric methylation values.", call. = FALSE)
  }
  if (any(is.na(methyl)) || any(!is.finite(methyl))) {
    stop("`methyl_col` must contain only finite, non-missing values.", call. = FALSE)
  }
  if (!is.numeric(w_size) || length(w_size) != 1L || is.na(w_size) || !is.finite(w_size) || w_size <= 0) {
    stop("`w_size` must be a single number > 0.", call. = FALSE)
  }
  if (!is.numeric(genome_size) || length(genome_size) != 1L || is.na(genome_size) || !is.finite(genome_size) || genome_size <= 0) {
    stop("`genome_size` must be a single number > 0.", call. = FALSE)
  }

  dplyr::bind_cols(position = as.numeric(position), methyl = as.numeric(methyl))
}
