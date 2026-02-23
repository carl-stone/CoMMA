methylRollingMean <- function(df,
                              position_col,
                              methyl_col,
                              w_size = 1000,
                              genome_size = 4641652,
                              method = "fast",
                              verbose = FALSE) {
  legacy_soft_deprecate("methylRollingMean", "methylRollingMedian")
  if (!is.character(method) || length(method) != 1L || is.na(method) || method != "fast") {
    stop("`method` must be 'fast'.", call. = FALSE)
  }
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("`verbose` must be a single TRUE/FALSE value.", call. = FALSE)
  }
  tstart <- Sys.time()
  # Take the first w_size of sites from the beginning of the chromosome and add them to the end so it wraps around
  df <- .validate_rolling_input(df, position_col, methyl_col, w_size, genome_size)
  nsites <- nrow(df)
  beginning_sites <- df |>
    dplyr::filter(position <= w_size) |>
    dplyr::mutate(position = position + genome_size)
  df <- dplyr::bind_rows(df, beginning_sites) |>
    dplyr::arrange(position)
  df <- data.matrix(df)
  out_df <- data.frame(position = as.numeric(rep(NA, nsites)),
                   mean_methyl = as.double(rep(NA, nsites)))
  out_df <- data.matrix(out_df)

  for (i in 1:nsites) {
    if(verbose && i %% 1000 == 0) {
      print(i)
    }
    out_df[[i,'position']] <- df[[i,'position']]
    out_df[[i,'mean_methyl']] <- mean(df[df[,1] >= df[[i,1]] & df[,1] < (df[[i,1]]+w_size), 2])
  }

  out_df <- tibble::as_tibble(out_df)
  if (verbose) print(Sys.time()-tstart)
  return(out_df)
}
