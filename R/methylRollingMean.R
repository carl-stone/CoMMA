#' Calculate rolling mean of methylation levels
#'
#' Calculates the rolling mean of methylation levels for a given window size. The function takes in a data frame with columns
#' for genomic position and methylation level, and returns a data frame with columns for position and the rolling mean
#' methylation level for that position.
#'
#' @param df A data frame with columns for genomic position and methylation level.
#' @param position_col The name of the column in \code{df} that contains the genomic position data.
#' @param methyl_col The name of the column in \code{df} that contains the methylation level data.
#' @param w_size The window size for calculating the rolling mean. Default value is 1000.
#' @param genome_size The size of the genome. Default value is 4641652.
#' @param verbose If \code{TRUE}, progress updates will be printed. Default value is \code{FALSE}.
#'
#' @return A data frame with columns for genomic position and the rolling mean methylation level for that position.
#'
#' @examples
#' df <- data.frame(position = c(1, 2, 3, 4, 5), methyl = c(0.1, 0.2, 0.3, 0.4, 0.5))
#' result <- methylRollingMean(df, position_col = "position", methyl_col = "methyl", w_size = 2)
#'
#' @import dplyr
methylRollingMean <- function(df, position_col, methyl_col, w_size=1000, genome_size=4641652, verbose=FALSE) {
  tstart <- Sys.time()
  # Take the first w_size of sites from the beginning of the chromosome and add them to the end so it wraps around
  df <- dplyr::bind_cols(position = df[[position_col]], methyl = df[[methyl_col]])
  nsites <- nrow(df)
  beginning_sites <- df %>%
    dplyr::filter(position <= w_size) %>%
    dplyr::mutate(position = position + genome_size)
  df <- dplyr::bind_rows(df, beginning_sites)
  df <- df[order(df[, 1]),]  # sort df by position
  df <- data.matrix(df)
  # TODO sort df by position
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

  out_df <- as_tibble(out_df)
  if (verbose) print(Sys.time()-tstart)
  return(out_df)
}
