#' Title Detrend Coverage of Methylation Sites by Position
#'
#' @param coverage a numeric vector of coverage values for each site.
#' @param position a numeric vector of positions of each GATC site.
#'
#' @return a numeric vector the same length as coverage of altered coverage
#'  values.
#' @export
#'
#' @examples
#' cov <- detrendCoverage(WT_average$cov, WT_average$Position)
detrendCoverage <- function(coverage, position) {
  # Coverage and position are both numeric vectors
  # Output vector of coverages
  cov_pos <- data.frame(position, coverage)
  colnames(cov_pos) <- c("position", "coverage")
  fit <- MASS::rlm(coverage ~ splines::ns(position, 3), data = cov_pos)
  cov_pos <- cov_pos %>%
    modelr::add_residuals(fit)
  cov_pos$new_cov <- cov_pos$resid + (median(cov_pos$coverage) - median(cov_pos$resid))
  return(cov_pos$new_cov)
}
