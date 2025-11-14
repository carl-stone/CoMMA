#' normalizeMethylation: Function to normalize and rescale methylation data
#'
#' @description The function takes a data frame with columns 'Position' and coverage starting with 'cov' and betas
#' starting with 'methyl' and normalizes the methylation beta values to number of methylated reads.
#' Coverage is normalized by quantile for each sample and transformed to normalized coverage by
#' fitting a linear model to normalized coverage ~ original coverage. The intercept and slope from the
#' linear model are then used to normalize each sample's number of methylated reads.
#'
#' @param df A data frame with columns 'Position' and coverage starting with 'cov' and betas starting
#' with 'methyl'.
#' @param alpha A numeric value to add to the denominator in the formula for methyl_df_normalized to
#' avoid division by zero.
#' @param normalize_position A logical value to determine if position should be normalized.
#' @param rescale A logical value to determine if the normalized methylation beta values should be rescaled.
#' @param plots A logical value to determine if diagnostic plots should be generated.
#' @return A list with three elements:
#' 1. methyl_df_normalized - a matrix with the normalized methylation data.
#' 2. models - a list of linear models fit to normalized coverage ~ original coverage.
#' 3. diagnostics - a list containing diagnostic plots (if plots=TRUE), and diagnostic metrics.
#' If plots = FALSE, diagnostics will only contain diagnostic metrics.
#' @details The function first drops rows with NA values, extracts the columns starting with 'cov' and
#' 'methyl', and converts the methylation beta values to number of methylated reads. Coverage is then
#' normalized by quantile for each sample, and normalized coverage is transformed by fitting a linear model
#' to normalized coverage ~ original coverage, using the intercept and slope from the linear model to normalize
#' each sample's number of methylated reads. If rescale=TRUE, the normalized methylation values are rescaled to
#' have minimum zero and maximum one. If plots = TRUE, the function generates diagnostic plots for each sample
#' to visually check the normalization.
#' @examples
#' data(all_samples)
#' normalized_data <- normalizeMethylation_dep(df = all_samples, alpha = 0.001, normalize_position = TRUE, rescale = FALSE, plots = TRUE)
#' @importFrom ggplot2 ggplot geom_point geom_abline geom_smooth autoplot
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom dplyr select starts_with
#' @importFrom tidyr drop_na
#' @export
normalizeMethylation_dep <- function(df,
                                 alpha = 0.001,
                                 normalize_position = TRUE,
                                 rescale = FALSE,
                                 plots = TRUE) {
  temp_sample <- new("MicrobeMethyl",
                     df,
                     sample_name = "temp",
                     assembly = NA_character_,
                     sample_metadata = list(),
                     site_metadata = data.frame())
  normalized <- normalize_microbe_methyl_sample(temp_sample,
                                                alpha = alpha,
                                                normalize_position = normalize_position,
                                                rescale = rescale,
                                                plots = plots)
  position_col <- get_position_column(df)
  output_df <- dplyr::bind_cols(df[position_col],
                                as.data.frame(normalized$coverage),
                                as.data.frame(normalized$betas))
  colnames(output_df)[1] <- "Position"
  diagnostics <- normalized$diagnostics
  if (!is.null(diagnostics$plots)) {
    return(list(data = output_df,
                models = diagnostics$models,
                model_plots = diagnostics$plots$model_plots,
                cov_plots = diagnostics$plots$coverage_plots,
                r_squared_cov = diagnostics$r_squared_cov,
                r_squared_methyl = diagnostics$r_squared_methyl,
                methyl_plots = diagnostics$plots$methyl_plots,
                methyl_read_plots = diagnostics$plots$methyl_read_plots,
                RMSE = diagnostics$RMSE))
  }
  return(list(data = output_df,
              models = diagnostics$models,
              r_squared_cov = diagnostics$r_squared_cov,
              r_squared_methyl = diagnostics$r_squared_methyl,
              RMSE = diagnostics$RMSE))
}

#' @keywords internal
normalize_microbe_methyl_sample <- function(sample,
                                            alpha = 0.001,
                                            normalize_position = TRUE,
                                            rescale = FALSE,
                                            plots = TRUE) {
  sample_df <- as.data.frame(sample)
  position_col <- get_position_column(sample_df)
  coverage_cols <- get_coverage_columns(sample_df)
  methyl_cols <- get_methyl_columns(sample_df)
  if (length(position_col) == 0) {
    stop("No positional column found for normalization.")
  }
  if (length(coverage_cols) == 0) {
    stop("No coverage columns detected for normalization.")
  }
  if (length(methyl_cols) == 0) {
    stop("No methylation columns detected for normalization.")
  }
  complete_cols <- unique(c(position_col, coverage_cols, methyl_cols))
  complete_idx <- stats::complete.cases(sample_df[, complete_cols, drop = FALSE])
  complete_df <- sample_df[complete_idx, , drop = FALSE]
  coverage_df <- as.matrix(complete_df[, coverage_cols, drop = FALSE])
  if (normalize_position) {
    for (i in seq_len(ncol(coverage_df))) {
      coverage_df[, i] <- detrendCoverage(coverage_df[, i], complete_df[[position_col]])
    }
  }
  methyl_df <- as.matrix(complete_df[, methyl_cols, drop = FALSE])
  methyl_reads <- round(coverage_df * methyl_df * 100, digits = -2) / 100
  coverage_quant <- preprocessCore::normalize.quantiles(coverage_df)
  colnames(coverage_quant) <- coverage_cols
  colnames(methyl_reads) <- methyl_cols
  normalized_cov <- matrix(NA_real_, nrow = nrow(sample_df), ncol = ncol(coverage_quant))
  normalized_beta <- matrix(NA_real_, nrow = nrow(sample_df), ncol = ncol(methyl_df))
  colnames(normalized_cov) <- coverage_cols
  colnames(normalized_beta) <- methyl_cols
  transformabx <- function(a, b, x) {
    a + b * x
  }
  models <- vector("list", ncol(methyl_df))
  r_squared_cov <- vector("list", ncol(methyl_df))
  r_squared_methyl <- vector("list", ncol(methyl_df))
  rmse_list <- vector("list", ncol(methyl_df))
  if (plots) {
    model_plots <- vector("list", ncol(methyl_df))
    cov_plots <- vector("list", ncol(methyl_df))
    methyl_plots <- vector("list", ncol(methyl_df))
    methyl_read_plots <- vector("list", ncol(methyl_df))
  }
  coverage_original <- coverage_df
  for (i in seq_len(ncol(methyl_df))) {
    fit <- stats::lm(coverage_quant[, i] ~ coverage_original[, i])
    models[[i]] <- fit
    r_squared_cov[[i]] <- summary(fit)$r.squared
    a <- fit$coefficients[1]
    b <- fit$coefficients[2]
    normalized_reads <- transformabx(a, b, methyl_reads[, i])
    normalized_beta_complete <- normalized_reads / (coverage_quant[, i] + alpha)
    rmse_fit <- stats::lm(normalized_beta_complete ~ methyl_df[, i])
    rmse_pred <- seq(0, 1, length.out = length(methyl_df[, i]))
    rmse_list[[i]] <- sqrt(mean((rmse_fit$coefficients[2] * rmse_pred + rmse_fit$coefficients[1] - rmse_pred)^2))
    r_squared_methyl[[i]] <- summary(rmse_fit)$r.squared
    if (isTRUE(rescale)) {
      normalized_beta_complete <- (normalized_beta_complete - rmse_fit$coefficients[1]) / rmse_fit$coefficients[2]
      normalized_beta_complete[normalized_beta_complete < 0] <- 0
      normalized_beta_complete[normalized_beta_complete > 1] <- 1
    }
    normalized_beta[complete_idx, i] <- normalized_beta_complete
    if (plots) {
      model_plots[[i]] <- ggplot2::autoplot(fit)
      cov_plots[[i]] <- ggplot2::ggplot(data = data.frame(coverage = coverage_original[, i],
                                                          norm_coverage = coverage_quant[, i]),
                                        ggplot2::aes(x = coverage, y = norm_coverage)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 1, intercept = 0, color = "gray") +
        ggplot2::geom_smooth(method = "lm", color = "red")
      methyl_plots[[i]] <- ggplot2::ggplot(data = data.frame(beta = methyl_df[, i],
                                                             norm_beta = normalized_beta_complete),
                                           ggplot2::aes(x = beta, y = norm_beta)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 1, intercept = 0, color = "gray") +
        ggplot2::geom_smooth(method = "lm", color = "red")
      methyl_read_plots[[i]] <- ggplot2::ggplot(data = data.frame(methyl_reads = methyl_reads[, i],
                                                                  norm_methyl_reads = normalized_reads),
                                                ggplot2::aes(x = methyl_reads, y = norm_methyl_reads)) +
        ggplot2::geom_point() +
        ggplot2::geom_abline(slope = 1, intercept = 0, color = "gray") +
        ggplot2::geom_smooth(method = "lm", color = "red")
    }
  }
  normalized_cov[complete_idx, ] <- coverage_quant
  diagnostics <- list(models = models,
                      r_squared_cov = r_squared_cov,
                      r_squared_methyl = r_squared_methyl,
                      RMSE = rmse_list)
  if (plots) {
    diagnostics$plots <- list(model_plots = model_plots,
                              coverage_plots = cov_plots,
                              methyl_plots = methyl_plots,
                              methyl_read_plots = methyl_read_plots)
  }
  list(coverage = normalized_cov,
       betas = normalized_beta,
       diagnostics = diagnostics,
       position = sample_df[[position_col]])
}

get_position_column <- function(df) {
  position_candidates <- c("Position", "position", "start")
  position_candidates[position_candidates %in% colnames(df)][1]
}

get_coverage_columns <- function(df) {
  coverage_cols <- grep("^cov", colnames(df), value = TRUE)
  if (length(coverage_cols) == 0) {
    coverage_cols <- intersect(c("coverage", "cov"), colnames(df))
  }
  coverage_cols
}

get_methyl_columns <- function(df) {
  methyl_cols <- grep("^methyl", colnames(df), value = TRUE)
  if (length(methyl_cols) == 0) {
    methyl_cols <- intersect(c("percentMethylation", "beta"), colnames(df))
  }
  methyl_cols
}
