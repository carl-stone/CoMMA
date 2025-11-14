library(testthat)
library(CoMMA)

build_test_sample <- function(sample_name, coverage_shift = 0) {
  df <- data.frame(
    Position = 1:4,
    cov_rep1 = c(10, 20, 30, 40) + coverage_shift,
    cov_rep2 = c(5, 15, 25, 35) + coverage_shift,
    methyl_rep1 = c(0.1, 0.2, 0.3, 0.4),
    methyl_rep2 = c(0.2, 0.3, 0.4, 0.5)
  )
  new("MicrobeMethyl",
      df,
      sample_name = sample_name,
      assembly = "asm",
      sample_metadata = list(),
      site_metadata = data.frame())
}

compute_expected_betas <- function(df, alpha = 0.001) {
  coverage_cols <- grep("^cov", colnames(df), value = TRUE)
  methyl_cols <- grep("^methyl", colnames(df), value = TRUE)
  coverage_mat <- as.matrix(df[, coverage_cols])
  methyl_mat <- as.matrix(df[, methyl_cols])
  coverage_quant <- preprocessCore::normalize.quantiles(coverage_mat)
  methyl_reads <- round(coverage_mat * methyl_mat * 100, digits = -2) / 100
  out <- matrix(NA_real_, nrow = nrow(methyl_mat), ncol = ncol(methyl_mat))
  for (i in seq_len(ncol(methyl_mat))) {
    fit <- stats::lm(coverage_quant[, i] ~ coverage_mat[, i])
    normalized_reads <- fit$coefficients[1] + fit$coefficients[2] * methyl_reads[, i]
    normalized_beta <- normalized_reads / (coverage_quant[, i] + alpha)
    rmse_fit <- stats::lm(normalized_beta ~ methyl_mat[, i])
    normalized_beta <- (normalized_beta - rmse_fit$coefficients[1]) / rmse_fit$coefficients[2]
    normalized_beta[normalized_beta < 0] <- 0
    normalized_beta[normalized_beta > 1] <- 1
    out[, i] <- normalized_beta
  }
  colnames(out) <- methyl_cols
  out
}


test_that("normalizeMethylation normalizes coverage and betas", {
  sample_one <- build_test_sample("sample_one")
  sample_two <- build_test_sample("sample_two", coverage_shift = 5)
  experiment <- MicrobeMethylExperiment(sample_one, sample_two, experiment_metadata = list())

  normalized <- normalizeMethylation(experiment, rescale = TRUE, plots = TRUE)
  normalized_sample <- normalized@samples[[1]]
  normalized_df <- as.data.frame(normalized_sample)
  coverage_cols <- grep("^cov", colnames(normalized_df), value = TRUE)
  methyl_cols <- grep("^methyl", colnames(normalized_df), value = TRUE)

  original_cov <- as.matrix(as.data.frame(sample_one)[, coverage_cols])
  expected_cov <- preprocessCore::normalize.quantiles(original_cov)
  expect_equal(as.matrix(normalized_df[, coverage_cols]), expected_cov)

  expected_betas <- compute_expected_betas(as.data.frame(sample_one))
  expect_equal(as.matrix(normalized_df[, methyl_cols]), expected_betas)

  normalization_meta <- normalized_sample@sample_metadata$normalization
  expect_type(normalization_meta, "list")
  expect_true(all(c("models", "r_squared_cov", "r_squared_methyl", "RMSE") %in% names(normalization_meta)))
  expect_true("plots" %in% names(normalization_meta))
  expect_length(normalization_meta$plots$model_plots, length(methyl_cols))
})
