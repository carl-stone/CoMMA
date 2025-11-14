#' @rdname normalizeMethylation
setGeneric("normalizeMethylation",
           function(object, alpha = 0.001, normalize_position = TRUE, rescale = FALSE, plots = TRUE) {
             standardGeneric("normalizeMethylation")
           }
)

#' @rdname normalizeMethylation
#' @aliases normalizeMethylation,MicrobeMethylExperiment-method
#' @importFrom dplyr bind_cols
#' @importFrom tidyr drop_na
#' @exportMethod normalizeMethylation
setMethod("normalizeMethylation", "MicrobeMethylExperiment",
          function(object, alpha = 0.001, normalize_position = TRUE, rescale = FALSE, plots = TRUE) {
            samples <- object@samples
            normalized_samples <- lapply(samples, function(sample) {
              normalized <- normalize_microbe_methyl_sample(sample,
                                                            alpha = alpha,
                                                            normalize_position = normalize_position,
                                                            rescale = rescale,
                                                            plots = plots)
              sample_df <- as.data.frame(sample)
              coverage_cols <- get_coverage_columns(sample_df)
              methyl_cols <- get_methyl_columns(sample_df)
              if (length(coverage_cols) > 0) {
                sample_df[, coverage_cols] <- normalized$coverage
              }
              if (length(methyl_cols) > 0) {
                sample_df[, methyl_cols] <- normalized$betas
              }
              sample@.Data <- sample_df
              sample@sample_metadata$normalization <- normalized$diagnostics
              sample
            })

            object@samples <- normalized_samples
            object
          }
)
