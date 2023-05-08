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
              df <- data.frame(Position = start(sample),
                               coverage = sample@num_reads,
                               beta = sample@beta_values)

              # Rest of the normalizeMethylation code as a function of df
              # ...
              # At the end of the function, return the updated MicrobeMethyl object
              return(updated_sample)
            })

            # Update the MicrobeMethylExperiment object with the normalized samples
            object@samples <- normalized_samples

            return(object)
          }
)
