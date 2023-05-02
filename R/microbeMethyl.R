#' microbeMethyl class
#'
#' An S4 class to represent adenine methylation (6mA) data from a single
#' bacterial sample, particularly in E. coli. This class stores methylation
#' data (beta values and number of reads), sample name, and sample metadata.
#'
#' @slot sample_name A character string representing the name of the sample.
#' @slot beta_values A numeric vector of beta values representing the
#'   methylation levels for each genomic position.
#' @slot num_reads An integer vector of read counts for each genomic
#'   position, representing the number of reads covering each position.
#' @slot metadata A list containing metadata for the sample, such as
#'   experimental conditions and any other relevant information.
#'
#' @name microbeMethyl-class
#' @aliases microbeMethyl
#' @docType class
#' @importFrom GenomicRanges GRanges
#' @exportClass microbeMethyl
setClass("microbeMethyl",
         slots = list(
           sample_name = "character",
           beta_values = "numeric",
           num_reads = "integer",
           metadata = "list"
         ),
         contains = "GRanges"
)

#' Get beta values from microbeMethyl or microbeMethylExperiment objects
#'
#' This generic function retrieves the beta values from microbeMethyl or
#' microbeMethylExperiment objects.
#'
#' @param object An object of class microbeMethyl or microbeMethylExperiment
#' @return A vector of beta values for a microbeMethyl object, or a list of
#'   beta value vectors for a microbeMethylExperiment object
#' @export
setGeneric("getBetaValues",
           function(object) {
             standardGeneric("getBetaValues")
           })

#' @title getBetaValues for microbeMethyl objects
#' @description Retrieve the beta values from a microbeMethyl object.
#' @param object A microbeMethyl object
#' @return A numeric vector of beta values for the given microbeMethyl object
#' @exportMethod getBetaValues
#' @aliases getBetaValues,microbeMethyl-method
setMethod("getBetaValues",
          signature = "microbeMethyl",
          definition = function(object) {
            return(object@beta_values)
          })

#' microbeMethylExperiment class
#'
#' An S4 class to represent an experiment containing multiple microbeMethyl
#' objects. This class stores a list of microbeMethyl objects and experiment
#' metadata.
#'
#' @slot samples A list of microbeMethyl objects, each representing a single
#'   sample from the same experiment.
#' @slot experiment_metadata A list containing metadata for the experiment,
#'   such as experimental conditions, date, and any other relevant information.
#'
#' @name microbeMethylExperiment-class
#' @aliases microbeMethylExperiment
#' @docType class
#' @exportClass microbeMethylExperiment
setClass("microbeMethylExperiment",
         slots = list(
           samples = "list",
           experiment_metadata = "list"
         )
)


#' @title getBetaValues for microbeMethylExperiment objects
#' @description Retrieve the beta values for all samples in a
#'   microbeMethylExperiment object.
#' @param object A microbeMethylExperiment object
#' @return A list of numeric vectors, where each vector contains the beta values
#'   for a sample in the microbeMethylExperiment object
#' @exportMethod getBetaValues
#' @aliases getBetaValues,microbeMethylExperiment-method
setMethod("getBetaValues",
          signature = "microbeMethylExperiment",
          definition = function(object) {
            beta_values_list <- lapply(object@samples, function(sample) {
              sample@beta_values
            })
            return(beta_values_list)
          })
