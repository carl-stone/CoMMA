#' MicrobeMethyl class
#'
#' An S4 class to represent adenine methylation (6mA) data from a single
#' bacterial sample, particularly in *E. coli*. This class stores methylation
#' data (beta values and number of reads), sample name, and sample metadata.
#'
#' @slot sample_name A character string representing the name of the sample.
#' @slot assembly A character vector with the chromosome name(s) contained in
#'   the object.
#' @slot sample_metadata A list containing metadata for the sample, such as
#'   experimental conditions and any other relevant information.
#' @slot site_metadata A data.frame containing additional site-level annotation
#'   that accompanies the methylation measurements.
#'
#' @name MicrobeMethyl-class
#' @aliases MicrobeMethyl
#' @docType class
#' @exportClass MicrobeMethyl
setClass(
  "MicrobeMethyl",
  contains = "data.frame",
  slots = c(
    sample_name = "character",
    assembly = "character",
    sample_metadata = "list",
    site_metadata = "data.frame"
  ),
  prototype = prototype(
    sample_name = NA_character_,
    assembly = NA_character_,
    sample_metadata = list(),
    site_metadata = data.frame()
  )
)

#' MicrobeMethyl constructor
#'
#' @description This is a constructor function for the
#'   \code{\linkS4class{MicrobeMethyl}} class.
#'
#' @param assembly A character vector indicating the chromosome name(s).
#' @param start A numeric vector with start positions for each methylated base.
#' @param end A numeric vector with end positions for each methylated base.
#' @param coverage A numeric vector with coverage values for each position.
#' @param percentMethylation A double vector with percent methylation
#'   (in decimals 0-1) for each position.
#' @param sample_name A character string indicating the sample name.
#' @param sample_metadata A list containing metadata for the sample.
#' @param site_metadata A data.frame containing site-level metadata.
#'
#' @return A MicrobeMethyl object
#' @export
#'
#' @examples
#' microbe_methyl_obj <- MicrobeMethyl(
#'   assembly = "NC_000913.3",
#'   start = c(620, 621, 727, 728),
#'   end = c(620, 621, 727, 728),
#'   coverage = c(41, 89, 110, 76),
#'   percentMethylation = c(0.99, 0.95, 1, 0.84),
#'   sample_name = "WT1",
#'   sample_metadata = list(condition = "WT")
#' )
MicrobeMethyl <- function(assembly, start, end,
                          coverage, percentMethylation, sample_name,
                          sample_metadata = list(),
                          site_metadata = data.frame()) {
  data <- data.frame(
    start = start,
    end = end,
    coverage = coverage,
    percentMethylation = percentMethylation
  )
  new(
    "MicrobeMethyl",
    data,
    sample_name = sample_name,
    assembly = assembly,
    sample_metadata = sample_metadata,
    site_metadata = site_metadata
  )
}

# Validator function for MicrobeMethyl object
# 1) All columns in .Data must be the same length.
#
setValidity("MicrobeMethyl", function(object) {
  if (length(unique(lapply(object@.Data, length))) != 1) {
    "All columns in methylation data must be the same length"
  } else if (!all(as.logical(lapply(object@.Data, is.numeric)))) {
    "All columns in methylation data must be numeric"
  } else {
    TRUE
  }
})

#' MicrobeMethylExperiment class
#'
#' An S4 class to represent an experiment containing multiple MicrobeMethyl
#' objects. This class stores a list of MicrobeMethyl objects and experiment
#' metadata.
#'
#' @slot samples A list of MicrobeMethyl objects, each representing a single
#'   sample from the same experiment.
#' @slot experiment_metadata A list containing metadata for the experiment,
#'   such as experimental conditions, date, and any other relevant information.
#'
#' @name MicrobeMethylExperiment-class
#' @aliases MicrobeMethylExperiment
#' @docType class
#' @exportClass MicrobeMethylExperiment
setClass(
  "MicrobeMethylExperiment",
  slots = c(
    samples = "list",
    experiment_metadata = "list"
  )
)

#' MicrobeMethylExperiment constructor
#'
#' @description This is a constructor function for the
#'   \code{\linkS4class{MicrobeMethylExperiment}} class.
#'
#' @param ... One or more MicrobeMethyl objects.
#' @param experiment_metadata A list containing metadata for the experiment,
#'   such as experimental conditions, date, and any other relevant information.
#'
#' @return A MicrobeMethylExperiment object
#' @export
#'
#' @examples
#' microbe_methyl_obj1 <- MicrobeMethyl(
#'   assembly = "NC_000913.3",
#'   start = c(620, 621, 727, 728),
#'   end = c(620, 621, 727, 728),
#'   coverage = c(41, 89, 110, 76),
#'   percentMethylation = c(0.99, 0.95, 1, 0.84),
#'   sample_name = "WT1"
#' )
#' microbe_methyl_obj2 <- MicrobeMethyl(
#'   assembly = "NC_000913.3",
#'   start = c(630, 631, 737, 738),
#'   end = c(630, 631, 737, 738),
#'   coverage = c(45, 90, 100, 80),
#'   percentMethylation = c(0.98, 0.92, 1, 0.82),
#'   sample_name = "WT2"
#' )
#' experiment_metadata <- list(date = "2023-05-02", condition = "wild_type")
#' experiment_obj <- MicrobeMethylExperiment(
#'   microbe_methyl_obj1, microbe_methyl_obj2,
#'   experiment_metadata = experiment_metadata
#' )
MicrobeMethylExperiment <- function(..., experiment_metadata) {
  microbe_methyl_objs <- list(...)

  # Check that all input objects are MicrobeMethyl objects
  for (obj in microbe_methyl_objs) {
    if (!is(obj, "MicrobeMethyl")) {
      stop("All input objects must be of class 'MicrobeMethyl'")
    }
  }

  # Create a named list of MicrobeMethyl objects using the sample_name as names
  samples <- setNames(
    microbe_methyl_objs,
    vapply(microbe_methyl_objs, slot, character(1), name = "sample_name")
  )

  # Create the MicrobeMethylExperiment object
  new(
    "MicrobeMethylExperiment",
    samples = samples,
    experiment_metadata = experiment_metadata
  )
}

#' Generic function for sampleNames
#'
#' @param object An object for which to extract sample names.
#'
#' @return A character vector of sample names.
#' @export
setGeneric(
  "sampleNames",
  function(object) {
    standardGeneric("sampleNames")
  }
)

#' Extract sample names from a MicrobeMethylExperiment object
#'
#' @description This method extracts the sample names from a
#'   \code{\linkS4class{MicrobeMethylExperiment}} object.
#'
#' @param object A MicrobeMethylExperiment object.
#'
#' @return A character vector of sample names.
#' @export
#'
#' @examples
#' sample_names <- sampleNames(experiment_obj)
#' print(sample_names)
setMethod("sampleNames", "MicrobeMethylExperiment", function(object) {
  names(object@samples)
})

#' Get percent methylation values from MicrobeMethyl or MicrobeMethylExperiment
#'
#' This generic function retrieves the percent methylation values (beta values)
#' from MicrobeMethyl or MicrobeMethylExperiment objects.
#'
#' @param object An object of class MicrobeMethyl or MicrobeMethylExperiment.
#'
#' @return A numeric vector of beta values for a MicrobeMethyl object, or a
#'   list of beta value vectors for a MicrobeMethylExperiment object.
#' @export
setGeneric("getBetaValues", function(object) standardGeneric("getBetaValues"))

#' @title getBetaValues for MicrobeMethyl objects
#' @description Retrieve the percent methylation values from a MicrobeMethyl
#'   object.
#' @param object A MicrobeMethyl object.
#' @return A numeric vector of percent methylation values for the given
#'   MicrobeMethyl object.
#' @exportMethod getBetaValues
#' @aliases getBetaValues,MicrobeMethyl-method
setMethod(
  "getBetaValues",
  signature = "MicrobeMethyl",
  definition = function(object) {
    object$percentMethylation
  }
)

#' @title getBetaValues for MicrobeMethylExperiment objects
#' @description Retrieve the percent methylation values for all samples in a
#'   MicrobeMethylExperiment object.
#' @param object A MicrobeMethylExperiment object.
#' @return A list of numeric vectors, where each vector contains the percent
#'   methylation values for a sample in the MicrobeMethylExperiment object.
#' @exportMethod getBetaValues
#' @aliases getBetaValues,MicrobeMethylExperiment-method
setMethod(
  "getBetaValues",
  signature = "MicrobeMethylExperiment",
  definition = function(object) {
    lapply(object@samples, function(sample) {
      sample$percentMethylation
    })
  }
)
