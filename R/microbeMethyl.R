#' MicrobeMethyl class
#'
#' An S4 class to represent adenine methylation (6mA) data from a single
#' bacterial sample, particularly in E. coli. This class stores methylation
#' data (beta values and number of reads), sample name, and sample metadata.
#'
#' @slot sample_name A character string representing the name of the sample.
#' @slot metadata A list containing metadata for the sample, such as
#'   experimental conditions and any other relevant information.
#' @slot assembly A character string with the chromosome name.
#'
#' @name MicrobeMethyl-class
#' @aliases MicrobeMethyl
#' @docType class
#' @exportClass MicrobeMethyl
setClass("MicrobeMethyl",
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
#' \code{\linkS4class{MicrobeMethyl}} class.
#'
#' @param assembly A character string indicating the chromosome name.
#' @param start A numeric vector with start positions for each methylated base.
#' @param end A numeric vector with end positions for each methylated base.
#' @param coverage A numeric vector with coverage values for each position.
#' @param percentMethylation A double vector with percent methylation
#' (in decimals 0-1) for each position.
#' @param sample_name A character string indicating the sample name.
#'
#' @return A MicrobeMethyl object
#' @export
#'
#' @examples
#' microbe_methyl_obj <- MicrobeMethyl(assembly = "NC_000913.3",
#'                                     start = c(620, 621, 727, 728),
#'                                     end = c(620, 621, 727, 728),
#'                                     coverage = c(41, 89, 110, 76),
#'                                     percentMethylation = c(0.99, 0.95, 1, 0.84),
#'                                     sample_name = "WT1")
#'
MicrobeMethyl <- function(assembly, start, end,
                          coverage, percentMethylation, sample_name) {
  data <- data.frame(start = start, end = end, coverage = coverage,
                     percentMethylation = percentMethylation)
  new("MicrobeMethyl", data,
      sample_name = sample_name, assembly = assembly,
      sample_metadata = list(), site_metadata = data.frame())
}

#' Build MicrobeMethyl object from bedMethyl file
#'
#' Reads a bedMethyl file (as defined by the ENCODE consortium) and creates
#' a MicrobeMethyl object containing methylation data and sample metadata for a
#' single sample.
#'
#' @param bedMethyl_file Path to the bedMethyl file (a tab-separated file).
#' @param sample_name A character string representing the name of the sample.
#' @param metadata A list containing metadata for the sample, such as
#'   experimental conditions and any other relevant information (optional).
#'
#' @return A MicrobeMethyl object.
#' @export
#' @importFrom data.table fread
buildMicrobeMethyl <- function(bedMethyl_file, sample_name, metadata) {
  # Read the bedMethyl file using data.table's fread function
  bedMethyl_data <- data.table::fread(bedMethyl_file)

  # Check if the input file has the required columns
  required_columns <- c("chrom", "start", "end", "name", "score", "strand",
                        "thickStart", "thickEnd", "itemRgb", "coverage", "percentMethylation")
  if (!all(required_columns %in% colnames(bedMethyl_data))) {
    stop("The input file is missing one or more required columns.")
  }

  # Create a MicrobeMethyl object
  MicrobeMethyl_obj <- MicrobeMethyl(assembly = unique(bedMethyl_data$chrom),
                                     start = bedMethyl_data$start,
                                     end = bedMethyl_data$end,
                                     coverage = bedMethyl_data$coverage,
                                     percentMethylation = bedMethyl_data$percentMethylation,
                                     sample_name = sample_name)

  MicrobeMethyl_obj@sample_metadata <- metadata

  return(MicrobeMethyl_obj)
}

#' Get beta values from MicrobeMethyl or MicrobeMethylExperiment objects
#'
#' This generic function retrieves the beta values from MicrobeMethyl or
#' MicrobeMethylExperiment objects.
#'
#' @param object An object of class MicrobeMethyl or MicrobeMethylExperiment
#' @return A vector of beta values for a MicrobeMethyl object, or a list of
#'   beta value vectors for a MicrobeMethylExperiment object
#' @export
setGeneric("getBetaValues",
           function(object) {
             standardGeneric("getBetaValues")
           })

#' @title getBetaValues for MicrobeMethyl objects
#' @description Retrieve the beta values from a MicrobeMethyl object.
#' @param object A MicrobeMethyl object
#' @return A numeric vector of beta values for the given MicrobeMethyl object
#' @exportMethod getBetaValues
#' @aliases getBetaValues,MicrobeMethyl-method
setMethod("getBetaValues",
          signature = "MicrobeMethyl",
          definition = function(object) {
            return(object@beta_values)
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
setClass("MicrobeMethylExperiment",
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
#' microbe_methyl_obj1 <- MicrobeMethyl(assembly = "NC_000913.3",
#'                                      start = c(620, 621, 727, 728),
#'                                      end = c(620, 621, 727, 728),
#'                                      coverage = c(41, 89, 110, 76),
#'                                      percentMethylation = c(0.99, 0.95, 1, 0.84),
#'                                      sample_name = "WT1")
#' microbe_methyl_obj2 <- MicrobeMethyl(assembly = "NC_000913.3",
#'                                      start = c(630, 631, 737, 738),
#'                                      end = c(630, 631, 737, 738),
#'                                      coverage = c(45, 90, 100, 80),
#'                                      percentMethylation = c(0.98, 0.92, 1, 0.82),
#'                                      sample_name = "WT2")
#' experiment_metadata <- list(date = "2023-05-02", condition = "wild_type")
#' experiment_obj <- MicrobeMethylExperiment(microbe_methyl_obj1, microbe_methyl_obj2,
#'                                           experiment_metadata = experiment_metadata)
#'
MicrobeMethylExperiment <- function(..., experiment_metadata) {
  microbe_methyl_objs <- list(...)

  # Check that all input objects are MicrobeMethyl objects
  for (obj in microbe_methyl_objs) {
    if (!is(obj, "MicrobeMethyl")) {
      stop("All input objects must be of class 'MicrobeMethyl'")
    }
  }

  # Create a named list of MicrobeMethyl objects using the sample_name as names
  samples <- setNames(microbe_methyl_objs, sapply(microbe_methyl_objs, function(x) x@sample_name))

  # Create the MicrobeMethylExperiment object
  new("MicrobeMethylExperiment",
      samples = samples,
      experiment_metadata = experiment_metadata)
}

#' Generic function for sampleNames
#'
#' @param object An object for which to extract sample names.
#'
#' @return A character vector of sample names.
#' @export
setGeneric("sampleNames",
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
#' # Using the example MicrobeMethylExperiment object created earlier
#' sample_names <- sampleNames(experiment_obj)
#' print(sample_names)
#'
setMethod("sampleNames", "MicrobeMethylExperiment", function(object) {
  return(names(object@samples))
})




#' @title getBetaValues for MicrobeMethylExperiment objects
#' @description Retrieve the beta values for all samples in a
#'   MicrobeMethylExperiment object.
#' @param object A MicrobeMethylExperiment object
#' @return A list of numeric vectors, where each vector contains the beta values
#'   for a sample in the MicrobeMethylExperiment object
#' @exportMethod getBetaValues
#' @aliases getBetaValues,MicrobeMethylExperiment-method
setMethod("getBetaValues",
          signature = "MicrobeMethylExperiment",
          definition = function(object) {
            beta_values_list <- lapply(object@samples, function(sample) {
              sample@beta_values
            })
            return(beta_values_list)
          })
