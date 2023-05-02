#' Build microbeMethyl object from bedMethyl file
#'
#' Reads a bedMethyl file (as defined by the ENCODE consortium) and creates
#' a microbeMethyl object containing methylation data and sample metadata for a
#' single sample.
#'
#' @param bedMethyl_file Path to the bedMethyl file (a tab-separated file).
#' @param sample_name A character string representing the name of the sample.
#' @param metadata A list containing metadata for the sample, such as
#'   experimental conditions and any other relevant information (optional).
#'
#' @return A microbeMethyl object.
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table fread
#' @importFrom S4Vectors Rle
buildMicrobeMethyl <- function(bedMethyl_file, sample_name, metadata = list()) {
  # Read the bedMethyl file using data.table's fread function
  bedMethyl_data <- data.table::fread(bedMethyl_file)

  # Check if the input file has the required columns
  required_columns <- c("chrom", "start", "end", "name", "score", "strand",
                        "thickStart", "thickEnd", "itemRgb", "coverage", "percentMethylation")
  if (!all(required_columns %in% colnames(bedMethyl_data))) {
    stop("The input file is missing one or more required columns.")
  }

  # Create a GRanges object from the chrom, start, and end columns
  granges <- GRanges(
    seqnames = Rle(bedMethyl_data$chrom),
    ranges = IRanges(start = bedMethyl_data$start, end = bedMethyl_data$end),
    strand = Rle(bedMethyl_data$strand)
  )

  # Extract coverage and percent methylation from the bedMethyl data
  num_reads <- bedMethyl_data$coverage
  beta_values <- bedMethyl_data$percentMethylation

  # Create a microbeMethyl object
  microbeMethyl_obj <- new("microbeMethyl",
                           sample_name = sample_name,
                           beta_values = beta_values,
                           num_reads = num_reads,
                           metadata = metadata,
                           GRanges = granges)

  return(microbeMethyl_obj)
}

#' Add a microbeMethyl object to a microbeMethylExperiment object
#'
#' Adds a microbeMethyl object to the list of samples in a
#' microbeMethylExperiment object.
#'
#' @param experiment A microbeMethylExperiment object.
#' @param sample A microbeMethyl object.
#'
#' @return The updated microbeMethylExperiment object.
#' @export
addSampleToExperiment <- function(experiment, sample) {
  if (!inherits(sample, "microbeMethyl")) {
    stop("The 'sample' argument must be a microbeMethyl object.")
  }

  if (!inherits(experiment, "microbeMethylExperiment")) {
    stop("The 'experiment' argument must be a microbeMethylExperiment object.")
  }

  experiment@samples[[length(experiment@samples) + 1]] <- sample
  return(experiment)
}

#' Combine multiple microbeMethyl objects into a microbeMethylExperiment object
#'
#' Takes a list of microbeMethyl objects and combines them into a single
#' microbeMethylExperiment object.
#'
#' @param microbeMethyl_list A list of microbeMethyl objects.
#' @param experiment_metadata A list containing metadata for the experiment (optional).
#'
#' @return A microbeMethylExperiment object containing all the input microbeMethyl objects.
#' @export
combineMicrobeMethyls <- function(microbeMethyl_list, experiment_metadata = list()) {
  # Check if all elements in the list are microbeMethyl objects
  if (!all(sapply(microbeMethyl_list, function(x) inherits(x, "microbeMethyl")))) {
    stop("All elements in 'microbeMethyl_list' must be microbeMethyl objects.")
  }

  # Create a microbeMethylExperiment object
  microbeMethyl_experiment <- new("microbeMethylExperiment",
                                  samples = microbeMethyl_list,
                                  experiment_metadata = experiment_metadata)

  return(microbeMethyl_experiment)
}
