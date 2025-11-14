#' Build MicrobeMethyl object from a bedMethyl file
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
buildMicrobeMethyl <- function(bedMethyl_file, sample_name, metadata = list()) {
  # Read the bedMethyl file using data.table's fread function
  bedMethyl_data <- data.table::fread(bedMethyl_file)

  # Check if the input file has the required columns
  required_columns <- c(
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "coverage", "percentMethylation"
  )
  if (!all(required_columns %in% colnames(bedMethyl_data))) {
    stop("The input file is missing one or more required columns.")
  }

  bed_df <- as.data.frame(bedMethyl_data)

  site_metadata_cols <- setdiff(
    names(bed_df),
    c("start", "end", "coverage", "percentMethylation")
  )

  site_metadata <- bed_df[, site_metadata_cols, drop = FALSE]

  MicrobeMethyl_obj <- MicrobeMethyl(
    assembly = unique(bed_df$chrom),
    start = bed_df$start,
    end = bed_df$end,
    coverage = bed_df$coverage,
    percentMethylation = bed_df$percentMethylation,
    sample_name = sample_name,
    sample_metadata = metadata,
    site_metadata = site_metadata
  )

  MicrobeMethyl_obj
}
