readBedMethyl <- function(input, sample_name, sample_metadata = list(), assembly_length = NULL,
                          genome_name = NULL, is_circular = FALSE, zero_based = FALSE) {
  # Read the data
  if (is.character(input)) {
    # Assume input is a file path
    bedMethyl_data <- data.table::fread(input)
  } else if (is.data.frame(input)) {
    bedMethyl_data <- input
  } else {
    stop("'input' must be a file path or a data frame")
  }

  # Check required columns
  required_columns <- c("chrom", "start", "end", "coverage", "percentMethylation", "strand")
  if (!all(required_columns %in% colnames(bedMethyl_data))) {
    missing_cols <- required_columns[!required_columns %in% colnames(bedMethyl_data)]
    stop("The input data is missing the following required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Adjust coordinates if zero-based
  if (zero_based) {
    # Convert from zero-based to one-based coordinates
    start_positions <- bedMethyl_data$start + 1
    end_positions <- bedMethyl_data$end + 1
  } else {
    start_positions <- bedMethyl_data$start
    end_positions <- bedMethyl_data$end
  }
  return(bedMethyl_data)
}
