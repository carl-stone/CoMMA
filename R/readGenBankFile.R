#' Read a GenBank File and Extract Genomic Information
#'
#' This function reads a GenBank file (assuming bacterial genome) and extracts genomic information,
#' including assembly lengths, features (genes, CDS, tRNA, rRNA, etc.), and sequence information.
#' It returns a list containing assembly information and features as \code{GRanges} objects.
#'
#' @param genbank_file A character string specifying the path to the GenBank file.
#' @param include_sequence Logical indicating whether to include the sequence data (default is \code{FALSE}).
#'
#' @return A list containing:
#' \item{assembly_info}{A \code{Seqinfo} object with assembly names, lengths, and circularity information.}
#' \item{features}{A \code{GRanges} object with genomic features (genes, CDS, tRNA, rRNA, etc.).}
#' \item{sequence}{An optional \code{DNAStringSet} object with the genome sequence (if \code{include_sequence = TRUE}).}
#' @export
#' @import GenomicRanges
#' @import Biostrings
#'
#' @examples
#' # Read a GenBank file
#' genbank_data <- readGenBankFile("path/to/genome.gbff")
#'
#' # Access assembly information
#' assembly_info <- genbank_data$assembly_info
#'
#' # Access genomic features
#' features <- genbank_data$features
#'
#' # If sequence data was included
#' if (!is.null(genbank_data$sequence)) {
#'   genome_sequence <- genbank_data$sequence
#' }
readGenBankFile <- function(genbank_file, include_sequence = FALSE) {
  if (!file.exists(genbank_file)) {
    stop("The GenBank file does not exist: ", genbank_file)
  }

  # Read the GenBank file into a character vector
  gb_lines <- readLines(genbank_file)

  # Initialize variables
  locus_line <- NULL
  origin_index <- NULL
  features_start <- NULL
  features_end <- NULL
  sequence <- NULL

  # Parse the LOCUS line
  for (i in seq_along(gb_lines)) {
    line <- gb_lines[i]
    if (startsWith(line, "LOCUS")) {
      locus_line <- line
    } else if (startsWith(line, "FEATURES")) {
      features_start <- i
    } else if (startsWith(line, "ORIGIN")) {
      origin_index <- i
      features_end <- i - 1
      break
    }
  }

  if (is.null(locus_line)) {
    stop("The GenBank file does not contain a LOCUS line.")
  }

  # Extract assembly information from the LOCUS line
  locus_parts <- unlist(strsplit(locus_line, "\\s+"))
  locus_parts <- locus_parts[locus_parts != ""]

  accession <- locus_parts[2]
  sequence_length <- as.numeric(locus_parts[3])
  topology <- locus_parts[5]
  is_circular <- tolower(topology) == "circular"
  genome_name <- NULL  # We'll extract it from the ORGANISM line

  # Extract genome name from the ORGANISM line
  organism_line_index <- grep("^  ORGANISM", gb_lines)
  if (length(organism_line_index) > 0) {
    organism_line <- gb_lines[organism_line_index[1]]
    organism_name <- trimws(sub("^  ORGANISM\\s+", "", organism_line))
    genome_name <- organism_name
  }

  # Create Seqinfo object
  assembly_info <- Seqinfo(
    seqnames = accession,
    seqlengths = sequence_length,
    isCircular = is_circular,
    genome = genome_name
  )

  # Extract features
  if (is.null(features_start) || is.null(features_end)) {
    stop("The GenBank file does not contain FEATURES section.")
  }

  features_lines <- gb_lines[features_start:features_end]
  features_gr <- parseGenBankFeatures(features_lines, accession)

  # Extract sequence if requested
  if (include_sequence) {
    if (is.null(origin_index)) {
      stop("The GenBank file does not contain ORIGIN section for sequence data.")
    }
    sequence_lines <- gb_lines[(origin_index + 1):length(gb_lines)]
    sequence <- parseGenBankSequence(sequence_lines)
    names(sequence) <- accession
    sequence <- DNAStringSet(sequence)
  } else {
    sequence <- NULL
  }

  # Return the extracted information as a list
  return(list(
    assembly_info = assembly_info,
    features = features_gr,
    sequence = sequence
  ))
}
