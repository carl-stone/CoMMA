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

#' Parse GenBank Features into GRanges
#'
#' This function parses the FEATURES section of a GenBank file and converts it into a \code{GRanges} object.
#'
#' @param features_lines A character vector containing the lines of the FEATURES section.
#' @param seqname The sequence name (e.g., accession number).
#'
#' @return A \code{GRanges} object with genomic features.
#' @export
#' @import GenomicRanges
#' @import IRanges
#' @importFrom dplyr bind_rows
parseGenBankFeatures <- function(features_lines, seqname) {
  # Initialize variables
  features_list <- list()
  current_feature <- NULL

  for (line in features_lines[-1]) {  # Skip the first line (FEATURES header)
    if (grepl("^\\s{5}\\S", line)) {
      # This is a new feature
      if (!is.null(current_feature)) {
        # Save the previous feature
        features_list <- c(features_list, list(current_feature))
      }
      # Start a new feature
      feature_key <- trimws(substr(line, 6, 20))
      location_string <- trimws(substr(line, 21, nchar(line)))
      current_feature <- list(
        key = feature_key,
        location = location_string,
        qualifiers = list()
      )
    } else if (grepl("^\\s{21}/", line)) {
      # This is a qualifier
      qualifier_line <- trimws(substr(line, 22, nchar(line)))
      qualifier_parts <- strsplit(qualifier_line, "=")[[1]]
      qualifier_name <- gsub("^/", "", qualifier_parts[1])
      qualifier_value <- ifelse(length(qualifier_parts) > 1, qualifier_parts[2], NA)
      # Remove surrounding quotes if present
      qualifier_value <- gsub('^"|"$', '', qualifier_value)
      # Append to qualifiers list
      current_feature$qualifiers[[qualifier_name]] <- qualifier_value
    } else if (grepl("^\\s{21}\\S", line)) {
      # Continuation of location or qualifier
      continuation_line <- trimws(substr(line, 22, nchar(line)))
      if (grepl("^/", continuation_line)) {
        # Continuation of a qualifier
        qualifier_line <- continuation_line
        qualifier_parts <- strsplit(qualifier_line, "=")[[1]]
        qualifier_name <- gsub("^/", "", qualifier_parts[1])
        qualifier_value <- ifelse(length(qualifier_parts) > 1, qualifier_parts[2], NA)
        # Remove surrounding quotes if present
        qualifier_value <- gsub('^"|"$', '', qualifier_value)
        # Append to qualifiers list
        current_feature$qualifiers[[qualifier_name]] <- qualifier_value
      } else {
        # Continuation of location
        current_feature$location <- paste0(current_feature$location, continuation_line)
      }
    }
  }
  # Add the last feature
  if (!is.null(current_feature)) {
    features_list <- c(features_list, list(current_feature))
  }

  # Now parse the features into GRanges
  gr_list <- lapply(features_list, function(feature) {
    # Parse the location string into ranges
    ranges_info <- parseLocationString(feature$location)
    if (is.null(ranges_info)) return(NULL)
    gr <- GRanges(
      seqnames = seqname,
      ranges = ranges_info$ranges,
      strand = ranges_info$strand,
      feature_type = feature$key,
      qualifiers = list(feature$qualifiers)
    )
    return(gr)
  })

  # Remove NULL entries
  gr_list <- gr_list[!sapply(gr_list, is.null)]

  # Combine all GRanges
  gr_all <- do.call(c, gr_list)

  # Flatten the qualifiers into metadata columns
  qualifiers_list <- lapply(gr_all$qualifiers, function(q) {
    as.data.frame(q, stringsAsFactors = FALSE)
  })
  qualifiers_df <- bind_rows(qualifiers_list)
  mcols(gr_all) <- cbind(mcols(gr_all), qualifiers_df)
  gr_all$qualifiers <- NULL  # Remove the qualifiers column

  return(gr_all)
}

#' Parse Location String from GenBank Feature
#'
#' This function parses a location string from a GenBank feature and returns ranges and strand information.
#'
#' @param location_string A character string representing the location.
#'
#' @return A list containing:
#' \item{ranges}{An \code{IRanges} object with the feature ranges.}
#' \item{strand}{A character representing the strand ("+", "-", or "*").}
parseLocationString <- function(location_string) {
  # Remove join and complement operators
  strand <- "+"
  loc_str <- location_string
  if (startsWith(loc_str, "complement(") && endsWith(loc_str, ")")) {
    strand <- "-"
    loc_str <- substr(loc_str, 12, nchar(loc_str) - 1)
  }
  if (startsWith(loc_str, "join(") && endsWith(loc_str, ")")) {
    loc_str <- substr(loc_str, 6, nchar(loc_str) - 1)
  }

  # Split by comma for multiple regions
  loc_parts <- strsplit(loc_str, ",")[[1]]

  # Parse each part into ranges
  starts <- integer()
  ends <- integer()

  for (part in loc_parts) {
    # Remove any < or > symbols indicating partial positions
    part <- gsub("[<>]", "", part)
    if (grepl("\\d+\\.\\.\\d+", part)) {
      # Range
      pos <- as.numeric(strsplit(part, "\\.\\.\\.?")[[1]])
      starts <- c(starts, pos[1])
      ends <- c(ends, pos[2])
    } else if (grepl("^\\d+$", part)) {
      # Single position
      pos <- as.numeric(part)
      starts <- c(starts, pos)
      ends <- c(ends, pos)
    } else {
      # Complex or unsupported location
      return(NULL)
    }
  }

  ranges <- IRanges(start = starts, end = ends)
  return(list(ranges = ranges, strand = strand))
}

#' Parse GenBank Sequence
#'
#' This function parses the sequence from the ORIGIN section of a GenBank file.
#'
#' @param sequence_lines A character vector containing the lines of the ORIGIN section.
#'
#' @return A character string representing the sequence.
parseGenBankSequence <- function(sequence_lines) {
  # Remove line numbers and spaces
  seq_chars <- gsub("[0-9\\s]", "", sequence_lines)
  sequence <- paste0(seq_chars, collapse = "")
  return(sequence)
}
