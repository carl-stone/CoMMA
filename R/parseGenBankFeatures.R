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
parseGenBankFeatures <- function(features_lines, seqname) {
  # Initialize variables
  features_list <- list()
  current_feature <- NULL
  current_qualifiers <- list()

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
      qualifiers = feature$qualifiers
    )
    return(gr)
  })

  # Combine all GRanges
  gr_all <- do.call(c, gr_list)

  # Flatten the qualifiers into metadata columns
  mcols(gr_all) <- cbind(mcols(gr_all), do.call(rbind, lapply(gr_all$qualifiers, as.data.frame, stringsAsFactors = FALSE)))
  gr_all$qualifiers <- NULL  # Remove the qualifiers column

  return(gr_all)
}
