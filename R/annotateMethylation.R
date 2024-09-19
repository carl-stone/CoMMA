#' Annotate Methylation Data with Features from a GenBank File
#'
#' This function reads a GenBank file, extracts genomic features, and annotates
#' a methylation \code{GRanges} object with overlapping genomic features and assembly information.
#'
#' @param methylation_gr A \code{GRanges} object containing methylation data.
#' @param genbank_file A character string specifying the path to the GenBank file.
#' @param feature_types A character vector specifying the types of features to include
#'   (e.g., \code{c("genes", "cds", "ncRNA")}).
#'
#' @return A \code{GRanges} object with methylation data annotated with genomic features
#'   and updated assembly information.
#' @export
#' @importFrom GenomicFeatures makeTxDbFromGenBank genes transcripts exons cds
#' @importFrom GenomicRanges seqlevels seqinfo findOverlaps
#' @importFrom GenomeInfoDb seqlevels<-
#' @importFrom S4Vectors mcols mcols<-
#' @importFrom rtracklayer import
#' @examples
#' # Read methylation data
#' methylation_gr <- readBedMethyl("path/to/methylation.bed", sample_name = "Sample1")
#'
#' # Annotate methylation data with GenBank features
#' annotated_methylation <- annotateMethylationWithGenBank(
#'   methylation_gr = methylation_gr,
#'   genbank_file = "path/to/genome.gbff",
#'   feature_types = c("genes", "cds", "ncRNA")
#' )
annotateMethylationWithGenBank <- function(methylation_gr, genbank_file, feature_types = c("genes", "cds", "tRNA")) {
  # Ensure necessary packages are installed
  if (!requireNamespace("GenomicFeatures", quietly = TRUE)) {
    stop("Package 'GenomicFeatures' is required but not installed.")
  }
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("Package 'GenomicRanges' is required but not installed.")
  }
  if (!requireNamespace("S4Vectors", quietly = TRUE)) {
    stop("Package 'S4Vectors' is required but not installed.")
  }
  if ("ncRNA" %in% feature_types && !requireNamespace("rtracklayer", quietly = TRUE)) {
    stop("Package 'rtracklayer' is required for 'ncRNA' features but not installed.")
  }

  # Read the GenBank file and create a TxDb object
  txdb <- GenomicFeatures::makeTxDbFromGFF(genbank_file)

  # Extract assembly information
  assembly_info <- seqinfo(txdb)

  # Update methylation_gr with assembly information
  seqlevels(methylation_gr, pruning.mode = "coarse") <- seqlevels(assembly_info)
  seqinfo(methylation_gr) <- assembly_info

  # Prepare features list
  features_list <- list()

  if ("genes" %in% feature_types) {
    genes_gr <- GenomicFeatures::genes(txdb)
    features_list$genes <- genes_gr
  }
  if ("transcripts" %in% feature_types) {
    transcripts_gr <- GenomicFeatures::transcripts(txdb)
    features_list$transcripts <- transcripts_gr
  }
  if ("cds" %in% feature_types) {
    cds_gr <- GenomicFeatures::cds(txdb)
    features_list$cds <- cds_gr
  }
  if ("exons" %in% feature_types) {
    exons_gr <- GenomicFeatures::exons(txdb)
    features_list$exons <- exons_gr
  }
  # For ncRNA, we need to get the non-coding RNA features
  if ("ncRNA" %in% feature_types) {
    # Import the GenBank file using rtracklayer
    genbank_gr <- rtracklayer::import(genbank_file, format = "genbank")
    # Filter for ncRNA features
    ncRNA_gr <- genbank_gr[genbank_gr$type == "ncRNA" | genbank_gr$type == "rRNA" | genbank_gr$type == "tRNA"]
    features_list$ncRNA <- ncRNA_gr
  }

  # Initialize metadata columns if not present
  if (is.null(mcols(methylation_gr)$feature_type)) {
    mcols(methylation_gr)$feature_type <- NA_character_
  }
  if (is.null(mcols(methylation_gr)$feature_id)) {
    mcols(methylation_gr)$feature_id <- NA_character_
  }
  if (is.null(mcols(methylation_gr)$feature_name)) {
    mcols(methylation_gr)$feature_name <- NA_character_
  }

  # Annotate methylation_gr with features
  for (feature_type in names(features_list)) {
    features_gr <- features_list[[feature_type]]
    if (is.null(features_gr)) next

    # Find overlaps
    overlaps <- GenomicRanges::findOverlaps(methylation_gr, features_gr)

    # Update metadata
    query_idx <- S4Vectors::queryHits(overlaps)
    subject_idx <- S4Vectors::subjectHits(overlaps)

    # Update feature_type
    mcols(methylation_gr)$feature_type[query_idx] <- feature_type
    # Update feature_id and feature_name
    feature_ids <- mcols(features_gr)$gene_id
    if (is.null(feature_ids)) {
      feature_ids <- mcols(features_gr)$ID  # Try 'ID' column
    }
    feature_names <- mcols(features_gr)$gene_name
    if (is.null(feature_names)) {
      feature_names <- mcols(features_gr)$Name  # Try 'Name' column
    }

    mcols(methylation_gr)$feature_id[query_idx] <- as.character(feature_ids[subject_idx])
    mcols(methylation_gr)$feature_name[query_idx] <- as.character(feature_names[subject_idx])
  }

  return(methylation_gr)
}

