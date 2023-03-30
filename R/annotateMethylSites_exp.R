#' Annotate Methylated Sites
#'
#' This function takes as input a data frame of methylated sites and a metadata data frame, and annotates each methylated site with its corresponding feature(s) based on the given location. If a site does not overlap with any feature in the metadata, it is labeled as "No_Feature".
#'
#' @param methyl_df A data frame containing information on methylated sites, with columns "Chr", "Start", and "End" indicating the chromosome and genomic coordinates of each site.
#' @param meta_df A data frame containing metadata on genomic features, with columns "Chr", "Left", "Right", "Type", and "Site". "Left" and "Right" specify the genomic coordinates of the feature, "Type" specifies the type of feature (e.g. exon, intron, promoter), and "Site" specifies a unique identifier for the feature.
#' @param location A character string indicating which column in methyl_df should be used to determine the location of each methylated site (e.g. "Start" or "End").
#'
#' @return A data frame containing the same columns as methyl_df, with additional columns added for each feature type present in meta_df. Each row represents a methylated site, and the value in each feature column indicates the unique identifier(s) of the feature(s) that the site overlaps with. If a site does not overlap with any feature, the "No_Feature" column will be set to 1.
#'
#' @examples
#' # Load example data
#' data(methyl_df)
#' data(meta_df)
#'
#' # Annotate methylated sites with genomic features
#' annotated_df <- annotateMethylSites(methyl_df, meta_df, "Start")
#'
#' @seealso
#' \code{\link{GenomicRanges-class}}
#'
#' @importFrom GenomicRanges GRanges
#'
#' @export
annotateMethylSites <- function(methyl_df, meta_df, location) {

  # Convert data frames to GRanges objects
  methyl_gr <- with(methyl_df, GRanges(Chromosome, IRanges(get(location), get(location))))
  meta_gr <- with(meta_df, GRanges(Chromosome, IRanges(Left, Right), Site = Site, Type = Type))

  # Find overlaps between methylated sites and genomic features
  olaps <- findOverlaps(methyl_gr, meta_gr)

  # Create a matrix of feature IDs for each methylated site
  mat <- matrix(NA, nrow=length(methyl_gr), ncol=length(levels(meta_df$Type))+1,
                dimnames=list(NULL, c("No_Feature", levels(meta_df$Type))))
  for (i in 1:length(olaps)) {
    row <- queryHits(olaps[i])
    col <- as.character(meta_gr[subjectHits(olaps[i]), "Type"])
    mat[row, col] <- meta_gr[subjectHits(olaps[i]), "Site"]
  }
  # Label sites that did not overlap with any feature as "No_Feature"
  mat[is.na(mat)] <- "1"

  # Merge the annotated matrix with the original data frame
  methyl_df_annotated <- cbind(methyl_df, mat)
  return(methyl_df_annotated)
}
