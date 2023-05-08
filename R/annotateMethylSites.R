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
#' df <- methyl_df
#' meta <- genome_sites
#' df_annotated <- annotateMethylSites(df, meta, location='Position')
#'
#' @seealso
#' \code{\link{findOverlaps}}, \code{\link{GenomicRanges-class}}
#'
#' @importFrom GenomicRanges findOverlaps
#'
#' @export
annotateMethylSites_dep <- function(methyl_df, meta_df, location) {
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[meta_df$Left <= position &
                     meta_df$Right >= position, ]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'No_Feature'] <- '1'
      next
    }
    sites_at_position <- meta_df[meta_df$Left <= position &
                                   meta_df$Right >= position, ]
    for (i in 1:nrow(sites_at_position)) {
      methyl_df[methyl_df[[location]] == position,
                toString(sites_at_position[i, 'Type'])] <-
        sites_at_position[i, 'Site']
    }
  }
  return(methyl_df)
}
