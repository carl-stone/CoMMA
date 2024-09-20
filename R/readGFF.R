#' Read a GFF file
#'
#' @description
#' This function reads a GFF file and returns a GRanges object. The genome name and whether the genome is circular can be specified. By default, the function will remove the "region" features in the GFF file, which are typically used to define the boundaries of a chromosome or contig.
#'
#' @param gff_path a character string of the path to the GFF file.
#' @param genome a character string specifying the genome name.
#' @param circular a logical value indicating whether the genome is circular. Default is TRUE.
#' @param keep_chr a logical value indicating whether to keep the "region" features in the GFF file. Default is FALSE. "region" features are typically used to define the boundaries of a chromosome or contig.
#'
#' @return a GRanges object.
#' @export readGFF
#'
#' @examples
readGFF <- function(gff_path,
                    genome = "eco",
                    circular = TRUE,
                    keep_chr = FALSE) {
  gff <- rtracklayer::import(gff_path,
                             sequenceRegionsAsSeqinfo = TRUE)
  genome(gff) <- genome
  isCircular(gff) <- circular
  #seqlengths(gff) <- width(gff[mcols(gff)$type == "region"])

  if (!keep_chr) {
    gff <- gff[mcols(gff)$type != "region"]
  }

  return(gff)
}
