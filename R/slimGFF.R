#' Select minimum genome annotation columns from a GFF GRanges object.
#'
#' @param gff A GRanges object with genome annotation.
#'
#' @return A GRanges object with range positions, type, gene, and locus tag.
#' @export slimGFF
#'
#' @examples
slimGFF <- function(gff) {
  gff_meta_cols <- c(
    "type",
    "gene",
    "locus_tag"
  )
  mcols(gff) <- mcols(gff)[, gff_meta_cols]
  return(gff)
}
