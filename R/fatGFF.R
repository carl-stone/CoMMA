#' @title Select maximum useful genome annotation columns from a GFF GRanges object.
#'
#' @param gff A GRanges object with genome annotation.
#'
#' @return A GRanges object with range positions, type, gene, locus tag, parent, product, original transcript ID, protein ID, translation table, and mobile element type.
#' @export fatGFF
#'
#' @examples
fatGFF <- function(gff) {
  gff_meta_cols <- c(
    "type",
    "gene",
    "locus_tag",
    #"Parent",
    "product",
    "protein_id",
    "mobile_element_type"
  )
  # For now including Parent throws an error, I think because it's the only list column
  mcols(gff) <- mcols(gff)[, gff_meta_cols]
  return(gff)
}
