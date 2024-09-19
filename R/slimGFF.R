slimGFF <- function(gff) {
  gff_meta_cols <- c(
    "type",
    "gene",
    "locus_tag"
  )
  mcols(gff) <- mcols(gff)[, gff_meta_cols]
  return(gff)
}
