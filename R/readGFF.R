readGFF <- function(gff_path,
                    genome = "eco",
                    circular = TRUE,
                    keep_chr = FALSE) {
  gff <- rtracklayer::import(gff_path)
  genome(gff) <- genome
  isCircular(gff) <- circular
  seqlengths(gff) <- width(gff[mcols(gff)$type == "region"])

  if (!keep_chr) {
    gff <- gff[mcols(gff)$type != "region"]
  }

  return(gff)
}
