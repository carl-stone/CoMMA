readGFF <- function(path,
                    genome = "eco",
                    circular = TRUE,
                    keep_chr = FALSE) {
  gff <- rtracklayer::import.gff(path,
                                 sequenceRegionsAsSeqinfo = TRUE,
                                 genome = genome)
  isCircular(gff) <- circular

  if (!keep_chr) {
    gff <- gff[-1]
    gff <- gff[, c("type", "ID", "Name", "gbkey", "gene", "gene_biotype",
                   "gene_synonym", "locus_tag", "Parent", "orig_transcript_id",
                   "product", "protein_id", "mobile_element_type",
                   "pseudo", "orig_protein_id", "exception", "transl_except")]
  }
  return(gff)
}

