buildMicrobeMethyl <- function(bed_methyl_path,
                               gff_path,
                               genome = "eco",
                               circular = TRUE) {
  # Load the genome annotation from GFF file
  gff <- readGFF(gff_path, genome = genome, circular = circular)

  # Reduce GFF to minimal metadata for now
  gff <- slimGFF(gff)

  # Load the methylated sites from BED file
  bed_df <- readBedMethyl(bed_methyl_path)
  genome(bed_df) <- genome
  isCircular(bed_df) <- circular

}


overlaps <- findOverlaps(bed_df, gff)
bed_hits <- queryHits(overlaps)
gff_hits <- subjectHits(overlaps)

