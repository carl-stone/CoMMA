buildMicrobeMethyl <- function(bed_methyl_path,
                               gff_path,
                               genome = "eco",
                               circular = TRUE,
                               metadata) {
  # Load the genome annotation from GFF file
  gff <- readGFF(gff_path, genome = genome, circular = circular)

  # Reduce GFF to minimal metadata for now
  gff <- slimGFF(gff)

  # Load the methylated sites from BED file
  bed_df <- readBedMethyl(bed_methyl_path)
  genome(bed_df) <- genome
  isCircular(bed_df) <- circular
  metadata(bed_df) <- metadata

  # Add metadata to methylation sites using plyranges
  bed_df <- join_overlap_left(bed_df, gff)
}
