buildMicrobeMethyl <- function(bedMethyl_path,
                               annotation_path) {
  # Read in the bedMethyl file
  bedMethyl <- readBedMethyl(bedMethyl_path)

  # Read in the annotation file
  annotation <- readGFF(annotation_path,
                        genome = "eco")

  # Build a GRanges object from bedMethyl
  mm <- GRanges(
    seqnames = unique(seqnames(annotation)),
    ranges = IRanges(start = bedMethyl$start,
                     end = bedMethyl$end),
    strand = bedMethyl$strand,
    coverage = bedMethyl$coverage,
    percent_methylation = bedMethyl$percentMethylation
  )

  seqinfo(mm) <- seqinfo(annotation)

  # Annotate the GRanges object
  o <- findOverlaps(annotation, mm)
  mcols(mm)$type <- NA
  mcols(mm)[queryHits(o)]$type <- mcols(annotation)$type[subjectHits(o)]


  return(mm)
}
