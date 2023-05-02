setGeneric("annotateMethylSites", function(x, bed_file, location) {
  standardGeneric("annotateMethylSites")
})

setMethod("annotateMethylSites", "data.frame", function(x, bed_file, location) {
  # Read the BED file into a data frame
  meta_df <- rtracklayer::import(bed_file)
  colnames(meta_df) <- c("Chr", "Left", "Right", "Type", "Site", "Strand")

  # Convert the x and meta_df into GRanges objects
  methyl_ranges <- GRanges(seqnames = x$Chr,
                           ranges = IRanges(start = x$Start, end = x$End))
  meta_ranges <- GRanges(seqnames = meta_df$Chr,
                         ranges = IRanges(start = meta_df$Left, end = meta_df$Right),
                         Type = meta_df$Type, Site = meta_df$Site)

  # Find overlaps between methylated sites and genomic features
  overlaps <- findOverlaps(methyl_ranges, meta_ranges)

  # Annotate the methylated sites with overlapping features
  for (i in 1:length(methyl_ranges)) {
    hits <- subjectHits(overlaps[queryHits(overlaps) == i])
    if (length(hits) == 0) {
      x[i, 'No_Feature'] <- '1'
    } else {
      for (hit in hits) {
        feature_type <- as.character(meta_ranges[hit]$Type)
        feature_site <- as.character(meta_ranges[hit]$Site)
        x[i, feature_type] <- feature_site
      }
    }
  }

  return(x)
})

setMethod("annotateMethylSites", "microbeMethyl", function(x, bed_file, location) {
  # Assuming the data frame is stored in the 'data' slot of the microbeMethyl object
  annotated_data <- annotateMethylSites(data@x, bed_file, location)
  data@x <- annotated_data
  return(x)
})

setMethod("annotateMethylSites", "microbeMethylExperiment", function(x, bed_file, location) {
  # Assuming the list of microbeMethyl objects is stored in the 'samples' slot of the microbeMethylExperiment object
  for (i in seq_along(samples@x)) {
    samples@x[[i]] <- annotateMethylSites(samples@x[[i]], bed_file, location)
  }
  return(x)
})
