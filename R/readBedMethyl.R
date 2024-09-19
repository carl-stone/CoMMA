readBedMethyl <- function(path,
                          header = TRUE,
                          GRanges = TRUE) {
  if (header == TRUE) {
    df <- read.table(path, header = TRUE, sep = "\t")
  } else {
    df <- read.table(path, header = FALSE, sep = "\t")
  }

  if (GRanges == TRUE) {
    df <- GRanges(seqnames = df[[1]],
                  ranges = IRanges(start = df[[2]],
                                   end = df[[3]]),
                  strand = df[[6]],
                  coverage = df[[10]],
                  percent_methylation = df[[11]])
  }
}
