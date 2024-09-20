#' Read a BedMethyl format file
#'
#' @param path a character string of the path to the file.
#' @param header a logical value indicating whether the file contains the names of the variables as its first line. Default is TRUE. Canonical bedMethyl files have no header.
#' @param add_header a logical value indicating whether to add a header to the data frame. Default is the opposite of `header`. This will add descriptive column names to the output. This is ignored if GRanges is TRUE.
#' @param GRanges a logical value indicating whether to return a GRanges object. Default is TRUE. If FALSE, a data frame will be returned.
#'
#' @return a GRanges object if GRanges is TRUE, otherwise a data frame.
#' @export
#'
#' @examples
readBedMethyl <- function(path,
                          header = TRUE,
                          add_header = !header,
                          GRanges = TRUE) {
  if (header == TRUE) {
    df <- read.table(path, header = TRUE, sep = "\t")
  } else {
    df <- read.table(path, header = FALSE, sep = "\t")
  }

  if (GRanges == TRUE) {
    df <- GRanges(
      seqnames = df[[1]],
      ranges = IRanges(start = df[[2]], end = df[[3]]),
      strand = df[[6]],
      coverage = df[[10]],
      percent_methylation = df[[11]]
    )
  }

  return(df)
}
