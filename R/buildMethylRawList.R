#' Build methylRawList object
#'
#' @param df data frame containing columns starting with "cov_" and "methyl_"
#'  for each sample (i.e. "cov_sample1")
#' @param assembly a string, chromosome assembly, only relevant if there are
#'  multiple chroms otherwise any string
#' @param context a string, sequence context around methylated base
#' @param resolution a string, either "base" or "region"
#'
#' @return a MethylRawList object ({methylKit} S4 class)
#' @export
#'
#' @examples
#' # Create a sample data frame with required columns
#' sample_data <- data.frame(
#'   Position = c(100, 200, 300),
#'   cov_sample1 = c(100, 150, 80),
#'   methyl_sample1 = c(0.9, 0.8, 0.95),
#'   cov_sample2 = c(40, 50, 60),
#'   methyl_sample2 = c(0.6, 0.5, 0.4)
#' )
#'
#' # Build a methylRawList object using the sample data frame
#' my_methylRawList <- buildMethylRawList(
#'   df = sample_data,
#'   assembly = "NC_000913.3",
#'   context = "GATC",
#'   resolution = "base"
#' )
#'
#' # Inspect the resulting methylRawList object
#' my_methylRawList

buildMethylRawList <- function(df,
                               assembly = "NC_000913.3",
                               context = "GATC",
                               resolution = "base") {
  coverage_df <- df %>%
    dplyr::select(starts_with("cov"))
  methyl_df <- df %>%
    dplyr::select(starts_with("methyl"))
  n_samples <- ncol(coverage_df)
  sample_ids <- strsplit(colnames(coverage_df), split = "_")
  sample_ids <- matrix(unlist(sample_ids), ncol = n_samples)[2,]
  methylRawList_obj <- new("methylRawList")

  for (i in 1:n_samples) {
    methylRaw_obj <- new("methylRaw")
    methylRaw_obj@sample.id = sample_ids[i]
    methylRaw_obj@assembly = assembly
    methylRaw_obj@context = context
    methylRaw_obj@resolution = resolution
    methylRaw_obj@.Data = list(rep(assembly, nrow(df)),
                               df$Position,
                               df$Position,
                               c(rep(c("+", "-"), floor(nrow(df / 2))), "+"),
                               coverage_df[,i],
                               methyl_df[,i] * coverage_df[,i],
                               (1 - methyl_df[,i]) * coverage_df[,i])
    methylRaw_obj@names = c("chr", "start", "end", "strand", "coverage", "numCs", "numTs")
    methylRaw_obj@row.names = rownames(df)
    methylRawList_obj[[i]] <- methylRaw_obj
  }

  return(methylRawList_obj)
}
