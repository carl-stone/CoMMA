#' Build methylRawList object
#'
#' @param df data frame containing columns starting with "cov_" and "methyl_"
#'  for each sample (i.e. "cov_sample1")
#' @param assembly a string, chromosome assembly, only relevant if there are
#'  multiple chroms otherwise any string
#' @param context a string, sequence context around methylated base
#' @param resolution a string, either "base" or "region"
#' @param strand_col optional column name in `df` that stores strand information
#'  to include in the resulting methylRaw objects
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
                               resolution = "base",
                               strand_col = NULL) {
  coverage_df <- df %>%
    dplyr::select(starts_with("cov"))
  methyl_df <- df %>%
    dplyr::select(starts_with("methyl"))
  n_samples <- ncol(coverage_df)
  sample_ids <- strsplit(colnames(coverage_df), split = "_")
  sample_ids <- matrix(unlist(sample_ids), ncol = n_samples)[2,]
  methylRawList_obj <- new("methylRawList")

  strand_values <- NULL
  if (!is.null(strand_col)) {
    if (!strand_col %in% names(df)) {
      stop("strand_col must match a column name in df")
    }
    strand_values <- df[[strand_col]]
  } else if ("strand" %in% names(df)) {
    strand_values <- df[["strand"]]
  }
  if (is.null(strand_values)) {
    strand_values <- rep("+", nrow(df))
  } else if (length(strand_values) != nrow(df)) {
    stop("strand column must have the same number of rows as df")
  }

  for (i in 1:n_samples) {
    methylRaw_obj <- new("methylRaw")
    methylRaw_obj@sample.id = sample_ids[i]
    methylRaw_obj@assembly = assembly
    methylRaw_obj@context = context
    methylRaw_obj@resolution = resolution
    methylRaw_obj@.Data = list(rep(assembly, nrow(df)),
                               df$Position,
                               df$Position,
                               strand_values,
                               coverage_df[,i],
                               methyl_df[,i] * coverage_df[,i],
                               (1 - methyl_df[,i]) * coverage_df[,i])
    methylRaw_obj@names = c("chr", "start", "end", "strand", "coverage", "numCs", "numTs")
    methylRaw_obj@row.names = rownames(df)
    methylRawList_obj[[i]] <- methylRaw_obj
  }

  return(methylRawList_obj)
}
