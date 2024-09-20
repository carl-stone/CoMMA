#' Create a SummarizedExperiment object with annotated methylation data
#'
#' @param bed_methyl_path a character string of the path to the bedMethyl file.
#' @param gff_path a character string of the path to the GFF file.
#' @param genome a character string specifying the genome name.
#' @param circular a logical value indicating whether the genome is circular. Default is TRUE.
#' @param metadata a list of metadata for the bedMethyl data, such as sample metadata.
#'
#' @return a SummarizedExperiment object with annotated methylation data.
#' @export buildMicrobeMethyl
#'
#' @examples
buildMicrobeMethyl <- function(bed_methyl_path_list,
                               gff_path,
                               genome = "eco",
                               circular = TRUE) {

  # Load the genome annotation from GFF file
  gff <- readGFF(gff_path, genome = genome, circular = circular)

  # Keep only useful genome annotations
  gff <- fatGFF(gff)

  # Load the methylated sites from bedMethyl files
  bed_df <- lapply(bed_methyl_path_list, readBedMethyl) |>
    GRangesList() |>
    as.data.frame(use.outer.mcols = TRUE)

  names(bed_df)[1:2] <- c("sample", "sample_name")

  # rowData
  row_data <- data.frame(
    start = bed_df[bed_df$sample == 1, "start"],
    end = bed_df[bed_df$sample == 1, "end"],
    strand = bed_df[bed_df$sample == 1, "strand"],
    seqnames = bed_df[bed_df$sample == 1, "seqnames"]
  ) |>
    as_granges()

  row_data <- plyranges::join_overlap_left(row_data, gff)
  # Some problem with fatGFF that I don't have with slimGFF
  row_data <- split(row_data, start(row_data))

  # colData
  col_data <- data.frame(
    sample = unique(bed_df$sample),
    sample_name = unique(bed_df$sample_name)
  )

  rownames(col_data) <- col_data$sample_name

  # get assay data
  assay_names <- c("coverage", "percent_methylation")
  assay_data <- bed_df[,c("sample", "start", assay_names)]

  assay_data <- assay_data |>
    pivot_wider(
      names_from = sample,
      values_from = c("coverage", "percent_methylation")
    )

  cov_data <- assay_data |>
    select(start, starts_with("coverage")) |>
    rename_with(~str_remove(., "coverage_"), starts_with("coverage")) |>
    column_to_rownames("start")

  methyl_data <- assay_data |>
    select(start, starts_with("percent_methylation")) |>
    rename_with(~str_remove(., "percent_methylation_"), starts_with("percent_methylation")) |>
    column_to_rownames("start")

  colnames(methyl_data) <- col_data$sample_name
  colnames(cov_data) <- col_data$sample_name

  # Create the RangedSummarizedExperiment object
  mm_data <- SummarizedExperiment(
    assays = list(
      coverage = cov_data,
      B = methyl_data
    ),
    rowRanges = row_data,
    colData = col_data
    )

  return(mm_data)
}
