#' Write methylation calls as a BED9 file
#'
#' Convert a methylation site table into a BED9 table with an optional UCSC/IGV
#' track header and write it to disk.
#'
#' ## Input schema
#' `site_table` must contain:
#' * `position_col`: 1-based genomic position (`integer`/`numeric`)
#' * `methyl_col`: methylation fraction in `[0, 1]` (`numeric`)
#' * `strand_col`: strand values (`"+"`, `"-"`, or `"."`)
#'
#' ## Return schema
#' Returns (invisibly) a `data.frame` with BED9 columns:
#' `chrom`, `chromStart`, `chromEnd`, `name`, `score`, `strand`,
#' `thickStart`, `thickEnd`, `itemRGB`.
#'
#' Input positions are validated as 1-based coordinates and then converted to
#' BED coordinates (0-based, half-open):
#' * `chromStart = position - 1`
#' * `chromEnd = position`
#'
#' @param site_table Data frame of methylation calls.
#' @param output_path File path to write.
#' @param position_col Column name containing 1-based site positions.
#' @param methyl_col Column name containing methylation beta values in `[0, 1]`.
#' @param strand_col Column name containing strand values.
#' @param chrom Chromosome/contig label to write in BED `chrom` column.
#' @param track_name Track name used in the optional header line.
#' @param description Track description used in the optional header line.
#' @param include_track_line If `TRUE`, write a UCSC-style track header before
#'   BED rows.
#'
#' @return Invisibly returns the BED data.frame written to `output_path`.
#' @export
#'
#' @examples
#' df <- data.frame(Position = c(10L, 25L), beta = c(0.2, 0.9), Strand = c("+", "-"))
#' out <- tempfile(fileext = ".bed")
#' writeBED(df, out, position_col = "Position", methyl_col = "beta", strand_col = "Strand")
writeBED <- function(site_table,
                     output_path,
                     position_col = "Position",
                     methyl_col = "beta",
                     strand_col = "Strand",
                     chrom = "U00096.3",
                     track_name = "methylation",
                     description = "Methylation",
                     include_track_line = TRUE) {
  if (!is.data.frame(site_table)) {
    stop("`site_table` must be a data.frame.", call. = FALSE)
  }
  required <- c(position_col, methyl_col, strand_col)
  missing <- setdiff(required, names(site_table))
  if (length(missing) > 0) {
    stop(sprintf("Missing required columns: %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  position <- site_table[[position_col]]
  methyl <- site_table[[methyl_col]]
  strand <- as.character(site_table[[strand_col]])

  if (!is.numeric(position) || any(is.na(position)) || any(position <= 0) || any(position != floor(position))) {
    stop("`position_col` must contain positive integer coordinates.", call. = FALSE)
  }
  if (!is.numeric(methyl) || any(is.na(methyl)) || any(methyl < 0 | methyl > 1)) {
    stop("`methyl_col` must contain numeric methylation fractions between 0 and 1.", call. = FALSE)
  }
  if (any(!(strand %in% c("+", "-", ".")))) {
    stop("`strand_col` must contain only '+', '-', or '.'.", call. = FALSE)
  }

  score <- round(methyl * 100)
  item_rgb <- rep("250,0,0", length(score))
  item_rgb[score <= 95] <- "222,0,28"
  item_rgb[score <= 90] <- "167,0,85"
  item_rgb[score <= 80] <- "125,0,128"
  item_rgb[score <= 60] <- "0,0,255"

  bed_file <- data.frame(
    chrom = rep(chrom, length(position)),
    chromStart = as.integer(position - 1L),
    chromEnd = as.integer(position),
    name = as.character(position),
    score = as.integer(score),
    strand = strand,
    thickStart = as.integer(position - 1L),
    thickEnd = as.integer(position),
    itemRGB = item_rgb,
    stringsAsFactors = FALSE
  )

  con <- file(output_path, open = "wt")
  on.exit(close(con), add = TRUE)

  if (isTRUE(include_track_line)) {
    cat(sprintf('track name=%s description="%s" itemRgb="On"\n', track_name, description), file = con)
  }

  utils::write.table(
    bed_file,
    file = con,
    sep = "\t",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  invisible(bed_file)
}
