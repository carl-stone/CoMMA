test_that("methylRollingMedian exact mode returns expected schema", {
  df <- data.frame(pos = c(1L, 3L, 5L), beta = c(0.2, 0.8, 0.4))

  out <- methylRollingMedian(df, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6, method = "exact")

  expect_named(out, c("position", "methyl", "med_methyl"))
  expect_equal(nrow(out), 6)
  expect_true(all(out$position == 1:6))
})

test_that("methylRollingMedian fast mode returns site-level output", {
  df <- data.frame(pos = c(1L, 3L, 5L), beta = c(0.2, 0.8, 0.4))

  out <- methylRollingMedian(df, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6, method = "fast")

  expect_named(out, c("position", "mean_methyl"))
  expect_equal(nrow(out), 3)
})

test_that("writeBED writes track line and BED rows", {
  site_table <- data.frame(Position = c(10L, 20L), beta = c(0.5, 0.95), Strand = c("+", "-"))
  out_path <- tempfile(fileext = ".bed")

  bed <- writeBED(site_table, out_path)

  expect_true(file.exists(out_path))
  expect_named(bed, c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "thickStart", "thickEnd", "itemRGB"))
  lines <- readLines(out_path)
  expect_match(lines[[1]], "^track name=methylation")
  expect_equal(length(lines), 3)
  expect_equal(bed$chromStart, c(9L, 19L))
  expect_equal(bed$chromEnd, c(10L, 20L))
  expect_equal(bed$thickStart, c(9L, 19L))
  expect_equal(bed$thickEnd, c(10L, 20L))

  bed_rows <- read.delim(
    out_path,
    sep = "\t",
    header = FALSE,
    skip = 1,
    stringsAsFactors = FALSE
  )
  expect_equal(bed_rows$V2, c(9L, 19L))
  expect_equal(bed_rows$V3, c(10L, 20L))
})

test_that("writeBED validates required columns", {
  bad <- data.frame(Position = 10L, Strand = "+")
  expect_error(writeBED(bad, tempfile(fileext = ".bed")), "Missing required columns")
})

test_that("internal helper .remove_mutated_sites supports numeric and data.frame inputs", {
  site_table <- data.frame(
    seqname = c("chr1", "chr1", "chr2"),
    pos = c(10L, 20L, 30L),
    strand = c("+", "+", "-"),
    mutated = c(FALSE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )

  out_default <- CoMMA:::.remove_mutated_sites(site_table, NULL)
  expect_false(any(out_default$mutated))

  out_numeric <- CoMMA:::.remove_mutated_sites(site_table, 10L)
  expect_false(any(out_numeric$pos == 10L))

  mutated_df <- data.frame(seqname = "chr2", pos = 30L)
  out_df <- CoMMA:::.remove_mutated_sites(site_table, mutated_df)
  expect_false(any(out_df$pos == 30L))
})


test_that("rolling helpers sort unsorted inputs before fast window calculations", {
  sorted_df <- data.frame(pos = c(1L, 3L, 5L), beta = c(0.2, 0.8, 0.4))
  unsorted_df <- sorted_df[c(3, 1, 2), ]

  out_sorted <- methylRollingMedian(sorted_df, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6, method = "fast")
  out_unsorted <- methylRollingMedian(unsorted_df, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6, method = "fast")

  expect_equal(out_unsorted, out_sorted)
})

test_that("rolling helpers validate required columns and numeric inputs", {
  good <- data.frame(pos = c(1L, 3L), beta = c(0.2, 0.8))

  expect_error(
    methylRollingMedian(good, position_col = "missing", methyl_col = "beta", w_size = 2, genome_size = 6),
    "Missing required columns"
  )
  expect_error(
    methylRollingMedian(data.frame(pos = c("a", "b"), beta = c(0.2, 0.8)), position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6),
    "position_col"
  )
  expect_error(
    methylRollingMedian(data.frame(pos = c(1L, 2L), beta = c("x", "y")), position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6),
    "methyl_col"
  )
})

test_that("rolling helpers validate window, genome size, and method", {
  df <- data.frame(pos = c(1L, 3L), beta = c(0.2, 0.8))

  expect_error(
    methylRollingMedian(df, position_col = "pos", methyl_col = "beta", w_size = 0, genome_size = 6),
    "w_size"
  )
  expect_error(
    methylRollingMedian(df, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 0),
    "genome_size"
  )
  expect_error(
    methylRollingMedian(df, position_col = "pos", methyl_col = "beta", w_size = 2, genome_size = 6, method = "slow"),
    "method"
  )
  expect_error(
    methylRollingMean(df, position_col = "pos", methyl_col = "beta", method = "exact"),
    "method"
  )
})
