test_that("annotateMethylSites legacy wrapper preserves expected output contract", {
  methyl_df <- data.frame(Position = c(100L, 300L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = c("Gene", "Promoter"),
    Site = c("geneA", "proA"),
    Left = c(90L, 95L),
    Right = c(200L, 110L),
    stringsAsFactors = FALSE
  )

  out <- annotateMethylSites(methyl_df, meta_df, location = "Position")

  expect_equal(out$Gene[out$Position == 100L], "geneA")
  expect_equal(out$Promoter[out$Position == 100L], "proA")
  expect_equal(out$No_Feature[out$Position == 300L], "1")
})

test_that("legacy soft-deprecated helpers emit lifecycle warnings", {
  methyl_df <- data.frame(Position = c(100L, 200L), methyl = c(0.1, 0.2), cov = c(10, 15))
  tss_meta <- data.frame(
    Type = "Transcription-Units",
    Strand = "+",
    Left = 95L,
    Right = 110L,
    stringsAsFactors = FALSE
  )

  expect_warning(
    annotateTSS(methyl_df, tss_meta, location = "Position", size = 20, long = FALSE),
    "soft-deprecated"
  )

  expect_warning(
    methylRollingMean(methyl_df, position_col = "Position", methyl_col = "methyl", w_size = 100, genome_size = 1000),
    "soft-deprecated"
  )

  expect_warning(
    calculateMethylSiteDepth(methyl_df, position_col = "Position", cov_col = "cov", w_size = 100),
    "soft-deprecated"
  )

  var_dataset <- data.frame(
    Coverage_Sample = c(5L, 5L, 6L, 6L),
    Percent_Methyl_Sample = c(40, 45, 60, 55),
    Ancestor_Mean = c(42, 42, 58, 58)
  )

  expect_warning(
    varByCoverage(var_dataset),
    "soft-deprecated"
  )
})
