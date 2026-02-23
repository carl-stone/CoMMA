test_that("annotateMethylSites annotates overlapping genomic features", {
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

test_that("annotateTSS returns long output with relative positions", {
  methyl_df <- data.frame(Position = c(100L, 205L, 500L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = c("Transcription-Units", "Transcription-Units"),
    Strand = c("+", "-"),
    Left = c(90L, 170L),
    Right = c(130L, 210L),
    stringsAsFactors = FALSE
  )

  out <- annotateTSS(methyl_df, meta_df, location = "Position", size = 40, long = TRUE)

  expect_true(all(c("TSS_strand", "RelPos") %in% names(out)))
  expect_true(any(out$Position == 100L & out$RelPos == 10))
  expect_true(any(out$Position == 205L & out$RelPos == 5))
  expect_false(any(out$Position == 500L))
})

test_that("annotateTSS can return wide output with NoTSS marker", {
  methyl_df <- data.frame(Position = c(100L, 500L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = "Transcription-Units",
    Strand = "+",
    Left = 90L,
    Right = 130L,
    stringsAsFactors = FALSE
  )

  out <- annotateTSS(methyl_df, meta_df, location = "Position", size = 20, long = FALSE)

  expect_true(any(grepl("^RelPos_", names(out))))
  expect_equal(out$NoTSS[out$Position == 500L], "X")
})
