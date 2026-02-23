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

test_that("annotateMethylSites supports non-Position location column names", {
  methyl_df <- data.frame(Loc = c(100L, 300L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = c("Gene"),
    Site = c("geneA"),
    Left = c(90L),
    Right = c(200L),
    stringsAsFactors = FALSE
  )

  out <- annotateMethylSites(methyl_df, meta_df, location = "Loc")

  expect_equal(out$Gene[out$Loc == 100L], "geneA")
  expect_equal(out$No_Feature[out$Loc == 300L], "1")
})

test_that("annotateMethylSites keeps duplicate locations isolated by row", {
  methyl_df <- data.frame(
    Position = c(100L, 100L, 350L),
    Sample = c("A", "B", "A"),
    Strand = c("+", "-", "+"),
    stringsAsFactors = FALSE
  )
  meta_df <- data.frame(
    Type = c("Gene"),
    Site = c("geneA"),
    Left = c(90L),
    Right = c(200L),
    stringsAsFactors = FALSE
  )

  out <- annotateMethylSites(methyl_df, meta_df, location = "Position")

  expect_equal(nrow(out), nrow(methyl_df))
  expect_equal(out$Sample, methyl_df$Sample)
  expect_equal(out$Strand, methyl_df$Strand)
  expect_equal(out$Gene[1:2], c("geneA", "geneA"))
  expect_equal(out$No_Feature[3], "1")
})

test_that("annotateMethylSites validates required columns", {
  methyl_df <- data.frame(Position = c(100L), stringsAsFactors = FALSE)
  meta_df <- data.frame(Type = "Gene", Left = 90L, Right = 200L, stringsAsFactors = FALSE)

  expect_error(
    annotateMethylSites(methyl_df, meta_df, location = "Position"),
    "missing required columns: Site"
  )

  expect_error(
    annotateMethylSites(methyl_df, meta_df = data.frame(), location = "Position"),
    "missing required columns"
  )

  expect_error(
    annotateMethylSites(methyl_df, data.frame(Type = "Gene", Site = "x", Left = 1L, Right = 2L), location = "Missing"),
    "missing required location column"
  )
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

test_that("annotateTSS supports non-Position location column names", {
  methyl_df <- data.frame(Loc = c(100L, 205L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = c("Transcription-Units", "Transcription-Units"),
    Strand = c("+", "-"),
    Left = c(90L, 170L),
    Right = c(130L, 210L),
    stringsAsFactors = FALSE
  )

  out <- annotateTSS(methyl_df, meta_df, location = "Loc", size = 40, long = TRUE)

  expect_true("Loc" %in% names(out))
  expect_true(any(out$Loc == 100L & out$RelPos == 10))
  expect_true(any(out$Loc == 205L & out$RelPos == 5))
})

test_that("annotateTSS long output preserves duplicate locations across context columns", {
  methyl_df <- data.frame(
    Position = c(100L, 100L),
    Sample = c("A", "B"),
    Strand = c("+", "-"),
    stringsAsFactors = FALSE
  )
  meta_df <- data.frame(
    Type = "Transcription-Units",
    Strand = "+",
    Left = 90L,
    Right = 130L,
    stringsAsFactors = FALSE
  )

  out <- annotateTSS(methyl_df, meta_df, location = "Position", size = 20, long = TRUE)

  expect_equal(nrow(out), 2)
  expect_setequal(out$Sample, c("A", "B"))
  expect_true(all(out$RelPos == 10))
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

test_that("annotateTSS long output handles empty matches", {
  methyl_df <- data.frame(Position = c(100L, 500L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = "Transcription-Units",
    Strand = "+",
    Left = 1000L,
    Right = 1100L,
    stringsAsFactors = FALSE
  )

  out <- annotateTSS(methyl_df, meta_df, location = "Position", size = 20, long = TRUE)

  expect_equal(nrow(out), 0)
  expect_true(all(c("Position", "NoTSS", "TSS_strand", "RelPos") %in% names(out)))
})

test_that("annotateTSS validates required columns", {
  methyl_df <- data.frame(Position = c(100L), stringsAsFactors = FALSE)
  meta_df <- data.frame(Type = "Transcription-Units", Left = 90L, Right = 130L)

  expect_error(
    annotateTSS(methyl_df, meta_df, location = "Position", size = 20),
    "missing required columns: Strand"
  )

  expect_error(
    annotateTSS(methyl_df, data.frame(Type = "Transcription-Units", Strand = "+", Left = 90L, Right = 130L), location = "Missing", size = 20),
    "missing required location column"
  )
})

test_that("annotation helpers run without attaching tidyverse", {
  expect_false("package:dplyr" %in% search())
  expect_false("package:tidyr" %in% search())
  expect_false("package:tibble" %in% search())

  methyl_df <- data.frame(Loc = c(100L, 500L), stringsAsFactors = FALSE)
  meta_df <- data.frame(
    Type = "Transcription-Units",
    Strand = "+",
    Left = 90L,
    Right = 130L,
    stringsAsFactors = FALSE
  )

  out_tss <- annotateTSS(methyl_df, meta_df, location = "Loc", size = 20, long = FALSE)
  expect_s3_class(out_tss, "data.frame")

  out_ann <- annotateMethylSites(
    methyl_df,
    meta_df = data.frame(Type = "Gene", Site = "geneA", Left = 90L, Right = 130L, stringsAsFactors = FALSE),
    location = "Loc"
  )
  expect_s3_class(out_ann, "data.frame")
})
