library(testthat)
library(GenomicRanges)

test_that("parseGenBankFeatures flattens qualifiers without plyr", {
  features_lines <- c(
    "FEATURES             Location/Qualifiers",
    "     gene            1..100",
    "                     /gene=\"geneA\"",
    "                     /note=\"example\""
  )

  result <- parseGenBankFeatures(features_lines, seqname = "NC_000001")

  expect_s4_class(result, "GRanges")
  expect_equal(mcols(result)$gene, "geneA")
  expect_equal(mcols(result)$note, "example")
})
