library(testthat)
library(CoMMA)


test_that("methylRaw objects keep vector lengths in sync", {
  sample_df <- data.frame(
    Position = c(100L, 200L, 300L, 400L),
    strand = c("+", "-", "+", "-"),
    cov_sample1 = c(10L, 20L, 15L, 25L),
    methyl_sample1 = c(0.5, 0.2, 0.8, 0.4),
    cov_sample2 = c(5L, 10L, 8L, 6L),
    methyl_sample2 = c(0.1, 0.9, 0.3, 0.7)
  )

  methyl_list <- buildMethylRawList(sample_df, strand_col = "strand")

  expect_s3_class(methyl_list, "methylRawList")
  expect_length(methyl_list, 2)

  for (i in seq_along(methyl_list)) {
    mr <- methyl_list[[i]]
    chr_len <- length(mr@.Data[[1]])
    start_len <- length(mr@.Data[[2]])
    end_len <- length(mr@.Data[[3]])
    strand_len <- length(mr@.Data[[4]])
    cov_len <- length(mr@.Data[[5]])

    expect_equal(chr_len, start_len)
    expect_equal(chr_len, end_len)
    expect_equal(chr_len, strand_len)
    expect_equal(chr_len, cov_len)
    expect_equal(mr@.Data[[4]], sample_df$strand)
  }
})
