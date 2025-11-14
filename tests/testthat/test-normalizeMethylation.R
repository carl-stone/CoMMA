context("normalizeMethylation_dep")

test_that("normalizeMethylation_dep works without modelr on search path", {
  if ("package:modelr" %in% search()) {
    detach("package:modelr", unload = TRUE, character.only = TRUE)
  }
  data("all_samples", package = "CoMMA")
  expect_no_error(
    normalizeMethylation_dep(
      df = all_samples,
      plots = FALSE,
      rescale = FALSE
    )
  )
})
