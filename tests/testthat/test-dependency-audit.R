test_that("static dependency audit reports no drift", {
  audit_script <- testthat::test_path("..", "..", "tools", "dependency_audit.R")
  source(audit_script, local = TRUE)

  audit <- comma_dependency_audit(repo_root = testthat::test_path("..", ".."))

  expect_length(audit$unused_imports, 0)
  expect_length(audit$undeclared_runtime_namespaces, 0)
  expect_length(audit$vignette_not_suggested, 0)

  expect_invisible(comma_assert_dependency_audit(audit))
})
