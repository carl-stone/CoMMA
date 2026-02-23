# Internal helpers for legacy compatibility shims --------------------------------

.comma_lifecycle_state <- new.env(parent = emptyenv())

legacy_soft_deprecate <- function(fn, replacement) {
  if (!exists(fn, envir = .comma_lifecycle_state, inherits = FALSE)) {
    warning(
      sprintf(
        "`%s()` is soft-deprecated and retained for compatibility; migrate to `%s()`.",
        fn,
        replacement
      ),
      call. = FALSE
    )
    assign(fn, TRUE, envir = .comma_lifecycle_state)
  }

  invisible(NULL)
}
