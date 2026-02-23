# Lightweight static dependency audit for CoMMA
#
# Usage from repo root:
#   source("tools/dependency_audit.R")
#   audit <- comma_dependency_audit()
#   print(audit)
#   comma_assert_dependency_audit(audit)

comma_read_desc_field <- function(field, desc_path = "DESCRIPTION") {
  desc <- as.list(utils::read.dcf(desc_path)[1, ])
  value <- desc[[field]]
  if (is.null(value) || !nzchar(value)) {
    return(character())
  }

  items <- strsplit(value, ",")[[1]]
  items <- trimws(items)
  items <- gsub("\\s*\\(.*?\\)", "", items)
  items <- items[nzchar(items)]
  unique(items)
}

comma_extract_namespace_packages <- function(files) {
  if (length(files) == 0) {
    return(character())
  }

  pattern <- "\\b([A-Za-z][A-Za-z0-9.]*):::{0,2}[A-Za-z][A-Za-z0-9._]*"
  refs <- unlist(lapply(files, function(path) {
    txt <- readLines(path, warn = FALSE)
    m <- gregexpr(pattern, txt, perl = TRUE)
    hits <- regmatches(txt, m)
    unlist(hits, use.names = FALSE)
  }), use.names = FALSE)

  if (length(refs) == 0) {
    return(character())
  }

  pkgs <- sub(":::{0,2}.*$", "", refs)
  unique(pkgs)
}

comma_importfrom_packages <- function(namespace_path = "NAMESPACE") {
  if (!file.exists(namespace_path)) {
    return(character())
  }

  ns <- readLines(namespace_path, warn = FALSE)
  import_lines <- grep("^importFrom\\(", ns, value = TRUE)
  if (length(import_lines) == 0) {
    return(character())
  }

  sub("^importFrom\\(([^,]+),.*$", "\\1", import_lines)
}

comma_dependency_audit <- function(repo_root = ".") {
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(repo_root)

  r_files <- Sys.glob("R/*.R")
  vignette_files <- Sys.glob("vignettes/*.Rmd")

  imports <- setdiff(comma_read_desc_field("Imports"), "R")
  suggests <- setdiff(comma_read_desc_field("Suggests"), "R")

  runtime_namespace_refs <- comma_extract_namespace_packages(r_files)
  vignette_namespace_refs <- comma_extract_namespace_packages(vignette_files)
  imported_by_namespace <- comma_importfrom_packages("NAMESPACE")

  used_imports <- sort(unique(c(runtime_namespace_refs, imported_by_namespace)))

  list(
    imports = sort(imports),
    suggests = sort(suggests),
    runtime_namespace_refs = sort(runtime_namespace_refs),
    vignette_namespace_refs = sort(vignette_namespace_refs),
    imported_by_namespace = sort(unique(imported_by_namespace)),
    unused_imports = sort(setdiff(imports, used_imports)),
    undeclared_runtime_namespaces = sort(setdiff(runtime_namespace_refs, c(imports, suggests))),
    vignette_not_suggested = sort(setdiff(vignette_namespace_refs, c(imports, suggests)))
  )
}

comma_assert_dependency_audit <- function(audit = comma_dependency_audit()) {
  issues <- character()

  if (length(audit$unused_imports) > 0) {
    issues <- c(
      issues,
      sprintf(
        "Unused Imports: %s",
        paste(audit$unused_imports, collapse = ", ")
      )
    )
  }

  if (length(audit$undeclared_runtime_namespaces) > 0) {
    issues <- c(
      issues,
      sprintf(
        "Runtime namespace refs missing from DESCRIPTION Imports/Suggests: %s",
        paste(audit$undeclared_runtime_namespaces, collapse = ", ")
      )
    )
  }

  if (length(audit$vignette_not_suggested) > 0) {
    issues <- c(
      issues,
      sprintf(
        "Vignette namespace refs missing from DESCRIPTION Imports/Suggests: %s",
        paste(audit$vignette_not_suggested, collapse = ", ")
      )
    )
  }

  if (length(issues) > 0) {
    stop(paste(issues, collapse = "\n"), call. = FALSE)
  }

  invisible(audit)
}
