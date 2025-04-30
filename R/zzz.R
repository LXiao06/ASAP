.pkg_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  optional_pkgs <- c(
    "ggplot2", "tidyr", "patchwork", "pbapply",
    "pbmcapply", "parallel", "reticulate"
  )

  for (pkg in optional_pkgs) {
    .pkg_env[[paste0("has_", pkg)]] <- requireNamespace(pkg, quietly = TRUE)
  }
}
