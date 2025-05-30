% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umap.R
\name{plot_single_umap}
\alias{plot_single_umap}
\title{Internal UMAP Plotting Function}
\usage{
plot_single_umap(
  data,
  dims = c("X1", "X2"),
  cols = NULL,
  pt.size = 0.5,
  stroke = 0.5,
  alpha = 1,
  highlight.alpha = NULL,
  label = FALSE,
  repel = FALSE,
  label.size = 4,
  highlight.by = NULL,
  highlight.value = NULL,
  cols.highlight = "#DE2D26",
  sizes.highlight = 1,
  background.value = NULL,
  na.value = "grey80"
)
}
\description{
Internal function for creating individual UMAP plots.
}
\keyword{internal}
