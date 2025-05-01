#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom Rcpp sourceCpp
#' @useDynLib ASAP, .registration = TRUE
#' @importClassesFrom tuneR Wave
#' @importFrom tuneR readWave
#' @importFrom rlang `%||%` .data
#' @importFrom tools file_path_sans_ext
#' @importFrom grDevices colorRampPalette dev.cur dev.off png rgb
#' @importFrom utils head tail
#' @importFrom graphics abline axis box image layout legend lines mtext par
#'             plot.new plot.window points rect text title
#' @importFrom dplyr mutate select filter arrange left_join bind_rows %>%
#'             group_by n_distinct group_split case_when all_of
#' @importFrom stats TukeyHSD aggregate aov approx as.formula dist gaussian
#'             median na.omit prcomp quantile setNames sd time fft
#' @importFrom seewave spec inputw ftwindow sfm sh th meanspec afilter
# #' @importFrom ggplot2 ggplot aes geom_line geom_boxplot facet_wrap labs
# #'             theme_minimal scale_color_brewer scale_fill_brewer ggtitle
# #'             theme element_text stat_summary
## usethis namespace: end
NULL

# Define global variables used in NSE contexts
utils::globalVariables(c(
  "day_post_hatch",
  "dph",
  "UMAP1",
  "UMAP2",
  "start_time",
  "end_time",
  "label",
  "filename",
  "duration"
))
