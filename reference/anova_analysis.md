# Perform ANOVA and Multiple Comparisons Analysis

Performs one-way ANOVA and Tukey's HSD test for multiple comparisons
across different segments. Provides both statistical results and
optional visualization.

## Usage

``` r
anova_analysis(stats_df, plot = TRUE)
```

## Arguments

- stats_df:

  A data frame containing columns:

  - segment_id: Numeric identifier for segments

  - label: Factor or character indicating groups to compare

  - mean: Numeric values for comparison

- plot:

  Logical, whether to create a boxplot visualization (default: TRUE)

## Value

A tibble containing ANOVA results with columns:

- segment_id: Segment identifier

- term: Source of variation (label or Residuals)

- df: Degrees of freedom

- sumsq: Sum of squares

- meansq: Mean squares

- statistic: F-statistic

- p.value: P-value

- significant: Logical indicating if p.value \< 0.05

## Details

The function performs two main analyses:

- One-way ANOVA for each segment

- Tukey's HSD test for multiple comparisons with adjusted p-values

The printed output includes:

- Tukey's HSD results with adjusted p-values

- Significance levels: \*\*\* (p\<0.001), \*\* (p\<0.01), \* (p\<0.05),
  ns (p\>=0.05)

- Optional boxplot visualization

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
results <- anova_results(stats_df)

# Without plot
results <- anova_results(stats_df, plot = FALSE)
} # }
```
