# Refine Spectral Entropy Using Temporal Template

Refines spectral entropy measurements by creating a temporal template
based on pitch goodness. The template is created using either quantile
thresholding or Hidden Markov Model (HMM) segmentation of the pitch
goodness values from a reference label.

## Usage

``` r
refine_sh(
  x,
  segment_type = c("motifs", "syllables", "segments"),
  reference_label,
  matrix = c("wiener", "shannon"),
  method = c("quantile", "hmm"),
  minimal_duration = 20,
  split_dips = TRUE,
  quantile_threshold = 0.5,
  random_seed = 222,
  plot = TRUE,
  plot_entropy_lim = NULL,
  color_palette = NULL,
  stats = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A Sap object containing spectral entropy and pitch goodness
  measurements

- segment_type:

  Character, type of segments to analyze: "motifs", "syllables", or
  "segments"

- reference_label:

  Character, label to use as reference for template creation

- matrix:

  Character, type of entropy matrix to refine: "wiener" or "shannon"

- method:

  Character, method for template creation: "quantile" or "hmm"

- minimal_duration:

  Numeric, minimum duration (in ms) for segments (default: 20)

- split_dips:

  Logical, whether to split segments at local minima (default: TRUE)

- quantile_threshold:

  Numeric, threshold for quantile method (default: 0.5)

- random_seed:

  Integer, seed for reproducibility in HMM (default: 222)

- plot:

  Logical, whether to plot results (default: TRUE)

- plot_entropy_lim:

  Numeric vector of length 2, limits for entropy plot

- color_palette:

  Function, custom color palette for plotting

- stats:

  Logical, whether to calculate segment statistics (default: TRUE)

- verbose:

  Logical, whether to print progress messages (default: TRUE)

- ...:

  Additional arguments passed to plotting functions

## Value

Invisibly returns the modified Sap object with refined spectral entropy
data

## Details

The function creates a temporal template based on pitch goodness values
from a reference label, which is then used to filter the spectral
entropy measurements. The template can be created using either a
quantile threshold or HMM segmentation. The filtered entropy values are
stored in the Sap object and can be visualized as a heatmap with the
template overlay.

## References

For HMM method: Visser, I., & Speekenbrink, M. (2010). depmixS4: An R
Package for Hidden Markov Models. Journal of Statistical Software,
36(7), 1-21.

## See also

[`refine_FF`](https://lxiao06.github.io/ASAP/reference/refine_FF.md) for
refining fundamental frequency
[`spectral_entropy`](https://lxiao06.github.io/ASAP/reference/spectral_entropy.md)
for calculating spectral entropy

## Examples

``` r
if (FALSE) { # \dontrun{
# Refine Wiener entropy using quantile method
sap <- refine_sh(sap,
                 segment_type = "motifs",
                 reference_label = "a",
                 matrix = "wiener",
                 method = "quantile",
                 plot = TRUE)

# Refine Shannon entropy using HMM method
sap <- refine_sh(sap,
                 segment_type = "motifs",
                 reference_label = "a",
                 matrix = "shannon",
                 method = "hmm",
                 plot = TRUE)
} # }
```
