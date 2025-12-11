# Refine Fundamental Frequency Detection

Refines fundamental frequency detection by identifying regions of
reliable pitch tracking using various methods and applying temporal
constraints.

## Usage

``` r
refine_FF(
  x,
  segment_type = c("motifs", "syllables", "segments"),
  reference_label,
  method = c("quantile", "hmm"),
  minimal_duration = 20,
  split_dips = TRUE,
  quantile_threshold = 0.5,
  hmm_trans_prob = 0.9,
  random_seed = 222,
  plot = TRUE,
  plot_freq_lim = NULL,
  color_palette = NULL,
  stats = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  A Sap object containing fundamental frequency and pitch goodness data

- segment_type:

  Character, type of segment to analyze: "motifs", "syllables", or
  "segments"

- reference_label:

  Character, label to use as reference for template creation

- method:

  Character, method for template creation: "quantile" or "hmm"

- minimal_duration:

  Numeric, minimum duration (in ms) for valid segments

- split_dips:

  Logical, whether to split segments at significant dips in goodness

- quantile_threshold:

  Numeric, threshold for quantile method (0-1). Controls sensitivity of
  pitch detection - higher values (e.g., 0.7) retain only high-quality
  regions, lower values (e.g., 0.3) include more regions but may
  introduce noise

- hmm_trans_prob:

  Numeric, transition probability for HMM method (0-1). Controls state
  persistence - higher values (e.g., 0.95) create fewer, longer
  segments; lower values (e.g., 0.8) create more, shorter segments

- random_seed:

  Numeric, seed for reproducibility

- plot:

  Logical, whether to plot the results

- plot_freq_lim:

  Numeric vector of length 2, frequency limits for plotting

- color_palette:

  Function, color palette for plotting

- stats:

  Logical, whether to calculate segment statistics

- verbose:

  Logical, whether to print access information

- ...:

  Additional arguments passed to plotting functions

## Value

Invisibly returns the modified Sap object with refined fundamental
frequency data

## Details

The function uses either quantile-based thresholding or Hidden Markov
Models to identify regions of reliable pitch tracking. It can optionally
split segments at local minima in pitch goodness and applies minimum
duration constraints.

### Method Selection

- **Quantile method**: Simple thresholding based on pitch goodness
  values. Works well for clean recordings with good signal-to-noise
  ratio.

- **HMM method**: Uses a Hidden Markov Model to identify segments. More
  robust to noise and can better detect natural boundaries in the
  signal.

### Parameter Selection Guidelines

#### Quantile Threshold (quantile_threshold)

- **0.3-0.4**: Liberal threshold that includes most pitch-tracked
  regions. Use for high-quality recordings where you want to maximize
  data retention.

- **0.5-0.6**: Balanced threshold (default). Works well for most
  recordings.

- **0.7-0.8**: Conservative threshold that only includes regions with
  very reliable pitch tracking. Use for noisy recordings where precision
  is more important than recall.

#### HMM Transition Probability (hmm_trans_prob)

- **0.8-0.85**: Creates more responsive segmentation with shorter
  segments. Use when analyzing rapid vocalizations with frequent
  transitions.

- **0.9**: Balanced setting (default). Works well for most
  vocalizations.

- **0.95-0.98**: Creates more stable segmentation with fewer, longer
  segments. Use for sustained vocalizations or when you want to minimize
  over-segmentation.

## Note

The HMM method requires the 'depmixS4' package to be installed.

## Examples

``` r
if (FALSE) { # \dontrun{
# Default quantile method with balanced threshold
sap <- refine_FF(sap,
                  reference_label = "BL",
                  method = "quantile",
                  quantile_threshold = 0.5,
                  plot = TRUE)

# Conservative quantile threshold for noisy recordings
sap <- refine_FF(sap,
                  reference_label = "BL",
                  method = "quantile",
                  quantile_threshold = 0.7)

# HMM method with default transition probability
sap <- refine_FF(sap,
                  reference_label = "BL",
                  method = "hmm",
                  hmm_trans_prob = 0.9)

# HMM method for sustained vocalizations
sap <- refine_FF(sap,
                  reference_label = "BL",
                  method = "hmm",
                  hmm_trans_prob = 0.95)
} # }

```
