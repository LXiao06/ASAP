# Filter Outlier Renditions from Trajectory Embeddings

Identifies and removes entire renditions (trials) from trajectory
embeddings when a substantial fraction of their time steps are
statistical outliers relative to the group distribution. Outlier
detection is performed per label, per time step, and per embedding
dimension using the IQR method.

## Usage

``` r
filter_trajectory_outliers(x, ...)

# Default S3 method
filter_trajectory_outliers(
  x,
  dims = c("PC1", "PC2"),
  iqr_multiplier = 4,
  min_outlier_fraction = 0.1,
  labels = NULL,
  verbose = TRUE,
  ...
)

# S3 method for class 'Sap'
filter_trajectory_outliers(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  dims = c("PC1", "PC2"),
  iqr_multiplier = 4,
  min_outlier_fraction = 0.1,
  labels = NULL,
  plot = FALSE,
  balanced = FALSE,
  ordered = FALSE,
  clusters = NULL,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  An object to filter: a trajectory embeddings data frame or SAP object

- ...:

  Additional arguments

- dims:

  Character vector of dimension columns to assess (default: c("PC1",
  "PC2"))

- iqr_multiplier:

  Multiplier applied to IQR to define outlier bounds (default: 4)

- min_outlier_fraction:

  Minimum fraction of time steps that must be outliers before the entire
  rendition is removed (default: 0.1, i.e., 10%)

- labels:

  Optional character vector to restrict filtering to specific labels

- verbose:

  Whether to print a removal summary (default: TRUE)

- segment_type:

  For SAP objects: Type of segments ('motifs', 'syllables', 'bouts',
  'segments')

- plot:

  Logical. If `TRUE`, plot a heatmap of the filtered trajectories
  (default: `FALSE`).

- balanced:

  Logical. Whether to balance segment numbers across groups for the
  heatmap (default: `FALSE`).

- ordered:

  Logical. Whether to order segments for the heatmap (default: `FALSE`).

- clusters:

  Integer vector of clusters to include for the heatmap (default:
  `NULL`).

## Value

For the default method: the filtered trajectory embeddings data frame
(same structure as input, with outlier renditions removed).

For SAP objects: the updated SAP object with filtered `traj.embeds`
stored back in `x$features[[feature_type]]$traj.embeds`.

## Details

The outlier detection algorithm works as follows:

- For each label separately, at each time step, compute Q1, Q3, and IQR
  across all renditions for each specified dimension.

- A time step is flagged as an outlier for a given rendition if its
  value in **any** dimension falls outside
  `[Q1 - iqr_multiplier * IQR, Q3 + iqr_multiplier * IQR]`.

- A rendition is removed entirely if the fraction of flagged time steps
  exceeds `min_outlier_fraction`.

Outlier bounds are computed within each label to avoid penalizing
conditions with naturally different mean trajectories.

## Examples

``` r
if (FALSE) { # \dontrun{
# Filter from data frame directly
traj_clean <- filter_trajectory_outliers(
  sap$features$motif$traj.embeds,
  dims = c("PC1", "PC2"),
  iqr_multiplier = 4,
  min_outlier_fraction = 0.1
)

# Filter from SAP object (updates traj.embeds in place)
sap <- filter_trajectory_outliers(sap)

# Stricter: remove renditions with > 5% outlier time points
sap <- filter_trajectory_outliers(sap, min_outlier_fraction = 0.05)

# Use filtered data directly with trajectory_dispersion
result <- trajectory_dispersion(traj_clean, dims = c("PC1", "PC2"))
} # }
```
