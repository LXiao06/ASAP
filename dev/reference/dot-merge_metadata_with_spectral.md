# Merge Segment Metadata with Spectral Features

Internal helper to merge exported segment metadata with newly extracted
spectral features before writing to a CSV.

## Usage

``` r
.merge_metadata_with_spectral(
  metadata_df,
  spectral_df,
  output_file,
  time_digits = 6,
  verbose = TRUE
)
```
