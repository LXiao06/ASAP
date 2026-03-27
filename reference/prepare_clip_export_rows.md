# Prepare Clip Data Frame for Export

Ensures required columns (`day_post_hatch`, `label`, `bird_id`) are
present, substituting sensible defaults for missing or NA values.

## Usage

``` r
prepare_clip_export_rows(clips)
```

## Arguments

- clips:

  Data frame of clip rows.

## Value

Standardised clip data frame ready for `export_clip_rows`.
