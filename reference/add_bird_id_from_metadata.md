# Add Bird ID from SAP Metadata

Joins the `bird_id` column from a SAP metadata table into a clip data
frame using shared key columns (`filename`, `day_post_hatch`, `label`).
Rows that cannot be matched are assigned `"unknown_bird"`.

## Usage

``` r
add_bird_id_from_metadata(clips, metadata)
```

## Arguments

- clips:

  Data frame of clip rows.

- metadata:

  SAP metadata data frame containing a `bird_id` column.

## Value

`clips` with a populated `bird_id` column.
