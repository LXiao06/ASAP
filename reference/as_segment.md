# Convert data to a segment object

Convert data to a segment object

## Usage

``` r
as_segment(x)
```

## Arguments

- x:

  A data frame to convert to a segment object

## Value

A validated segment object

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(
  filename = "song.wav",
  day_post_hatch = 1,
  label = "a",
  start_time = 0,
  end_time = 1,
  duration = 1
)
segment <- as_segment(df)
} # }
```
