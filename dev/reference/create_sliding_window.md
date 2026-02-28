# Create Sliding Windows for Segment

Internal function to generate overlapping time windows for a segment.

## Usage

``` r
create_sliding_window(i, x, window_size = 0.1, step_size = 0.005, ...)
```

## Arguments

- i:

  Rendition number

- x:

  Data frame with segment information

- window_size:

  Window size in seconds

- step_size:

  Step size in seconds

- ...:

  Additional arguments

## Value

Data frame containing sliding window information
