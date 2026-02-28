# Internal Constructor for Segment Object

An internal constructor function for creating a segment object with a
standardized structure for storing acoustic segment information.

## Usage

``` r
new_segment(x = data.frame())
```

## Arguments

- x:

  A data frame containing segment information (default is an empty data
  frame)

## Value

A segment object of class "segment" with the specified segment
information

## Details

Creates a segment object with a predefined structure. If no data frame
is provided, it initializes an empty segment with the following columns:

- filename: Name of the source audio file

- day_post_hatch: Numeric day post-hatch

- label: Categorical label for the segment

- start_time: Numeric start time of the segment

- end_time: Numeric end time of the segment

- duration: Numeric duration of the segment

When a data frame is provided, the function:

- Validates the presence of required columns

- Converts columns to appropriate data types

- Calculates segment duration if not provided

## Note

This is an internal function not intended to be called directly by
users. It is used within the package for creating segment objects.
