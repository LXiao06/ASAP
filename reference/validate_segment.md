# Validate Segment Object

An internal validation function to ensure the integrity and correctness
of a segment object.

## Usage

``` r
validate_segment(x)
```

## Arguments

- x:

  The segment object to validate

## Value

Returns `TRUE` if the object passes all validation checks. Throws an
informative error if any validation fails.

## Details

Performs comprehensive checks on the segment object:

- Verifies the object is of class 'segment'

- Checks for invalid time representations

- Ensures segment times are non-negative

- Validates that end times are not before start times

Specific validation checks include:

- Confirming the object is a segment class

- Checking that end times are greater than or equal to start times

- Ensuring all time values are non-negative

## Note

This is an internal validation function not intended for direct user
calls. It is used automatically during segment object creation and
manipulation.
