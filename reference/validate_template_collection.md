# Validate Template Collection Object

An internal validation function to ensure the integrity and correct
structure of a template collection object.

## Usage

``` r
validate_template_collection(x)
```

## Arguments

- x:

  The template collection object to validate

## Value

Returns `TRUE` if the object passes all validation checks. Throws an
informative error if any validation fails.

## Details

Performs comprehensive checks on the template collection:

- Verifies the object is of class 'template_collection'

- Checks for the presence of required components

- Validates the structure of the template_info data frame

- Ensures template_list and template_matches are lists

Required columns in template_info include:

- template_name

- start_time

- end_time

- duration

- freq_min

- freq_max

- threshold

- clip_name

- clip_path

- source_file

- source_file_path

- creation_date

## Note

This is an internal validation function not intended for direct user
calls. It is used automatically during template collection creation and
manipulation.
