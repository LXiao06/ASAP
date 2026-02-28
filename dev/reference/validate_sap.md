# Validate a Sound Analysis Pro (SAP) Object

Internal function to perform comprehensive validation of a SAP object,
ensuring it meets the required structure and contains all necessary
components.

## Usage

``` r
validate_sap(x)
```

## Arguments

- x:

  An object to be validated as a SAP object

## Value

Returns `TRUE` if the object passes all validation checks. Throws an
informative error if any validation fails.

## Details

This is an internal validation function used primarily during SAP object
creation and manipulation. It performs rigorous checks on the object's
structure:

- Verifies the object is of class 'Sap'

- Checks for the presence of all required components

- Validates the type and structure of each component

- Ensures specific objects have the correct class

Required components include:

- metadata: A data frame with file and recording information

- base_path: A character string of the base directory

- motifs: A segment object for motif-level analysis

- bouts: A segment object for bout-level analysis

- syllables: A segment object for syllable-level analysis

- segments: A segment object for segment-level analysis

- templates: A template_collection object

- features: A list containing feature categories

- misc: A list of miscellaneous information

- version: A character string representing the SAP object version

## Note

This is an internal function not intended to be called directly by
users. It is used automatically during SAP object creation and
manipulation.
