# Internal Constructor for Sound Analysis Pro (SAP) Object

An internal constructor function for creating a new SAP object with
predefined default components and structure.

## Usage

``` r
new_sap(
  metadata = data.frame(),
  base_path = character(),
  templates = new_template_collection(),
  motifs = new_segment(),
  bouts = new_segment(),
  syllables = new_segment(),
  segments = new_segment(),
  features = list(motif = list(), bout = list(), syllable = list(), segment = list()),
  misc = list(),
  version = "1.0.1"
)
```

## Arguments

- metadata:

  A data frame containing file and recording metadata (default is an
  empty data frame)

- base_path:

  A character string representing the base directory path (default is an
  empty character vector)

- templates:

  A template_collection object for storing analysis templates (default
  is a new empty template collection)

- motifs:

  A segment object for motif-level analysis (default is a new empty
  segment object)

- bouts:

  A segment object for bout-level analysis (default is a new empty
  segment object)

- syllables:

  A segment object for syllable-level analysis (default is a new empty
  segment object)

- segments:

  A segment object for general segment analysis (default is a new empty
  segment object)

- features:

  A list containing feature categories for different analysis levels
  (default is an empty list with predefined categories)

- misc:

  A list for storing miscellaneous information (default is an empty
  list)

- version:

  A character string representing the SAP object version (default is
  "1.0.1")

## Value

A SAP object of class "Sap" with the specified components

## Details

This is an internal constructor function used to create a standardized
SAP object. It ensures that all SAP objects have a consistent structure
with predefined components and default empty values.

The function creates a structure with the following components:

- metadata: Information about recorded files

- base_path: Directory path for the recordings

- motifs: Motif-level segment information

- bouts: Bout-level segment information

- syllables: Syllable-level segment information

- segments: General segment information

- templates: Collection of analysis templates

- features: Feature lists for different analysis levels

- misc: Miscellaneous metadata and information

- version: Package version

## Note

This is an internal function not intended to be called directly by
users. It is used within the package for creating SAP objects.
