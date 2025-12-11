# Refine Motif Boundaries Using Segment Alignment

Adjusts motif boundaries based on underlying segments and optional
label-specific time adjustments. Integrates segment information to
determine precise motif onsets/offsets.

## Usage

``` r
refine_motif_boundaries(x, adjustments_by_label = NULL, verbose = TRUE)
```

## Arguments

- x:

  A SAP object containing 'segments' and 'motifs' data

- adjustments_by_label:

  Optional named list of time adjustments (in seconds) to apply to motif
  end limits, where names correspond to labels

- verbose:

  Logical flag for printing progress messages (default: TRUE)

## Value

Returns a modified SAP object with updated motifs containing:

- motif_onset - Precise start time based on segments

- motif_offset - Precise end time based on segments

- first_seg_index - First segment index in motif

- last_seg_index - Last segment index in motif

- motif_duration - Calculated motif duration

## Details

Key operations:

1.  Applies label-specific time adjustments to motif end limits

2.  Identifies segments contained within adjusted motif boundaries

3.  Calculates precise motif timing based on contained segments

4.  Preserves original motif structure while adding new timing columns

## Examples

``` r
if (FALSE) { # \dontrun{
# Apply 0.1s extension to "BL" motifs
sap <- refine_motif_boundaries(
  sap,
  adjustments_by_label = list(BL = 0.1)
)
} # }
```
