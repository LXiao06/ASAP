# Manual Syllable Labeling for Bird Song Analysis

Provides interactive or map-based assignment of alphabetic labels to
automatically clustered syllables. Supports both interactive console
input and predefined label maps, with persistent storage of label
assignments.

## Usage

``` r
manual_label(
  x,
  data_type = c("segment", "syllable"),
  label_map = NULL,
  interactive = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  A Sap object containing bird song analysis data

- data_type:

  Type of data to label: "segment" or "syllable"

- label_map:

  Optional data frame for non-interactive labeling with columns:

  - cluster: numeric cluster IDs

  - syllable: corresponding letter labels (a-z)

- interactive:

  Logical; force interactive mode even if label_map exists

- verbose:

  Logical; print progress messages and summaries

## Value

Returns modified Sap object with:

- Updated syllable labels in x\$syllables

- Stored label map in feature embeddings attributes

- All labels converted to lowercase

## Details

The function provides two main methods for syllable labeling:

1.  Interactive Console Mode:

    - Prompts for letter input for each cluster

    - Validates input (single letters a-z only)

    - Allows quitting mid-process ('q')

    - Automatically stores label map for future use

2.  Label Map Mode:

    - Uses predefined mapping of clusters to letters

    - Validates completeness and format of mapping

    - Updates stored label map

Label Storage and Retrieval:

- Label maps are stored in the Sap object attributes

- Automatically retrieved in subsequent runs

- Can be overridden with new map or interactive session

- Accessible via attr(x\$features\$\<data_type\>\$feat.embeds,
  "label_map")

## Label Map Format

The label_map data frame must have:

    data.frame(
      cluster = numeric_vector,  # Cluster IDs from auto_label
      syllable = character_vector  # Single letters a-z
    )

## Tips for Interactive Labeling

- View clusters first using plot_cluster(..., label_type = "pre")

- Consider similar acoustic features when assigning same letter

- Use 'q' to exit without saving if mistakes are made

- Review label summary after completion

- Store label_map externally for reproducibility

## Common Workflows

1.  First-time labeling:

        sap <- auto_label(sap)  # Generate clusters
        sap <- manual_label(sap, data_type = "syllable")  # Label interactively

2.  Using stored labels:

        sap <- manual_label(sap, data_type = "syllable")  # Uses stored map

3.  Updating labels:

        # Extract and modify existing map
        map <- attr(sap$features$syllable$feat.embeds, "label_map")
        map$syllable[map$cluster == 3] <- "d"  # Change cluster 3 to 'd'
        sap <- manual_label(sap, data_type = "syllable", label_map = map)

## See also

- auto_label() for automatic cluster generation

- plot_cluster() for visualizing clusters and syllables

## Examples

``` r
if (FALSE) { # \dontrun{
# Interactive labeling (first time)
sap <- manual_label(sap, data_type = "syllable")

# Using stored labels (subsequent runs)
sap <- manual_label(sap, data_type = "syllable")

# Force new interactive session
sap <- manual_label(sap,
                   data_type = "syllable",
                   interactive = TRUE)

# Using predefined label map
label_map <- data.frame(
  cluster = 1:5,
  syllable = c("a", "b", "b", "c", "a")
)
sap <- manual_label(sap,
                   data_type = "syllable",
                   label_map = label_map)

# Access stored label map
stored_map <- attr(sap$features$syllable$feat.embeds, "label_map")
print(stored_map)

# Typical workflow:
# 1. Run automatic clustering
sap <- auto_label(sap)

# 2. View clusters
plot_cluster(sap,
            data_type = "syllable",
            label_type = "pre")

# 3. Manually label clusters
sap <- manual_label(sap, data_type = "syllable")

# 4. View syllable labels
plot_cluster(sap,
            data_type = "syllable",
            label_type = "post")
} # }
```
