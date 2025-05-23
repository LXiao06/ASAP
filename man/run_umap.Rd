% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/umap.R
\name{run_umap}
\alias{run_umap}
\alias{run_umap.default}
\alias{run_umap.Sap}
\title{Run UMAP Dimensionality Reduction}
\usage{
run_umap(x, ...)

\method{run_umap}{default}(
  x,
  metadata_cols = NULL,
  scale = TRUE,
  n_neighbors = 15,
  n_components = 2,
  min_dist = 0.1,
  seed = 222,
  n_threads = NULL,
  verbose = TRUE,
  ...
)

\method{run_umap}{Sap}(
  x,
  segment_type = c("motifs", "syllables", "bouts", "segments"),
  data_type = c("spectral_feature", "spectrogram", "traj_mat"),
  label = NULL,
  scale = TRUE,
  n_neighbors = 20,
  n_components = 2,
  min_dist = 0.1,
  seed = 222,
  n_threads = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{An object to analyze, either a data frame or SAP object}

\item{...}{Additional arguments passed to specific methods}

\item{metadata_cols}{Column indices for metadata (for default method)}

\item{scale}{Whether to scale features before UMAP}

\item{n_neighbors}{Number of neighbors (default: 15)}

\item{n_components}{Number of output dimensions (default: 2)}

\item{min_dist}{Minimum distance parameter (default: 0.1)}

\item{seed}{Random seed for reproducibility}

\item{n_threads}{Number of computation threads}

\item{verbose}{Whether to print progress messages}

\item{segment_type}{For SAP objects: Type of segments to analyze ('motifs', 'syllables', 'bouts', 'segments')}

\item{data_type}{For SAP objects: Type of feature data ('spectral_feature','spectrogram', 'traj_mat')}

\item{label}{For SAP objects: Specific label to filter data}
}
\value{
For default method: Matrix of UMAP coordinates
For SAP objects: Updated SAP object with UMAP coordinates stored in features slot
}
\description{
A generic function to perform UMAP dimensionality reduction on feature data.
}
\details{
This generic function supports UMAP analysis through two methods:
\itemize{
\item Default method for feature data frames
\item SAP object method for organized song features
}
}
\examples{
\dontrun{
# Run UMAP on feature data frame
coords <- run_umap(features,
                   metadata_cols = c(1:5),
                   n_neighbors = 15)

# Run UMAP on SAP object
sap_obj <- run_umap(sap_object,
                    segment_type = "motifs",
                    data_type = "spectral_feature")

# UMAP with specific parameters
coords <- run_umap(features,
                   metadata_cols = 1:3,
                   scale = TRUE,
                   n_neighbors = 20,
                   seed = 123)

# UMAP with label filtering
sap_obj <- run_umap(sap_obj,
                    segment_type = "syllables",
                    data_type = "spectral_feature",
                    label = "a")
}

}
