% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectrogramming.R
\name{visualize_song}
\alias{visualize_song}
\alias{visualize_song.default}
\alias{visualize_song.Sap}
\title{Visualize Song Data}
\usage{
visualize_song(x, ...)

\method{visualize_song}{default}(
  x,
  start_time_in_second = NULL,
  end_time_in_second = NULL,
  fft_window_size = 1024,
  overlap = 0.5,
  dark_mode = TRUE,
  legend = FALSE,
  keep.par = FALSE,
  verbose = TRUE,
  ...
)

\method{visualize_song}{Sap}(
  x,
  template_clips = FALSE,
  indices = NULL,
  n_samples = NULL,
  random = TRUE,
  start_time_in_second = NULL,
  end_time_in_second = NULL,
  fft_window_size = 1024,
  overlap = 0.75,
  keep.par = TRUE,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{An object to visualize, either a file path or SAP object}

\item{...}{Additional arguments passed to specific methods}

\item{start_time_in_second}{Numeric start time in seconds}

\item{end_time_in_second}{Numeric end time in seconds}

\item{fft_window_size}{Size of FFT window (default: 512 for default, 1024 for SAP)}

\item{overlap}{Overlap between windows (default: 0.5 for default, 0.75 for SAP)}

\item{dark_mode}{For default method: Use dark theme (default: TRUE)}

\item{legend}{For default method: Show spectrogram legend (default: FALSE)}

\item{keep.par}{Preserve plotting parameters}

\item{verbose}{Print processing messages}

\item{template_clips}{Logical. For SAP objects: whether to visualize original songs (FALSE)
or template clips (TRUE) (default: FALSE)}

\item{indices}{For SAP objects: Numeric vector of specific indices to visualize}

\item{n_samples}{For SAP objects: Number of samples to visualize if indices is NULL.
Default is 6 or max available}

\item{random}{For SAP objects: Randomly sample songs if TRUE}
}
\value{
Generates spectrogram plot(s) and returns the input object invisibly.
}
\description{
Creates spectrograms and visualizes acoustic data from WAV files or SAP objects.
}
\details{
For WAV files:
\itemize{
\item Creates single spectrogram using FFmpeg's FFT
\item Customizable time range and FFT settings
\item Optional dark mode and legend
}

For SAP objects:
\itemize{
\item Creates multi-panel spectrograms
\item Supports random or sequential sampling
\item Maintains plotting state for sequential viewing
\item Adds day and label information to plots
}
}
\examples{
\dontrun{
# Visualize a single WAV file
visualize_song("path/to/song.wav",
               start_time_in_second = 10,
               end_time_in_second = 20)

# Basic visualization from SAP object
visualize_song(sap_object, n_sample = 4)

# Visualize specific indices with custom FFT settings
visualize_song(sap_object,
               indices = c(1, 3, 5),
               fft_window_size = 2048,
               overlap = 0.8)

# Sequential visualization with time ranges
visualize_song(sap_object,
               n_sample = 6,
               random = FALSE,
               start_time_in_second = rep(0, 6),
               end_time_in_second = rep(5, 6))
}

}
