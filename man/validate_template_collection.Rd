% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/objects.R
\name{validate_template_collection}
\alias{validate_template_collection}
\title{Validate Template Collection Object}
\usage{
validate_template_collection(x)
}
\arguments{
\item{x}{The template collection object to validate}
}
\value{
Returns \code{TRUE} if the object passes all validation checks.
Throws an informative error if any validation fails.
}
\description{
An internal validation function to ensure the integrity and
correct structure of a template collection object.
}
\details{
Performs comprehensive checks on the template collection:
\itemize{
\item Verifies the object is of class 'template_collection'
\item Checks for the presence of required components
\item Validates the structure of the template_info data frame
\item Ensures template_list and template_matches are lists
}

Required columns in template_info include:
\itemize{
\item template_name
\item start_time
\item end_time
\item duration
\item freq_min
\item freq_max
\item threshold
\item clip_name
\item clip_path
\item source_file
\item source_file_path
\item creation_date
}
}
\note{
This is an internal validation function not intended for direct user calls.
It is used automatically during template collection creation and manipulation.
}
\keyword{internal}
