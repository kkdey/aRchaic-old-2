% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filter_signatures_wo_location.R
\name{filter_signatures_wo_location}
\alias{filter_signatures_wo_location}
\title{Function to extract counts for only mutatiosignatures without position information}
\usage{
filter_signatures_wo_location(counts)
}
\arguments{
\item{counts}{The matrix of the signature counts obtained by \code{aggregate_bin_counts} or
\code{club_signature_counts} functions.}
}
\value{
Returns a filtered matrix with same number of rows (samples) but reduced set of columns that only contains
the counts for each mutational signature information without the location information.
}
\description{
In case the user wants to disregard the position information for the mutational
signatures, this function allows the user to filter the data by removing the positional infromation and
aggregating the counts for each mutation signature across all the position informations.
}
\keyword{filter-signatures}

