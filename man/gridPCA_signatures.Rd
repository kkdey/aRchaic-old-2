% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gridPCA_signatures.R
\name{gridPCA_signatures}
\alias{gridPCA_signatures}
\title{Function to perform grid plots of PCA for signature counts}
\usage{
gridPCA_signatures(signature_counts, labs, normalize = TRUE, cols = c("red",
  "blue", "darkgoldenrod1", "cyan", "firebrick", "green", "hotpink",
  "burlywood", "darkkhaki", "yellow", "darkgray", "deepskyblue", "brown4",
  "darkorchid", "magenta", "azure1", "azure4"))
}
\arguments{
\item{signature_counts}{The matrix of the signature counts}

\item{labs}{The factor labels for the samples used for coloring in the PC plot}

\item{normalize}{If TRUE, we normalize by the total number of mutations in that sample (analogous to library
size normalization in RNA-seq)}

\item{cols}{The palette used for labeling different data sources in the PCA plot.}
}
\value{
Returns grid plot of the PCA plots
}
\description{
This fucntion performs PCA on normalized signature counts data
obtained from \code{aggregate_bin_counts} or \code{club_signature_counts}
and plots the scatter plots for PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3 in a grid plot.
}
\keyword{PCA_signatures}

