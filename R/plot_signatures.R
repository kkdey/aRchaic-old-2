#' @title Function to check if the signatures should be clubbed
#'
#' @description Validity check whether the C->T and G->A signatures should be clubbed. This function
#' reads in the matrix of signature counts produced by the \code{aggregate_bin_counts} and plots the
#' C->T and G->A mutations for each sample.
#'
#' @param mat The matrix of counts of all signatures as produced by \code{aggregate_bin_counts}.
#' @return Returns a plot of number of C->T versus number of G->A mutations for each sample.
#'
#' @keywords plot_signatures
#'
#' @export
#'
#'

