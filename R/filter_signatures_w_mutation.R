#' @title Function to extract counts for only mutation with no flanking bases.
#'
#' @description In case the user wants to disregard the position information and flanking bases of the mutational
#' signatures, this function allows the user to filter the data by keeping only the mutation information and
#' aggregating the counts for each mutation across all the flanking bases and the positional information.
#'
#' @param counts The matrix of the signature counts obtained by \code{aggregate_bin_counts} or
#'  \code{club_signature_counts} functions.
#' @param flanking_bases The number of flanking bases for the mutation signatures
#'
#' @return Returns a filtered matrix with same number of rows (samples) but reduced set of columns that only contains
#' the counts for each mutation without the flanking bases and location information.
#' @keywords filter-signatures, mutation
#' @export
#'

filter_signatures_w_mutation <- function(counts, flanking_bases=2){

  leftflank <- grep("left", colnames(counts))
  rightflank <- grep("right", colnames(counts))
  if(length(leftflank) > 0 | length(rightflank) > 0){
    counts <- counts[, - c(leftflank, rightflank)]
  }

  names <- colnames(counts);

  names_mod <- as.vector(sapply(names, function(x) paste(strsplit(x, "")[[1]][(flanking_bases+1):(flanking_bases+4)], collapse="")))
  counts_w_mutation <- numeric();
  for(num in 1:dim(counts)[1]){
    counts_w_mutation <- rbind(counts_w_mutation,
                               tapply(counts[num,],
                                      factor(names_mod), sum))
  }

  rownames(counts_w_mutation) <- rownames(counts)

  if(sum((rowSums(counts) - rowSums(counts_w_mutation))^2) !=0){
    stop("the filtered matrix and the original matrix do not match")
  }

  return(counts_w_mutation)
}

