#' @title Function to extract signature counts of mutations occurring only upto a certain position along the genome
#'
#' @description This function can be used by the user to extract the counts data for mutations occurring
#' from the end of the read upto a certain maximum distance.
#'
#' @param mat The matrix of the signature counts obtained by \code{aggregate_bin_counts} or
#'  \code{club_signature_counts} functions.
#' @param max_pos The highest position along the read from the start for which the counts are extracted.
#' @param flanking_bases The number of flanking bases in the signatures
#'
#' @return Returns a filtered matrix of signature counts of mutations upto a certain distance from ends of the read.
#' @keywords filter-signatures
#' @export
#'


filter_signatures_by_location <-  function(mat, max_pos = 23, flanking_bases=2){
  input_pos <- 1:max_pos;
  leftflank <- grep("left", colnames(mat))
  rightflank <- grep("right", colnames(mat))
  if(length(leftflank) > 0 | length(rightflank) > 0){
    mat <- mat[, - c(leftflank, rightflank)]
  }

  pos <- as.numeric(sapply(as.character(colnames(mat)), function(l)
  {
    sym <- strsplit(as.character(l), "")[[1]]
    return(paste(sym[((4+2*flanking_bases)+4):length(sym)], collapse=""))
  }))

  reduced_dat <- mat[,which(!is.na(match(pos, input_pos)))]
  return(reduced_dat)
}

