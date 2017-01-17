#' @title Function to extract counts for only mutatiosignatures without position information
#'
#' @description In case the user wants to disregard the position information for the mutational
#' signatures, this function allows the user to filter the data by removing the positional infromation and
#' aggregating the counts for each mutation signature across all the position informations.
#'
#' @param counts The matrix of the signature counts obtained by \code{aggregate_bin_counts} or
#'  \code{club_signature_counts} functions.
#'
#' @return Returns a filtered matrix with same number of rows (samples) but reduced set of columns that only contains
#' the counts for each mutational signature information without the location information.
#' @keywords filter-signatures
#' @export
#'



filter_signatures_wo_location <- function(counts){

  leftflank <- grep("left", colnames(counts))
  rightflank <- grep("right", colnames(counts))
  if(length(leftflank) > 0 | length(rightflank) > 0){
    counts <- counts[, - c(leftflank, rightflank)]
  }
  names <- colnames(counts);

  names_mod <- as.vector(sapply(names, function(x) strsplit(x, "_")[[1]][1]));

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

