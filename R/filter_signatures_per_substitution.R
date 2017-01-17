#' @title Function to extract signature counts only for a specific substitution pattern
#'
#' @description This function can be used by the user to extract the counts data only for mutation signatures with
#' a specific substitution pattern (for example C -> T).
#'
#' @param counts The matrix of the signature counts obtained by \code{aggregate_bin_counts} or
#'  \code{club_signature_counts} functions.
#' @param pattern The substitution pattern for which the user intends to extract the counts.
#' @param flanking_bases The number of flanking bases for the mutation signatures.
#' @param use_prop If TRUE, we record the proportion of the signatures with that pattern
#' compared to signatures with all patterns for each sample. Defaults to FALSE.
#'
#' @return Returns a filtered matrix of signature counts for a specific substitution pattern.
#' @keywords filter-signatures
#' @export
#'

filter_signatures_per_substitution <- function(counts, pattern, flanking_bases=2, use_prop=FALSE){

   leftflank <- grep("left", colnames(counts))
   rightflank <- grep("right", colnames(counts))
   if(length(leftflank) > 0 | length(rightflank) > 0){
      counts <- counts[, - c(leftflank, rightflank)]
   }

   mutation_sigs <- colnames(counts)
   sub_pattern <- sapply(mutation_sigs, function(x) paste(strsplit(as.character(x), "")[[1]][(flanking_bases+1):(flanking_bases+4)], collapse=""))
   indices <- which(!is.na(match(sub_pattern, pattern)))
   if(use_prop){
     prop <- t(apply(counts, 1, function(x) return(x/sum(x))))
     filtered_prop <- prop[,indices];
     return(filtered_prop)
   }else{
     filtered_counts <- counts[,indices];
     return(filtered_counts)
   }
}


