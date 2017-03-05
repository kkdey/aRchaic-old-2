#' @title Club signatures from the two strands to remove strand bias
#'
#' @description Function for clubbing signatures from the two strands to remove
#' strnad bias. An example is G->A and C->T are clubbed into C->T to remove the bias
#' due to whether 5' or 3' ends of the strands has been sequenced.
#'
#' @param signature_counts The matrix of counts of all signatures as produced by \code{aggregate_bin_counts}.
#' @param flanking_bases The number of flanking bases. Defaults to 2.
#'
#' @return Returns a matrix of clubbed signatures. The default choice of the mutation
#' signatures are C->T, C->A, C->G, T->A, T->C and T->G. The other sets of mutations
#' like G->A are clubbed to C->T to remove strand bias.
#'
#' @keywords aggregate_counts
#'
#' @export
#'


club_signature_counts <- function(signature_counts, flanking_bases=2){

  # leftflank <- grep("left", colnames(signature_pos_str_counts))
  # rightflank <- grep("right", colnames(signature_pos_str_counts))
  # signature_counts_flank <- signature_pos_str_counts[,c(leftflank, rightflank)]
  # signature_counts <- signature_pos_str_counts[, -c(leftflank, rightflank)]

  signature_set <- colnames(signature_counts)

  signature_set_split <- do.call(rbind, lapply(signature_set, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  indices_G <-  which(signature_set_split[,(flanking_bases+1)]=="G");
  indices_A <-  which(signature_set_split[,(flanking_bases+1)]=="A");

  indices <- c(indices_G, indices_A);

  signature_set_2 <- signature_set
  signature_set_2[indices] <- signatureclub(signature_set[indices], flanking_bases=flanking_bases)

  # signature_set_3 <- signature_set_2
  #
  # temp <- sapply(signature_set_2[indices], function(x) {
  #                                               bases <- strsplit(as.character(x), split="")[[1]]
  #                                               left_bases <- bases[1:(flanking_bases)]
  #                                               right_bases <- bases[(flanking_bases+4+1):(4+2*flanking_bases)]
  #                                               other_bases <- bases[(4+2*flanking_bases+1):nchar(x)]
  #                                               if(other_bases[2] == "-") {other_bases[2] = "+"} else {other_bases[2] = "-"}
  #                                               newsig <- paste0(c(rev(right_bases), bases[(flanking_bases+1):(flanking_bases+4)], rev(left_bases), other_bases), collapse="")
  #                                               return(newsig)
  #                                          })
  # signature_set_3[indices] <-  temp
  #
  #
  signature_counts_pooled <- do.call(rbind, lapply(1:dim(signature_counts)[1], function(x) tapply(signature_counts[x,], signature_set_2, sum)))
  rownames(signature_counts_pooled) <- rownames(signature_counts)
  temp_split <- do.call(rbind, lapply(colnames(signature_counts_pooled), function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  if(length(which(temp_split[,(flanking_bases+1)]=="G")) !=0 || length(which(temp_split[,(flanking_bases+1)]=="A"))!=0){
    stop("G->A conversion did not fully work; aborting")
  }
  return(signature_counts_pooled)
}


