#' @title Function to check if the signatures should be clubbed
#'
#' @description Validity check whether the C->T and G->A signatures should be clubbed. This function
#' reads in the matrix of signature counts produced by the \code{aggregate_bin_counts} and plots the
#' C->T and G->A mutations for each sample.
#'
#' @param signature_counts The matrix of counts of all signatures as produced by \code{aggregate_bin_counts}.
#' @return Returns a plot of number of C->T versus number of G->A mutations for each sample.
#'
#' @keywords plot_signatures
#'
#' @export
#'
#'

club_signature_validation_plot <- function(mat, cex=1, col="red",
                                  pch=20, xlab="C->T", ylab= "G->A",
                                  lty=1, lwd=1){
  signature_set <- colnames(signature_counts)
  new_signature_set <- signatureclub(signature_set)

  signature_set_split <- do.call(rbind, lapply(signature_set, function(x) strsplit(as.character(x), split="")[[1]]))

  indices1 <-  which(signature_set_split[,3]=="C" & signature_set_split[,6]=="T")
  indices2 <-  which(signature_set_split[,3]=="G" & signature_set_split[,6]=="A")

  C_to_T_counts <- rowSums(signature_counts[,indices1]);
  G_to_A_counts <- rowSums(signature_counts[,indices2]);

  par(mar=c(5,4,4,4))
  plot(C_to_T_counts, G_to_A_counts, xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col,
       lty=lty, lwd=lwd)
  abline(0,1)
}



gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- chartr(pattern[i], replacement[i], x, ...)
  x
}

signatureclub <- function(signature_set){
  from <- c('A','T','G','C')
  to <- c('t','a','c','g');
  signature_set_mod <- array(0, length(signature_set));
  for(m in 1:length(signature_set)){
    signature_set_mod[m] <- toupper(gsub2(from, to, signature_set[m]))
  }
  return(signature_set_mod)
}
