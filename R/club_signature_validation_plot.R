#' @title Function to check if the signatures should be clubbed
#'
#' @description Validity check whether the C->T and G->A signatures should be clubbed. This function
#' reads in the matrix of signature counts produced by the \code{aggregate_bin_counts} and plots the
#' C->T and G->A mutations for each sample.
#'
#' @param signature_counts The matrix of counts of all signatures as produced by \code{aggregate_bin_counts}.
#' @param flanking_bases The number of flanking bases. Defaults to 2.
#' @param cex The size of the dots used for plotting
#' @param col The color of the dot
#' @param pch The shape of the dot
#' @param xlab The X-axis label
#' @param ylab The Y-axis label
#' @param lty The line type of the dot
#' @param lwd The line width of the dot
#' @param mar The margin of the plot
#' @param log if TRUE, the lof of the counts are plotted, else plot on original scale. Default is FALSE.
#' @return Returns a plot of number of C->T versus number of G->A mutations for each sample.
#'
#' @keywords plot_signatures
#'
#' @export
#'
#'

club_signature_validation_plot <- function(signature_counts, flanking_bases = 2, cex=1, col="red",
                                  pch=20, xlab="C->T", ylab= "G->A",
                                  lty=1, lwd=1, mar=c(5,4,4,4), log=FALSE){
  signature_set <- colnames(signature_counts)
  new_signature_set <- signatureclub(signature_set)

  signature_set_split <- do.call(rbind, lapply(signature_set, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  indices1 <-  which(signature_set_split[,3]=="C" & signature_set_split[,6]=="T")
  indices2 <-  which(signature_set_split[,3]=="G" & signature_set_split[,6]=="A")

  C_to_T_counts <- rowSums(signature_counts[,indices1]);
  G_to_A_counts <- rowSums(signature_counts[,indices2]);

  if(!log){
     par(mar=mar)
     plot(C_to_T_counts, G_to_A_counts, xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col,
       lty=lty, lwd=lwd, main="Plot of C->T and G->A counts")
     abline(0,1)
  }else{
    par(mar=mar)
    plot(log(C_to_T_counts+1), log(G_to_A_counts+1), xlab=xlab, ylab=ylab, pch=pch, cex=cex, col=col,
         lty=lty, lwd=lwd, main = "Plot of C->T and G->A counts (log scale)")
    abline(0,1)
  }

  ll <- list("CtoT" = C_to_T_counts, "GtoA" = G_to_A_counts)
  return(ll)
}



gsub2 <- function(pattern, replacement, x, ...) {
  for(i in 1:length(pattern))
    x <- chartr(pattern[i], replacement[i], x, ...)
  x
}

signatureclub <- function(signature_set, flanking_bases){
  from <- c('A','T','G','C')
  to <- c('t','a','c','g');
  signature_set_mod <- array(0, length(signature_set));
  for(m in 1:length(signature_set)){
    temp <- toupper(gsub2(from, to, signature_set[m]))
    temp_split <- strsplit(as.character(temp), split="")[[1]]
    side1 <- temp_split[1:flanking_bases]
    side2 <- temp_split[(5+flanking_bases):(4+2*flanking_bases)]
    temp_split[(5+flanking_bases):(4+2*flanking_bases)] <- rev(side1)
    temp_split[1:flanking_bases] <- rev(side2)
    sign <- temp_split[6+2*flanking_bases]
    if(sign=="-"){sign = "+"}else if(sign=="+"){sign = "-"}
    leftb <- temp_split[8+2*flanking_bases]
    rightb <- temp_split[10+2*flanking_bases]
    temp_split[10+2*flanking_bases] <- leftb
    temp_split[8+2*flanking_bases] <- rightb
    temp_split[6+2*flanking_bases] <- sign
    temp_new <- paste0(temp_split, collapse = "")
    signature_set_mod[m] <- temp_new
  }
  return(signature_set_mod)
}
