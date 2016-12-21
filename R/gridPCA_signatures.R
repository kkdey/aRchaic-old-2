#' @title Function to perform grid plots of PCA for signature counts
#'
#' @description This fucntion performs PCA on normalized signature counts data
#' obtained from \code{aggregate_bin_counts} or \code{club_signature_counts}
#' and plots the scatter plots for PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3 in a grid plot.
#'
#' @param signature_counts The matrix of the signature counts
#' @param labs The factor labels for the samples used for coloring in the PC plot
#' @param normalize If TRUE, we normalize by the total number of mutations in that sample (analogous to library
#' size normalization in RNA-seq)
#' @param cols The palette used for labeling different data sources in the PCA plot.
#'
#' @return Returns grid plot of the PCA plots
#' @keywords PCA_signatures
#' @import grid
#' @import gridExtra
#' @import limma
#' @import ggplot2
#' @importFrom stats prcomp predict
#' @importFrom graphics abline legend  lines  par  plot
#' @export


gridPCA_signatures <- function(signature_counts,
                               labs,
                               normalize=TRUE,
                               cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                        "hotpink","burlywood","darkkhaki","yellow","darkgray","deepskyblue",
                                        "brown4","darkorchid","magenta", "azure1","azure4"))
{
  if(length(labs) != dim(signature_counts)[1]){
    stop("the length of the labels vector must equal the number of rows in the data")
  }
  if(normalize){
      voom_signature_counts <- t(limma::voom(t(signature_counts))$E);
      pr <- prcomp(voom_signature_counts)
  }

  else{
      pr <- prcomp(signature_counts)
  }

  pc_data_frame <- data.frame("PC"=pr$x,
                              "labels"= labs)

  graphList <- vector(mode="list");


  graphList[[1]] <- ggplot2::qplot(PC.PC1, PC.PC2,
                          data=pc_data_frame,
                          colour=factor(labels))  + scale_colour_manual(values = cols)

  graphList[[2]] <- ggplot2::qplot(PC.PC1, PC.PC3,
                         data=pc_data_frame,
                         colour=factor(labels)) + scale_colour_manual(values = cols)

  graphList[[3]] <- ggplot2::qplot(PC.PC2, PC.PC3,
                         data=pc_data_frame,
                         colour=factor(labels)) + scale_colour_manual(values = cols)

  a <- do.call("grid.arrange",
               args = list(grobs=graphList,
                           ncol = 2,
                           nrow = 2))

  return(pr)
}
