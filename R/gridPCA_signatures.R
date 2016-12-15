#' @title Function to perform grid plots of PCA for signature counts
#'
#' @description This fucntion performs PCA on normalized signature counts data
#' obtained from \code{aggregate_bin_counts} or \code{club_signature_counts}
#' and plots the scatter plots for PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3 in a grid plot.
#'
#' @param signature_counts The matrix of the signature counts
#' @param labs The factor labels for the samples used for coloring in the PC plot
#'
#' @return Returns grid plot of the PCA plots
#' @keywords PCA_signatures
#' @import grid
#' @import gridExtra
#' @import limma
#' @export


gridPCA_signatures <- function(counts,
                               labs)
{
  if(length(labs) != dim(counts)[1]){
    stop("the length of the labels vector must equal the number of rows in the data")
  }
  voom_signature_counts <- t(limma::voom(t(counts))$E);
  pr <- prcomp(voom_signature_counts)


  pc_data_frame <- data.frame("PC"=pr$x,
                              "labels"= labs)

  graphList <- vector(mode="list");
  library(ggplot2)


  graphList[[1]] <- qplot(PC.PC1, PC.PC2,
                          data=pc_data_frame,
                          colour=labels)

  graphList[[2]] <-qplot(PC.PC1, PC.PC3,
                         data=pc_data_frame,
                         colour=labels)

  graphList[[3]] <-qplot(PC.PC2, PC.PC3,
                         data=pc_data_frame,
                         colour=labels)

  library(grid)
  library(gridExtra)
  a <- do.call("grid.arrange",
               args = list(grobs=graphList,
                           ncol = 2,
                           nrow = 2))

  return(pr)
}
