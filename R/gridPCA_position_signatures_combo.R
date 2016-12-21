#' @title Function to perform grid plots of PCA for signature counts with position information from multiple data sources
#'
#' @description This fucntion performs PCA on signature counts or proportions data (both normalized or unnormalized)
#' from multiple data soirces, and stored in a list of counts table as obtained from
#' \code{aggregate_bin_counts} or \code{club_signature_counts}. The data is clubbed into positional bins and
#' scatter plots for PC1 vs PC2, PC2 vs PC3 and PC1 vs PC3 are displayed in a grid.
#'
#' @param signature_list A list of signature counts or signature proportion matrices from multiple sources of data
#' @param labs The factor labels for the samples used for coloring in the PC plot
#' @param source_names The names of the list of signature counts (equals the number of data sources)
#' @param input_pos The number of positions to be taken into account for building the PCA
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
#'


gridPCA_position_signatures_combo <- function(signature_list,
                               labs,
                               source_names=NULL,
                               input_pos = 1:20,
                               normalize=TRUE,
                               cols = c("red","blue","darkgoldenrod1","cyan","firebrick", "green",
                                        "hotpink","burlywood","yellow","darkgray","deepskyblue","darkkhaki",
                                        "brown4","darkorchid","magenta","yellow", "azure1","azure4"))
{
  if(!is.null(source_names)){
    if(length(source_names) != length(signature_list)){
      stop("the length of the source names must equal to the length of the list of signatures")
    }
    names(signature_list) <- source_names;
  }

  mat <- as.numeric();
  labs <- c();
  for(len in 1:length(ll)){

    pos <- as.numeric(sapply(as.character(colnames(ll[[len]])), function(l)
    {
      sym <- strsplit(as.character(l), "")[[1]]
      return(paste(sym[10:length(sym)], collapse=""))
    }))
   # cat(min(pos), "\n")
    pos <- pos - min(pos)
    reduced_dat <- ll[[len]][,which(!is.na(match(pos, input_pos)))]
    pos2 <- pos[which(!is.na(match(pos, input_pos)))];
    pos2fac <- factor(pos2, levels=input_pos)

    for(m in 1:dim(ll[[len]])[1]){
      mat <- rbind(mat, tapply(reduced_dat[m,], pos2fac, sum));
      labs <- c(labs, names(ll)[len])
    }
  }

  out <- gridPCA_signatures(mat, factor(labs), normalize = normalize, cols=cols)
  ll <- list("PC_project" = out,
             "position_matrix"=mat,
             "labs"=labs)
  return(ll)
}




