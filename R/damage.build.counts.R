#' @title Summarize the mutationFileFormat data into bins along the read length
#'
#' @description Build a summary counts table by defining bins along the read length
#' based on the distance from the ends of the read and then sumaarizing the mutationFile
#' format data into bins by counting the number of occurrences of reads into each
#' bin.
#'
#' @param file a MutationFile format data obtained from applying
#' \code{generate_summary_data} on the BAM files.
#'
#' @param breaks The breaks along the read length used to form the bins. The default os
#' is NULL in which case, we choose the default option of having evenly separated
#' bins defined at spacings of 5 along the read length.
#'
#' @param type Two designs of output matrices available, according to user preference.
#'
#' @return Returns a matrix with columns representing the signature and the bins
#' in which the signature belongs to and a column representing the counts of the
#' occurrences of that signature in that bin  across the entire BAM file.
#'
#' @keywords summarize_counts
#' @importFrom utils read.csv
#' @importFrom plyr summarise ddply
#'
#' @export



damage_build_bin_counts =  function(file,
                                breaks=NULL,
                                type=1)
{
  file <-  read.csv(file=file, header=FALSE);
  file[which(file[,3]==-1), 3] <- 0
  file[which(file[,2]==-1), 2] <- 0

  min_dist_from_end <- apply(file[,2:3], 1, function(x) return(min(x)))


  if(is.null(breaks)){
    bins <- c(-1, 5, 10, 15, 20, max(min_dist_from_end))
  }else{
    message("breaks values provided : adjusting")
    bins <- c(intersect(breaks, (-1):max(min_dist_from_end)),max(min_dist_from_end))
  }

  bin_values <- .bincode(min_dist_from_end, bins, TRUE)

  modified_file <- cbind.data.frame(file[,1], file[,4], file[,5], file[,6], file[,7], bin_values)
  colnames(modified_file) <- c("pattern", "base.lsb", "base.rsb", "strand", "counts", "bin_values")

  library(plyr)
  df1 <- plyr::ddply(modified_file, .(pattern, strand, base.lsb, base.rsb, strand, bin_values), plyr::summarise, newvar = sum(counts))
  colnames(df1) <- c("pattern", "strand", "base.lsb", "base.rsb", "bin", "counts")

  if(type==2){
    df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,2], "_", df1[,3], "_", df1[,4], "_", df1[,5]), df1[,6])
    colnames(df2) <- c("pattern-strand-breaks-bin", "counts")
    out <- df2
  }else{
    out <- df1
  }

  return(out)
}
