#' @title Plot mutation signature patterns across read length
#'
#' @description This function reads in the MutationFile format data comprising
#' of the type of  mutation signatures and its location information from the
#' ends of reads, and plots a line plot of the frequency of a pattern provided
#' by the user, along the read length.
#'
#' @param file a MutationFile format data obtained from applying
#' \code{generate_summary_data} on the BAM files.
#'
#' @param pattern the mutation signature pattern, whose frequency along the
#' read, the user wants to plot.
#'
#' @param plot_type Three options are - "left", "right" and "both" demonstrating
#' whether the user wants the frequency plot of the pattern at the left end of the
#' read, right end of the read or both ends.
#'
#' @param use_prop If equal to 0, plots the pattern of counts along the reads. If equals to 1, plots
#' the proportion of the specific pattern changes against all other mutations along the reads. If equals to 2,
#' plots the proportion of the specific pattern against all mutations of the same pattern along the reads.
#'
#' @return Returns a frequency plot of the pattern across the read length.
#'
#' @keywords pattern, mutation signature
#'
#' @importFrom stats loess
#' @export
#'


pattern_plot <- function(file,
                         pattern,
                         plot_type=c("left", "right", "both"),
                         use_prop=0)
{
  data <- read.csv(file, header=FALSE)
  if(use_prop==1){
    data[,4] <- data[,4]/sum(data[,4])
  }
  pattern_data <- data[grep(pattern, data[,1]),]
  pattern_data[which(pattern_data[,3]==-1), 3] <- 0

  pattern_data <- data[grep(pattern, data[,1]),]

  pattern_data[which(pattern_data[,3]==-1), 3] <- 0

  if(use_prop==2){
    pattern_data[,4] <- pattern_data[,4]/sum(pattern_data[,4])
  }

  tab_pattern_left <- tapply(pattern_data[,4], pattern_data[,2], sum);

  lo_left <- loess(as.numeric(tab_pattern_left)~as.numeric(names(tab_pattern_left)))

  tab_pattern_right <- tapply(pattern_data[,4], pattern_data[,3], sum);

  lo_right <- loess(as.numeric(tab_pattern_right)~as.numeric(names(tab_pattern_right)))

  if(plot_type=="left"){
    par(mfrow=c(1,1))
    tmp <- predict(lo_left)
    plot(tmp, col='red', lwd=2, type="l",
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from left)",
         main=paste0("Damage plot across reads: ", pattern)
    )
    return(tmp)
  }

  if(plot_type=="right"){
    par(mfrow=c(1,1))
    tmp <- predict(lo_right)
    plot(predict(lo_right), col='red', lwd=2, type="l",
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from right)",
         main=paste0("Damage plot across reads: ", pattern),
         xlim=rev(range(as.numeric(names(tab_pattern_right)))))
    return(tmp)
  }

  if(plot_type=="both"){
    par(mfrow=c(1,2))
    tmp1 <- predict(lo_left);
    tmp2 <- predict(lo_right);
    plot(predict(lo_left), col='red', lwd=2, type="l",
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from left)",
         main=paste0("Damage plot across reads: ", pattern)
    )
    plot(predict(lo_right), col='red', lwd=2, type="l",
         ylab=paste0("predicted no. of damages: ", pattern),
         xlab="read positions (from right)",
         main=paste0("Damage plot across reads: ", pattern),
         xlim=rev(range(as.numeric(names(tab_pattern_right)))))
    tmp_list <- list("left"=tmp1, "right"=tmp2)
    return(tmp_list)
  }
}
