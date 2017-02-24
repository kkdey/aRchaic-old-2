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
#' @param sample_name The name of the BAM file to be used in title. Defaults to NULL.
#' @param cols The color vector of the mutational patterns plotted. The number of colors must be more than the number of
#' mutations used.
#' @param legend_cex The size of the legend which highlight the names of the patterns
#'
#' @return Returns a frequency plot of the pattern across the read length.
#'
#' @keywords pattern, mutation signature
#'
#' @importFrom stats loess
#' @importFrom graphics abline legend  lines  par  plot
#' @importFrom utils read.csv
#' @export
#'


pattern_plot_full <- function(file,
                         pattern,
                         plot_type=c("left", "right", "both"),
                         sample_name = NULL,
                         strand = "both",
                         cols= c("black", "red", "green", "blue", "orange", "purple"),
                         legend_cex = 0.5)
{
  data <- read.csv(file, header=FALSE)
  totsum <- sum(data[,3]);

  if(strand == "+"){
    data <- data[data[,6] == "+",]
  }
  if(strand =="_"){
    data <- data[data[,6] == "-",]
  }

  flag <- as.numeric();

  for(i in 1:length(pattern)){
    if(length(grep(pattern[i], data[,1])) > 0){
    flag <- c(flag, i);
    pattern_data <- data[grep(pattern[i], data[,1]),]
    pattern_data[,7] <- pattern_data[,7]/totsum;
    if(length(which(pattern_data[,3]==-1)) !=0){
    pattern_data[which(pattern_data[,3]==-1), 3] <- 0
    }

    tab_pattern_left <- tapply(pattern_data[,7], pattern_data[,2], sum);

    lo_left <- loess(as.numeric(tab_pattern_left)~as.numeric(names(tab_pattern_left)))

    tab_pattern_right <- tapply(pattern_data[,7], pattern_data[,3], sum);

    lo_right <- loess(as.numeric(tab_pattern_right)~as.numeric(names(tab_pattern_right)))

    if(i==1){

    if(plot_type=="left"){
      if(is.null(sample_name)){
        temp <- predict(lo_left)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from left)",
             main=paste0("Damage plot across reads: ")
        )
      }
      if(!is.null(sample_name)){
        temp <- predict(lo_left)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from left)",
             main=paste0("Damage plot across reads: ", sample_name)
        )}
    }

    if(plot_type=="right"){
      if(is.null(sample_name)){
        temp <- predict(lo_right)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from right)",
             main=paste0("Damage plot across reads: "),
             xlim=rev(range(as.numeric(names(tab_pattern_right)))))
      }
      if(!is.null(sample_name)){
        temp <- predict(lo_right)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from right)",
             main=paste0("Damage plot across reads: ", sample_name),
             xlim=rev(range(as.numeric(names(tab_pattern_right)))))
      }
    }


    if(plot_type=="both"){
      if(is.null(sample_name)){
        par(mfrow=c(1,2))
        temp <- predict(lo_left)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from left)",
             main=paste0("Damage plot across reads: ")
        )
        temp <- predict(lo_right)
        temp[temp < 0] =0
        plot(predict(lo_right), col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from right)",
             main=paste0("Damage plot across reads: "),
             xlim=rev(range(as.numeric(names(tab_pattern_right)))))
      }
      if(!is.null(sample_name)){
        par(mfrow=c(1,2))
        temp <- predict(lo_left)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from left)",
             main=paste0("Damage plot across reads: ", sample_name)
        )
        temp <- predict(lo_right)
        temp[temp < 0] =0
        plot(temp, col=cols[i], lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from right)",
             main=paste0("Damage plot across reads: ", sample_name),
             xlim=rev(range(as.numeric(names(tab_pattern_right)))))
      }
    }

    }else{

      if(plot_type=="left"){
        if(is.null(sample_name)){
          temp <- predict(lo_left)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
        }
        if(!is.null(sample_name)){
          temp <- predict(lo_left)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
        }
      }

      if(plot_type=="right"){
        if(is.null(sample_name)){
          temp <- predict(lo_right)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
        }
        if(!is.null(sample_name)){
          temp <- predict(lo_right)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
        }
      }



      if(plot_type=="both"){
        if(is.null(sample_name)){
          par(mfrow=c(1,2))
          temp <- predict(lo_left)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
          temp <- predict(lo_right)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
        }
        if(!is.null(sample_name)){
          temp <- predict(lo_left)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
          temp <- predict(lo_right)
          temp[temp < 0] =0
          lines(temp, col=cols[i], lwd=2)
        }
      }


     }
    }
  }
  if(plot_type == "left"){
    legend("topright", legend=pattern[flag], fill=cols[flag], cex=legend_cex);
  }else if (plot_type == "right"){
    legend("topleft", legend=pattern[flag], fill=cols[flag], cex=legend_cex);
  }
}

