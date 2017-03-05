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
#' @param strand Filters by strand if strand is provided to be "+" or "-", thereby fixing a strand.
#'               The default is "both", which means mutations from both strands are included.
#' @param numBases Number of bases from the ends of the read to be plotted.
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
                         plot_type=c("left", "right"),
                         sample_name = NULL,
                         strand = "both",
                         numBases = 20,
                         cols= c("black", "red", "green", "blue", "orange", "purple"),
                         legend_cex = 0.5,
                         main_cex = 0.7)
{
  data <- read.csv(file, header=FALSE)
  totsum <- sum(data[,7]);

  ######## putting G->A and C->T t the beginning ######################

  if(plot_type == "left"){
    index2 <- which(pattern == "G->A")
    if(length(index2) > 0){
      cols <- c(cols[index2], cols[-index2])
      pattern <- c(pattern[index2], pattern[-index2])
    }

    index1 <- which(pattern == "C->T")
    if(length(index1) > 0){
      cols <- c(cols[index1], cols[-index1])
      pattern <- c(pattern[index1], pattern[-index1])
    }
  }

  if(plot_type == "right"){
    index2 <- which(pattern == "C->T")
    if(length(index2) > 0){
      cols <- c(cols[index2], cols[-index2])
      pattern <- c(pattern[index2], pattern[-index2])
    }

    index1 <- which(pattern == "G->A")
    if(length(index1) > 0){
      cols <- c(cols[index1], cols[-index1])
      pattern <- c(pattern[index1], pattern[-index1])
    }
  }


###########   if strand information, focus on only one strand ###########


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

    if(min(pattern_data[,3]) == 1){ pattern_data[,3] = pattern_data[,3] - 1}

    pattern_data <- pattern_data[pattern_data[,2] < numBases | pattern_data[,3] < numBases, ]

    tab_pattern_left <- tapply(pattern_data[pattern_data[,2] < numBases,7], pattern_data[pattern_data[,2] < numBases,2], sum);

    lo_left <- loess(as.numeric(tab_pattern_left)~as.numeric(names(tab_pattern_left)))

    tab_pattern_right <- tapply(pattern_data[pattern_data[,3] < numBases,7], pattern_data[pattern_data[,3] < numBases,3], sum);

    lo_right <- loess(as.numeric(tab_pattern_right)~as.numeric(names(tab_pattern_right)))

    if(i==1){

    if(plot_type=="left"){
      if(is.null(sample_name)){
        temp <- predict(lo_left)
        temp[temp < 0] =0
        if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
        plot(temp, col=cols[i], ylim= c(0,max(temp)+0.0005),
             lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from left)",
             xaxs="i",yaxs="i",
             main=paste0("Damage plot (5' end) : "),
             cex.main = main_cex
        )
      }
      if(!is.null(sample_name)){
        temp <- predict(lo_left)
        temp[temp < 0] =0
        if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
        plot(temp, col=cols[i], ylim= c(0,max(temp)+0.0005),
             xlim = c(1,numBases), lwd=2, type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from left)",
             xaxs="i",yaxs="i",
             main=paste0("Damage plot (5' end): ", sample_name),
             cex.main = main_cex
        )}
    }

    if(plot_type=="right"){
      if(is.null(sample_name)){
        temp <- predict(lo_right)
        temp[temp < 0] =0
        if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
        plot(temp, col=cols[i], lwd=2, ylim= c(0,max(temp)+0.0005),
             type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from right)",
             main=paste0("Damage plot (3' end): "),
             cex.main = main_cex,
             xaxs="i",yaxs="i",
             #xlim=rev(range(as.numeric(names(tab_pattern_right))))
             xlim = c(numBases, 1))
      }
      if(!is.null(sample_name)){
        temp <- predict(lo_right)
        temp[temp < 0] =0
        if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
        plot(temp, col=cols[i], lwd=2, ylim= c(0,max(temp)+0.0005), type="l",
             ylab=paste0("predicted no. of damages: "),
             xlab="read positions (from right)",
             main=paste0("Damage plot (3' end): ", sample_name),
             cex.main = main_cex,
             xaxs="i",yaxs="i",
            # xlim=rev(range(as.numeric(names(tab_pattern_right))))
             xlim = c(numBases, 1))
      }
    }

    }else{

      if(plot_type=="left"){
        if(is.null(sample_name)){
          temp <- predict(lo_left)
          temp[temp < 0] =0
          if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
          lines(temp, col=cols[i], lwd=2)
        }
        if(!is.null(sample_name)){
          temp <- predict(lo_left)
          temp[temp < 0] =0
          if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
          lines(temp, col=cols[i], lwd=2)
        }
      }

      if(plot_type=="right"){
        if(is.null(sample_name)){
          temp <- predict(lo_right)
          temp[temp < 0] =0
          if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
          lines(temp, col=cols[i], lwd=2)
        }
        if(!is.null(sample_name)){
          temp <- predict(lo_right)
          temp[temp < 0] =0
          if(max(temp) < 1e-03){ temp <- temp*(1e-03/max(temp))}
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

