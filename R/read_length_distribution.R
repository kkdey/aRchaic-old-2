#' @title Read length distribution plots corresponding to mutationFileFormat files
#'
#' @description For each file (of mutationFileFormat type) in a directory, we plot
#' the read length distribution under three settings - all reads, reads containing
#' C to T mutations and reads containing C to T towards the ends of the fragments.
#'
#' @param dir The directory housing all the mutationFileFormat files.
#' @param pattern An optional regular expression. Only file names which match the regular expression will be returned.
#' @param end_break The position on the read from the end of the fragment such that if
#' C->T change occurs beyond that end break, it is considered in the third class of
#' reads (potentially damaged reads)
#' @param plot_layout The layout of the plot with multiple images (deafult is a 3x3 layout)
#' @param cols The color profiles of the three plots. Defaults to "red", "green" and "blue".
#' @param title The title of the plot. Defaults to name of the sample.
#' @param cex_legend The cex of the legend.
#' @param cex_main The cex of the title of the figure.
#'
#' @return The read length distribution plots for each of the mutationFileFormat files
#' in the folder
#'
#' @keywords read_length
#'
#' @export

read_length_distribution <- function(dir,
                                     pattern = NULL,
                                     end_break = 5,
                                     plot_layout = c(1,1),
                                     cols = c("red", "green", "blue"),
                                     title = NULL,
                                     cex_legend = 0.5,
                                     cex.main = 0.5){

  ###  read in all the files in the directory

  files <- list.files(dir, pattern = pattern)
  tab <- list()
  for(l in files){
    tab[[l]] <- read.csv(paste0(dir, l), header=FALSE)
    cat("We are at sample", l, "\n")
  }

  ancient_names <- as.character(sapply(files, function(x) return(strsplit(x, "[_]")[[1]][1])))

  par(mfrow=plot_layout)
  for(l in 1:length(tab)){
    tab1 <- tab[[l]]
    read_length <- tab1$V2 + tab1$V3
    indices <- grep("C->T", tab1$V1)
    indices <- c(indices, grep("G->A", tab1$V1))
    read_length_CtoT <- tab1[indices, ]$V2 + tab1[indices, ]$V3 ## read lengths for all C to T
    indices2 <- which(tab1$V2 < end_break | tab1$V3 < end_break)
    indices_matched <- intersect(indices, indices2)
    read_length_3 <- (tab1[indices_matched, ]$V2 + tab1[indices_matched, ]$V3)
    if(is.null(title)){
      title <- ancient_names[l]
    }
    plot(table(read_length)/sum(table(read_length)), type="o", col=cols[1],
         main=paste0(title), cex.main = cex.main, xlab="read pos", ylab="prop of occur")
    points(table(read_length_CtoT)/sum(table(read_length_CtoT)), type="o", col=cols[2])
    points(table(read_length_3)/ sum(table(read_length_3)), type="o", col=cols[3])
    legend("topright", fill=c("red", "green", "blue"),
           legend = c("all reads ", "reads with C->T", paste0("reads with C->T ends < ", end_break)),
           cex=cex_legend)
    cat("we are at sample", l, "\n")
  }
}
