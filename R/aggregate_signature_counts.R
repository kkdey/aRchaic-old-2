#' @title Aggregate mutationFileFormat for multiple files into a matrix of counts of signatures.
#'
#' @description For each file of the directory, performs the \code{damage.build.counts}
#' function and aggregates the data into counts matrix format with the number of rows
#' equal to the number of files in the directory and the columns corresponding to the
#' number of different mutation signatures recorded. The files must be in .csv format.
#'
#' @param dir The directory containing the files that are to be read and aggregated over.
#' @param breaks The breaks used for the \code{damage.build.counts} function.
#' @param flanking_bases The number of flanking bases. Defaults to 2.
#' @param pattern An optional regular expression. Only file names which match the regular expression will be returned.
#' @param output_rda If non-NULL, the output matrix is stored in the path provided (should end with .rda).
#'
#' @return Returns a matrix of aggregate counts of each signature for each file arranged along a
#' row of the matrix. The number of rows correspond to the number of files in the directory
#' and the number of columns equalling the number of distinct mutation signatures.
#'
#' @keywords aggregate_counts
#'
#' @export


aggregate_signature_counts <- function(dir,
                                       breaks = NULL,
                                       flanking_bases = 1,
                                       pattern = NULL,
                                       output_rda = NULL){
  ancient_files <- list.files(dir, pattern = pattern)
  signature_ancient <- vector(mode="list")
  signature_counts_ancient <- vector(mode="list")


  for(num in 1:length(ancient_files)){
    tmp_dat <- damage_build_bin_counts(paste0(dir, ancient_files[num]),
                                       breaks=breaks,
                                       type=2)
    signature_counts_ancient[[num]] <- tmp_dat[,2];
    signature_ancient[[num]] <- as.character(tmp_dat[,1]);
    cat("Reading file ", ancient_files[num], "\n")
  }

  merged_signature_ancient <- signature_ancient[[1]]

  if(length(ancient_files) > 2){
    for(num in 2:length(ancient_files)){
      merged_signature_ancient <- union(merged_signature_ancient, signature_ancient[[num]])
    }
  }

  ancient_counts <- matrix(0, length(ancient_files), length(merged_signature_ancient))

  for(num in 1:length(ancient_files)){
    ancient_counts[num, match(signature_ancient[[num]], merged_signature_ancient)] <- signature_counts_ancient[[num]]
  }

  ancient_files_filt <- as.character(sapply(ancient_files, function(x) strsplit(x, ".csv")[[1]][1]))

  rownames(ancient_counts) <- ancient_files_filt

  signature_split <- do.call(rbind, lapply(merged_signature_ancient, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases+6)]))

  indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])

  indices2 <- numeric()
  for(m in 1:(4+2*flanking_bases)){
    indices2 <- c(indices2, which(signature_split[,m]=="N"));
  }

  indices3 <- numeric()
  indices3 <- c(indices3, which(signature_split[,(4+2*flanking_bases + 4)]=="N"))
  indices3 <- c(indices3, which(signature_split[,(4+2*flanking_bases + 6)]=="N"))


  indices <- union(indices1, union(indices2, indices3))

  if(length(indices) > 0) {
    ancient_counts_filtered <- ancient_counts[,-indices]
  }else{
    ancient_counts_filtered <- ancient_counts
  }

  rownames(ancient_counts_filtered) <- ancient_files_filt
  if(length(indices) > 0){
    colnames(ancient_counts_filtered) <- merged_signature_ancient[-indices]
  }else{
    colnames(ancient_counts_filtered) <- merged_signature_ancient
  }

  if(is.null(output_rda)){
    return(ancient_counts_filtered)
  }else{
    save(ancient_counts_filtered, output_rda)
  }
}

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
}

