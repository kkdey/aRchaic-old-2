#' @title Aggregate bin summary counts into a matrix of counts format
#'
#' @description For each file of the directory, performs the \code{damage.build.counts}
#' function and aggregates the data into counts matrix format with the number of rows
#' equal to the number of files in the directory and the columns corresponding to the
#' number of different mutation signatures recorded. The files must be in .csv format.
#'
#' @param dir The directory containing the files that are to be read and aggregated over.
#' @param breaks The breaks used for the \code{damage.build.counts} function.
#' @param flanking_bases The number of flanking bases. Defaults to 2.
#'
#' @return Returns a matrix of aggregate counts of each signature for each file arranged along a
#' row of the matrix. The number of rows correspond to the number of files in the directory
#' and the number of columns equalling the number of distinct mutation signatures.
#'
#' @keywords aggregate_counts
#'
#' @export
#'


aggregate_bin_counts <- function(dir,
                                 breaks = NULL,
                                 flanking_bases = 2){

  leftflank_files <- list.files(dir, pattern="leftflank")
  rightflank_files <- list.files(dir, pattern="rightflank")
  ancient_files <- setdiff(list.files(dir), union(leftflank_files, rightflank_files))


  signature_ancient <- vector(mode="list")
  signature_counts_ancient <- vector(mode="list")


  for(num in 1:length(ancient_files)){
    tmp_dat <- damage_build_bin_counts(paste0(dir, ancient_files[num]),
                                       breaks=breaks,
                                       type=2)
    signature_counts_ancient[[num]] <- tmp_dat[,2];
    signature_ancient[[num]] <- as.character(tmp_dat[,1]);
    cat("Reading file ", num, "\n")
  }

  merged_signature_ancient <- signature_ancient[[1]]

  for(num in 2:length(ancient_files)){
    merged_signature_ancient <- union(merged_signature_ancient, signature_ancient[[num]])
  }

  ancient_counts <- matrix(0, length(ancient_files), length(merged_signature_ancient))

  for(num in 1:length(ancient_files)){
    ancient_counts[num, match(signature_ancient[[num]], merged_signature_ancient)] <- signature_counts_ancient[[num]]
  }

  ancient_files_filt <- as.character(sapply(ancient_files, function(x) strsplit(x, ".csv")[[1]][1]))

  rownames(ancient_counts) <- ancient_files_filt

  signature_split <- do.call(rbind, lapply(merged_signature_ancient, function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  indices1 <- which(signature_split[,3]==signature_split[,6])

  indices2 <- numeric()
  for(m in 1:(4+2*flanking_bases)){
    indices2 <- c(indices2, which(signature_split[,m]=="N"));
  }

  indices <- union(indices1, indices2)

  ancient_counts_filtered <- ancient_counts[,-indices]
  rownames(ancient_counts_filtered) <- ancient_files_filt
  colnames(ancient_counts_filtered) <- merged_signature_ancient[-indices]

  message("Reading the left flanking bases")
  left_flank <- numeric()
  for(tmp in leftflank_files){
    dat <- read.csv(paste0(dir, tmp), header = FALSE)
    dat1 <- dat[match(c("A", "G", "C", "T"), dat[,1]),]
    left_flank <- rbind(left_flank, t(dat1[,2]))
  }

  rownames(left_flank) <- ancient_files_filt
  colnames(left_flank) <- paste0(c("A", "G", "C", "T"), "_", "left")


  message("Reading the right flanking bases")
  right_flank <- numeric()
  for(tmp in rightflank_files){
    dat <- read.csv(paste0(dir, tmp), header = FALSE)
    dat1 <- dat[match(c("A", "G", "C", "T"), dat[,1]),]
    right_flank <- rbind(right_flank, t(dat1[,2]))
  }

  rownames(right_flank) <- ancient_files_filt
  colnames(right_flank) <- paste0(c("A", "G", "C", "T"), "_", "right")

  pooled_data <- cbind(ancient_counts_filtered, left_flank, right_flank)
  return(pooled_data)
}
