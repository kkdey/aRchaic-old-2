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
  ancient_files <- list.files(dir)
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

  indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])

  indices2 <- numeric()
  for(m in 1:(4+2*flanking_bases)){
    indices2 <- c(indices2, which(signature_split[,m]=="N"));
  }

  indices <- union(indices1, indices2)

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

  return(ancient_counts_filtered)
}
