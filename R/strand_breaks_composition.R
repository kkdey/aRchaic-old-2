#' @title Computes strand break composition for mutation signatures at the extremes
#'
#' @description This function computes the composition of nucleotides at the strand breaks positions corresponding
#' to mutation signatures involving strand breaks.
#'
#' @param file The file name for which the strand break composition is sought
#' @param flanking_bases The number of flanking bases in the signatures
#'
#' @return Returns a list with two elements, each being a vector of counts of nucleotides mapped to the left and right strand breaks
#' @keywords filter-signatures
#' @export
#'


strand_breaks_composition <- function(file, flanking_bases=2){
  data <- read.csv(file, header=FALSE)
  tab <- as.numeric()
  for(m in 1:flanking_bases){
    tmp_indices <- which(data[,2] <=(min(data[,2])+(m-1)))
    tmp_data <- data[tmp_indices,]
    sig <- sapply(as.character(tmp_data[,1]), function(x) strsplit(x,"")[[1]][(flanking_bases -m +1)])
    sig <- sig[!is.na(match(sig, c("A", "C", "G", "T")))]
    tab <- rbind(tab, table(factor(sig, levels=c("A", "C", "G", "T"))))
  }

  sum_tab_left <- colSums(tab)

  tab <- as.numeric()
  for(m in 1:flanking_bases){
    tmp_indices <- which(data[,3] <=(min(data[,3])+(m-1)))
    tmp_data <- data[tmp_indices,]
    sig <- sapply(as.character(tmp_data[,1]), function(x) strsplit(x,"")[[1]][(flanking_bases + 5 + m -1)])
    sig <- sig[!is.na(match(sig, c("A", "C", "G", "T")))]
    tab <- rbind(tab, table(factor(sig, levels=c("A", "C", "G", "T"))))
  }

  sum_tab_right <- colSums(tab)

  ll <- list("left-strand-break"=sum_tab_left, "right-strand-break"=sum_tab_right)
  return(ll)
}
