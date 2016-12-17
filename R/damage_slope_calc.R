#' @title Calculates slope of the damage/mutational pattern from ends of reads
#'
#' @description Calculates the slope of the damage patterns from ends of the reads for any mutational signature.
#'
#' @param file a MutationFile format data obtained from applying
#' \code{generate_summary_data} on the BAM files.
#'
#' @param breaks The breaks along the read length used to form the bins. The default os
#' is NULL in which case, we choose the default option of having evenly separated
#' bins defined at spacings of 5 along the read length.
#'
#' @param type two options - slope or relative ratio. If equal to slope (which is also the default),
#' calculates the slope of the line connecting each position on the read with respect to the previous
#' position on the read. If it is relative ratio, calculates relative ratio of damages at each position
#' on the read with respect to the last position.
#'
#' @return Returns a matrix of slope or relative ratio values at each position from the start and end of the reads
#'
#' @keywords summrize_counts
#'
#' @export

damage_slope_calc <-  function(file, breaks = c(1, 3, 5, 7, 10, 15), type="slope"){
  levels <- c("C->T", "C->G", "C->A", "T->G", "T->A", "T->C",
              "G->A", "G->C", "G->T", "A_>G", "A->T", "A->C");
  data <- read.csv(file, header=FALSE)
  counts_gen <- numeric()
  flag <- character();

  for(num in 1:length(breaks)){
    if(num ==1){
      data_red <- data[which(data[,2] <= breaks[num]),];
      mut_prof <- substring(data_red[,1],3,6)
      out <- tapply(data_red[,4],  factor(mut_prof, levels = levels), sum)
      out[is.na(out)] = 0;
      flag <- c(flag, paste0(-1, "-", breaks[num]));
      counts_gen <- rbind(counts_gen, out)
    }else{
      data_red <- data[(data[,2] <= breaks[num] & data[,2] > breaks[(num-1)]),];
      mut_prof <- substring(data_red[,1],3,6)
      out <- tapply(data_red[,4], factor(mut_prof, levels = levels), sum)
      out[is.na(out)] = 0;
      flag <- c(flag, paste0(breaks[(num-1)], "-", breaks[num]));
      counts_gen <- rbind(counts_gen, as.numeric(out))
    }
  }
  rownames(counts_gen) <- flag;

  slope_counts_gen <- matrix(0, dim(counts_gen)[1], dim(counts_gen)[2])
  slope_counts_gen[1,] <- rep(1,dim(counts_gen)[2])
  for(m in 2:dim(counts_gen)[1]){
    if(type=="ratio"){
      slope_counts_gen[m, ] <- (counts_gen[m,]+1e-15)/(counts_gen[1,]+1e-15);
    }else if (type=="slope"){
      slope_counts_gen[m, ] <- (counts_gen[m,] - counts_gen[(m-1),])/(breaks[m] - breaks[(m-1)]);
    }
  }
  colnames(slope_counts_gen) <- colnames(counts_gen)
  rownames(slope_counts_gen) <- rownames(counts_gen)
  ll1 <- list("counts"=counts_gen, "slope-counts"=slope_counts_gen);

  counts_gen <- numeric()
  flag <- character();
  for(num in 1:length(breaks)){
      if(num ==1){
        data_red <- data[which(data[,3] <= breaks[num]),];
        mut_prof <- substring(data_red[,1],3,6)
        out <- tapply(data_red[,4],  factor(mut_prof, levels = levels), sum)
        out[is.na(out)] = 0;
        flag <- c(flag, paste0(-1, "-", breaks[num]));
        counts_gen <- rbind(counts_gen, out)
      }else{
        data_red <- data[(data[,3] <= breaks[num] & data[,3] > breaks[(num-1)]),];
        mut_prof <- substring(data_red[,1],3,6)
        out <- tapply(data_red[,4], factor(mut_prof, levels = levels), sum)
        out[is.na(out)] = 0;
        flag <- c(flag, paste0(breaks[(num-1)], "-", breaks[num]));
        counts_gen <- rbind(counts_gen, as.numeric(out))
      }
    }
    rownames(counts_gen) <- flag;

    slope_counts_gen <- matrix(0, dim(counts_gen)[1], dim(counts_gen)[2])
    slope_counts_gen[1,] <- rep(1,dim(counts_gen)[2])
    if(type=="ratio"){
      slope_counts_gen[m, ] <- (counts_gen[m,]+1e-15)/(counts_gen[1,]+1e-15);
    }else if(type=="slope"){
      slope_counts_gen[m, ] <- (counts_gen[m,] - counts_gen[(m-1),])/(breaks[m] - breaks[(m-1)]);
    }
    colnames(slope_counts_gen) <- colnames(counts_gen)
    rownames(slope_counts_gen) <- rownames(counts_gen)

    ll2 <- list("counts"=counts_gen, "slope-counts"=slope_counts_gen);
    ll <- list("left"= ll1, "right"= ll2)
    return(ll)
}
