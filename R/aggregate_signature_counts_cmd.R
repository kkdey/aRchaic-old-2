
######  Define the args to list all command line arguments

defaults <- list(flanking_bases = 1, output_rda = "output.rda", breaks = NULL)

####  parse the command line arguments into separate variables

for (arg in commandArgs(TRUE))
  eval(parse(text=arg))

for (nm in names(defaults))
  assign(nm, mget(nm, ifnotfound=list(defaults[[nm]]))[[1]])


###  redefining the variables not defined through command line

if(is.null(dir)){stop("The name of the directory hosting the
                      mutationFileFormat files must be provided")}


####  list all the files in the input directory
ancient_files <- list.files(dir)


signature_ancient <- vector(mode="list")
signature_counts_ancient <- vector(mode="list")

### The function to compute the counts of mutation signatures per file

damage_build_bin_counts_cmd =  function(input_file,
                                        breaks=NULL,
                                        type=2)
{
  file <-  read.csv(file=input_file, header=FALSE);

  ### assigning 0 to the positions with -1 labeled base

  file[which(file[,3]==-1), 3] <- 0
  file[which(file[,2]==-1), 2] <- 0

  ### calculating minimum distance of the signature from the end of the read

  min_dist_from_end <- apply(file[,2:3], 1, function(x) return(min(x)))

  ### if break values not provided, a default option is used, else they are
  ### modified to match with the minimum distances obtained from the file

  if(is.null(breaks)){
    bins <- c(-1, 5, 10, 15, 20, max(min_dist_from_end))
  }else{
    message("breaks values provided : adjusting")
    bins <- c(intersect(breaks, (-1):max(min_dist_from_end)),max(min_dist_from_end))
  }

  library(plyr)

  ### assign rach signature to some bin based on minimum distance from the end
  ### of the read
  ###  The .bincode function below is drawn from the plyr package

  bin_values <- .bincode(min_dist_from_end, bins, TRUE)

  ### report the following features of a mutation signature-
  ##  - pattern (mutation + flanking bases)
  ##  - left strand break
  ##  - right strand break
  ##  - which strand the signature comes from
  ##  - the counts of number of times that mutation signature occurs in the file
  ##  - bin_values : which bin the mutation signature falls in


  modified_file <- cbind.data.frame(file[,1], file[,4], file[,5], file[,6],
                                    file[,7], bin_values)
  colnames(modified_file) <- c("pattern", "base.lsb", "base.rsb",
                               "strand", "counts", "bin_values")


  ### pooling in the mutation signatures based on bins and patterns

  df1 <- plyr::ddply(modified_file, .(pattern, strand, base.lsb, base.rsb,
                                      strand, bin_values), plyr::summarise,
                      newvar = sum(counts))
  colnames(df1) <- c("pattern", "strand", "base.lsb",
                     "base.rsb", "bin", "counts")

  ### if type 2, then different properties of mutation signature are merged

  if(type==2){
    df2 <- cbind.data.frame(paste0(df1[,1], "_", df1[,2], "_", df1[,3],
                                   "_", df1[,4], "_", df1[,5]), df1[,6])
    colnames(df2) <- c("pattern-strand-breaks-bin", "counts")
    out <- df2
  }else{
    out <- df1
  }

  return(out)
}


###  looping through all the files and computing the mutation signature counts


for(num in 1:length(ancient_files)){
    tmp_dat <- damage_build_bin_counts_cmd(paste0(dir, ancient_files[num]),
                                           breaks=breaks,
                                           type=2)
    signature_counts_ancient[[num]] <- tmp_dat[,2];
    signature_ancient[[num]] <- as.character(tmp_dat[,1]);
    cat("Reading file ", num, "\n")
  }

merged_signature_ancient <- signature_ancient[[1]]

if(length(ancient_files) > 2){
  for(num in 2:length(ancient_files)){
    merged_signature_ancient <- union(merged_signature_ancient,
                                      signature_ancient[[num]])
  }
}

ancient_counts <- matrix(0, length(ancient_files),
                         length(merged_signature_ancient))

for(num in 1:length(ancient_files)){
    ancient_counts[num, match(signature_ancient[[num]], merged_signature_ancient)] <- signature_counts_ancient[[num]]
}

ancient_files_filt <- as.character(sapply(ancient_files,
                                          function(x) strsplit(x, ".csv")[[1]][1]))
rownames(ancient_counts) <- ancient_files_filt
signature_split <- do.call(rbind, lapply(merged_signature_ancient,
                                         function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

signature_split <- do.call(rbind, lapply(sig,
                                         function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases+6)]))

###  the indices for signatures with same base mutations

indices1 <- which(signature_split[,(flanking_bases+1)]==signature_split[,(flanking_bases+4)])


### removing all signatures with "N" base (null base) reported among mutation
###  signatures.

indices2 <- numeric()
for(m in 1:(4+2*flanking_bases)){
    indices2 <- c(indices2, which(signature_split[,m]=="N"));
}

### removing all signatures with "N" in the left or right base adjacent to the
### read start or end

indices3 <- numeric()
indices3 <- c(indices3, which(signature_split[,(4+2*flanking_bases + 4)]=="N"))
indices3 <- c(indices3, which(signature_split[,(4+2*flanking_bases + 6)]=="N"))



### pooling the indices with same base mutation and mutation signatures
### comprising of base "N"


indices <- union(indices1, union(indices2, indices3))


###  removal of the indices defined above which correspond to bad mutation signatures

if(length(indices) > 0) {
    ancient_counts_filtered <- ancient_counts[,-indices]
}else{
    ancient_counts_filtered <- ancient_counts
}

###  rows labeled by names of the ancient files (each file corresponds to each sample)
rownames(ancient_counts_filtered) <- ancient_files_filt

###  the mutation signature labels are assigned to the columns of the data

if(length(indices) > 0){
    colnames(ancient_counts_filtered) <- merged_signature_ancient[-indices]
}else{
    colnames(ancient_counts_filtered) <- merged_signature_ancient
}


##############  saving the data as a RDA file  ###########################


if(is.null(output_rda)){
  return(ancient_counts_filtered)
}else{
  save(ancient_counts_filtered, file=output_rda)
}

