#' @title Builds damage logo plots along with mutational profile across read for the clusters
#'
#' @description Damage Logo plots for each cluster representing the substitution frequency of the 6 substitutional
#' patterns (adjusting for strand bias) and the flanking bases arranged as per relative frequency on either side of
#' the substitution along with the probability of the cluster mutational profile along the read.
#'
#' @param theta The theta matrix obtained from running the grade of membership model that stores for each cluster, the
#' probability distribution over all the mutational signatures.
#' @param sig_names The mutational signature names. Defaults to the rownames of the theta matrix above.
#' @param ic.scale A binary variable indicating whether the height of the bars for substitution and flanking bases should be
#'        adjusted by the information criterion.
#' @param max_pos The maximum distance from the end of the read upto which mutations are considered.
#' @param flanking_bases The number of flanking bases of the mutational signature.
#' @param yscale_change A binary variable indicating whether the Y axis scale should be adjusted based on the size of the
#'        logos, defaults to TRUE.
#' @param xaxis A binary indicating whether the X axis of the logo plot should be shown
#' @param yaxis A binary indicating whether the Y axis of the logo plot should be shown
#' @param xaxis_fontsize The fontsize of the X axis ticks.
#' @param xlab_fontsize The fontsize of the X axis labels.
#' @param y_fontsize The fontsize of the Y axis ticks.
#' @param mut_width Thw width of the bar for the mutation at the center.
#' @param start The starting point of the stacking of logos on the Y axis. Should be close to 0, defau;ts to 0.0001.
#' @param pop_names The title of the plot. Defaults to the cluster labels.
#' @param logoport_x the X-axis position of the plot window for the logo plot
#' @param logoport_y the Y-axis position of the plot window for the logo plot
#' @param logoport_width the width of the plot window for the logo plot
#' @param logoport_height the width of the plot window for the logo plot
#' @param lineport_x the X-axis position of the plot window for the mutational profile line plot.
#' @param lineport_y the Y-axis position of the plot window for the mutational profile line plot.
#' @param lineport_width the width of the plot window for the mutational profile line plot.
#' @param lineport_height the width of the plot window for the mutational profile line plot.
#' @return Returns logo plots for each cluster
#'
#' @export

damageLogo_pos <- function(theta,
                       sig_names = NULL,
                       ic.scale=TRUE,
                       max_pos = 15,
                       flanking_bases=2,
                       yscale_change = FALSE,
                       xaxis=TRUE,
                       yaxis=TRUE,
                       xaxis_fontsize=10,
                       xlab_fontsize=15,
                       y_fontsize=15,
                       mut_width=2,
                       start=0.0001,
                       pop_names=paste0("Cluster ",1:dim(theta)[2]),
                       logoport_x = 0.4,
                       logoport_y= 0.5,
                       logoport_width= 0.4,
                       logoport_height= 0.5,
                       lineport_x = 0.95,
                       lineport_y=0.85,
                       lineport_width=0.25,
                       lineport_height=0.3){
  if(is.null(sig_names))
    sig_names <- rownames(theta)

  prob_mutation <- filter_signatures_only_location(t(theta), max_pos = max_pos, flanking_bases = flanking_bases)
  prob_mutation <- t(apply(prob_mutation, 1, function(x) return(x/sum(x))))

  sig_split <- do.call(rbind,
                       lapply(sig_names,
                              function(x) strsplit(as.character(x), split="")[[1]][1:(4+2*flanking_bases)]))

  ncol_sig <- (4+2*flanking_bases)

  if(flanking_bases%%1 != 0){
    stop("flanking bases not evenly distributed")
  }


  sub_pattern <- sapply(1:dim(sig_split)[1],
                        function(x) paste(sig_split[x,(flanking_bases+1):(flanking_bases+4)], collapse=""))

  new_sig_split <- cbind(sig_split[,1:flanking_bases], sub_pattern, sig_split[,((ncol_sig - flanking_bases +1):ncol_sig)])
  colnames(new_sig_split) = NULL

  prop_patterns_list <- list()

  for(l in 1:dim(theta)[2]){
    prop_patterns_list[[l]] <- numeric();
    for(j in 1:ncol(new_sig_split)){
      temp <- tapply(theta[,l], factor(new_sig_split[,j], levels=c("A", "C", "G", "T", "X",
                                                                   "C->T", "C->A", "C->G",
                                                                   "T->A", "T->C", "T->G")), sum)

      temp[is.na(temp)]=0
      prop_patterns_list[[l]] <- cbind(prop_patterns_list[[l]], temp)
    }
  }

  ic <- damage.ic(prop_patterns_list)

  grob_list <- list()
  for(l in 1:length(prop_patterns_list)){
  damageLogo.pos.skeleton(pwm = prop_patterns_list[[l]],
                           probs = prob_mutation[l,],
                           ic = ic[,l],
                           max_pos = max_pos,
                           ic.scale = ic.scale,
                           yscale_change = yscale_change,
                           xaxis=xaxis,
                           yaxis=yaxis,
                           xaxis_fontsize=xaxis_fontsize,
                           xlab_fontsize=xlab_fontsize,
                           y_fontsize=y_fontsize,
                           mut_width=mut_width,
                           start=start,
                           pop_name = pop_names[l],
                           logoport_x = logoport_x,
                           logoport_y= logoport_y,
                           logoport_width= logoport_width,
                           logoport_height= logoport_height,
                           lineport_x = lineport_x,
                           lineport_y= lineport_y,
                           lineport_width=lineport_width,
                           lineport_height=lineport_height)

  }
}

pwm2ic<-function(pwm) {
  npos<-ncol(pwm)
  ic<-numeric(length=npos)
  for (i in 1:npos) {
    ic[i]<- log(nrow(pwm), base=2) + sum(sapply(pwm[, i], function(x) {
      if (x > 0) { x*log2(x) } else { 0 }
    }))
  }
  ic
}


damage.ic<-function(pwm) {
  npos<-ncol(pwm[[1]])
  ic<- matrix(0, npos, length(pwm))

  for(i in 1:npos){
    mat <- numeric()
    for(j in 1:length(pwm)){
      mat <- cbind(mat, pwm[[j]][,i])
    }
    mat_clean <- mat[rowSums(mat) != 0,]
    ic[i,] <- pwm2ic(mat_clean)
  }

  return(ic)
}


###########   skeleton damage logo    ###################



###################################################################
#####################  Damage Logos  ##############################
###################################################################



###################  letter  A   #################################


letterA <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6)
  y <- c(0,10,10,0,0,3,3,0,0,4,7.5,4,4)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- c(rep(1,9),rep(2,4))
  }else{
    id <- c(rep(id,9),rep(id+1,4))
  }

  fill <- c("green","white")

  list(x=x,y=y,id=id,fill=fill)
}

#grid.newpage()

#grid.polygon(x=a_let$x, y=a_let$y, gp=gpar(fill=a_let$fill,
#             col="transparent"))

################   letter  T   ##################################

letterT <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c(0,10,10,6,6,4,4,0)
  y <- c(10,10,9,9,0,0,9,9)
  x <- 0.1*x
  y <- 0.1*y

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,8)
  }else{
    id <- rep(id,8)
  }

  fill <- "red"

  list(x=x,y=y,id=id,fill=fill)
}

#################### letter  C #####################################

letterC <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  x <- x.pos + wt*x
  y <- y.pos + ht*y

  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }

  fill <- "blue"

  list(x=x,y=y,id=id,fill=fill)
}


##################  letter G  ######################################

letterG <- function(x.pos,y.pos,ht,wt,id=NULL){
  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



  if (is.null(id)){
    id <- c(rep(1,length(x)),rep(2,length(x.add)))
  }else{
    id <- c(rep(id,length(x)),rep(id+1,length(x.add)))
  }

  x <- c(rev(x),x.add)
  y <- c(rev(y),y.add)

  x <- x.pos + wt*x
  y <- y.pos + ht*y


  fill <- c("orange","orange")

  list(x=x,y=y,id=id,fill=fill)

}

letterX <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c( 0, 0.4, 0, 0.2, 0.5, 0.8, 1, 0.6, 1, 0.8, 0.5, 0.2)
  x <- 0.05 + 0.90*x
  y <- c( 0, 0.5, 1, 1, 0.6, 1, 1, 0.5, 0, 0, 0.4, 0)

  if (is.null(id)){
    id <- rep(1,length(x))
  }else{
    id <- rep(id,length(x))
  }

  fill <- "pink"

  x <- x.pos + wt*x
  y <- y.pos + ht*y


  ll <- list("x"= x,
             "y"= y,
             "id" = id,
             "fill" = fill)
  return(ll)
}





################  letter C to T  ###############################

letter_C_to_T <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1 <- 0.4*x
  y1 <- 1*y


  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y


  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id + c(rep(1,length(x1)), rep(3, length(x2)),
                      rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), rep(3, length(x2)),
                 rep(4, length(x3)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("blue","red","grey80")

  list(x=x,y=y,id=id_pool,fill=fill)

}

###############  letter C to G  ###################################

letter_C_to_G <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1_pool <- 0.4*x
  y1_pool <- 1*y

  id1_pool <- rep(1,length(x1_pool))



  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



  id2_pool <- c(rep(2,length(x)),rep(3,length(x.add)))


  x2_pool <- 0.6 + 0.4*c(rev(x),x.add)
  y2_pool <- c(rev(y),y.add)


  x3_pool <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3_pool <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1_pool, x2_pool, x3_pool)
  ypool <- c(y1_pool, y2_pool, y3_pool)

  if(!is.null(id)){
    id_pool <- id + c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }else{
    id_pool <- c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill <- c("blue","orange", "orange", "grey80")

  list(x=x,y=y,id=id_pool,fill=fill)
}

##############  letter C to A  ####################################

letter_C_to_A <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x1 <- 0.4*x
  y1 <- 1*y


  x <- 0.1* (c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6))
  y <- 0.1*(c(0,10,10,0,0,3,3,0,0,4,7.5,4,4))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y

  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id +  c(rep(1,length(x1)), c(rep(2,9),rep(3,4)),
                       rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), c(rep(2,9),rep(3,4)),
                 rep(4, length(x3)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("blue","green", "white", "grey80")

  list(x=x,y=y,id=id_pool,fill=fill)
}

################  letter T to G  ################################


letter_T_to_G <- function(x.pos,y.pos,ht,wt,id=NULL){

  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1_pool <- 0.4*x
  y1_pool <- 1*y
  id1_pool <- rep(1,length(x1_pool))


  x3_pool <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3_pool <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)


  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))

  h1 <- max(y.l1)
  r1 <- max(x.l1)

  h1 <- 0.4
  x.add <- c(r1,0.5,0.5,r1-0.2,r1-0.2,r1,r1)
  y.add <- c(h1,h1,h1-0.1,h1-0.1,0,0,h1)



  id2_pool <- c(rep(2,length(x)),rep(3,length(x.add)))


  x2_pool <- 0.6 + 0.4*c(rev(x),x.add)
  y2_pool <- c(rev(y),y.add)

  if(!is.null(id)){
    id_pool <- id + c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }else{
    id_pool <- c(id1_pool, id2_pool, rep(4, length(x3_pool)))
  }

  xpool <- c(x1_pool, x2_pool, x3_pool)
  ypool <- c(y1_pool, y2_pool, y3_pool)

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("red","orange","orange", "grey80")

  list(x=x,y=y,id=id_pool,fill=fill)
}

################ letter  T to C  ################################

letter_T_to_C <- function(x.pos,y.pos,ht,wt,id=NULL){

  angle1 <- seq(0.3+pi/2,pi,length=100)
  angle2 <- seq(pi,1.5*pi,length=100)
  x.l1 <- 0.5 + 0.5*sin(angle1)
  y.l1 <- 0.5 + 0.5*cos(angle1)
  x.l2 <- 0.5 + 0.5*sin(angle2)
  y.l2 <- 0.5 + 0.5*cos(angle2)

  x.l <- c(x.l1,x.l2)
  y.l <- c(y.l1,y.l2)

  x <- c(x.l,rev(x.l))
  y <- c(y.l,1-rev(y.l))

  x.i1 <- 0.5 +0.35*sin(angle1)
  y.i1 <- 0.5 +0.35*cos(angle1)
  x.i1 <- x.i1[y.i1<=max(y.l1)]
  y.i1 <- y.i1[y.i1<=max(y.l1)]
  y.i1[1] <- max(y.l1)

  x.i2 <- 0.5 +0.35*sin(angle2)
  y.i2 <- 0.5 +0.35*cos(angle2)

  x.i <- c(x.i1,x.i2)
  y.i <- c(y.i1,y.i2)

  x1 <- c(x.i,rev(x.i))
  y1 <- c(y.i,1-rev(y.i))

  x <- c(x,rev(x1))
  y <- c(y,rev(y1))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y


  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1 <- 0.4*x
  y1 <- 1*y


  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id +c(rep(1,length(x1)), rep(3, length(x2)),
                     rep(4, length(x3)))
  }else{
    id_pool <- c(rep(1,length(x1)), rep(3, length(x2)),
                 rep(4, length(x3)))
  }

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  fill=c("blue","red","grey80")

  list(x=x,y=y,id=id_pool,fill=fill)

}

#################  letter T to A  ####################################

letter_T_to_A <- function(x.pos,y.pos,ht,wt,id=NULL){
  x <- c( 0.4, 0.4, 0, 0, 1, 1, 0.6, 0.6)
  y <- c( 0, 0.85, 0.85, 1, 1, 0.85, 0.85, 0)
  x1 <- 0.4*x
  y1 <- 1*y


  x <- 0.1* (c(0,4,6,10,8,6.8,3.2,2,0,3.6,5,6.4,3.6))
  y <- 0.1*(c(0,10,10,0,0,3,3,0,0,4,7.5,4,4))
  x2 <- 0.6 + 0.4*x
  y2 <- 1*y


  x3 <- c(0.42, 0.42, 0.55, 0.55, 0.60, 0.55, 0.55)
  y3 <- c(0.45, 0.55, 0.55, 0.60, 0.50, 0.40, 0.45)

  xpool <- c(x1, x2, x3)
  ypool <- c(y1, y2, y3)

  if(!is.null(id)){
    id_pool <- id + c( rep(3, length(x1)), c(rep(1,9),rep(2,4)),
                       rep(4, length(x3)))
  }else{
    id_pool <- c( rep(3, length(x1)), c(rep(1,9),rep(2,4)),
                  rep(4, length(x3)))
  }

  fill=c("green","white", "red", "grey80")

  x <- x.pos + wt*xpool
  y <- y.pos + ht*ypool

  list(x=x,y=y,id=id_pool,fill=fill)
}

################# add letters to logo plot #########################


addLetter <- function(letters,which,x.pos,y.pos,ht,wt){

  if (which == "A"){
    letter <- letterA(x.pos,y.pos,ht,wt)
  }else if (which == "C"){
    letter <- letterC(x.pos,y.pos,ht,wt)
  }else if (which == "G"){
    letter <- letterG(x.pos,y.pos,ht,wt)
  }else if (which == "X"){
    letter <- letterX(x.pos,y.pos,ht,wt)
  }else if (which == "T"){
    letter <- letterT(x.pos,y.pos,ht,wt)
  }else if (which == "C->T"){
    letter <- letter_C_to_T(x.pos,y.pos,ht,wt)
  }else if (which == "C->A"){
    letter <- letter_C_to_A(x.pos,y.pos,ht,wt)
  }else if (which == "C->G"){
    letter <- letter_C_to_G(x.pos,y.pos,ht,wt)
  }else if (which == "T->A"){
    letter <- letter_T_to_A(x.pos,y.pos,ht,wt)
  }else if (which == "T->G"){
    letter <- letter_T_to_G(x.pos,y.pos,ht,wt)
  }else if (which == "T->C"){
    letter <- letter_T_to_C(x.pos,y.pos,ht,wt)
  }else{
    stop("which must be one of A,C,G,T,X, C->T,
         C->G, C->A, T->A, T->C, T->G")
  }

  letters$x <- c(letters$x,letter$x)
  letters$y <- c(letters$y,letter$y)

  lastID <- ifelse(is.null(letters$id),0,max(letters$id))
  letters$id <- c(letters$id,lastID+letter$id)
  letters$fill <- c(letters$fill,letter$fill)
  letters
  }

damageLogo.pos.skeleton <- function(pwm,
                                probs,
                                ic,
                                max_pos,
                                ic.scale=TRUE,
                                xaxis=TRUE,
                                yaxis=TRUE,
                                xaxis_fontsize=10,
                                xlab_fontsize=15,
                                y_fontsize=15,
                                mut_width=2,
                                start=0.0001,
                                yscale_change=FALSE,
                                pop_name = NULL,
                                logoport_x=0.4,
                                logoport_y= 0.5,
                                logoport_width= 0.4,
                                logoport_height= 0.5,
                                lineport_x = 0.95,
                                lineport_y=0.85,
                                lineport_width=0.25,
                                lineport_height=0.3){

  if (class(pwm) == "data.frame"){
    pwm <- as.matrix(pwm)
  }else if (class(pwm) != "matrix"){
    stop("pwm must be of class matrix or data.frame")
  }

  if (any(abs(1 - apply(pwm,2,sum)) > 0.01))
    stop("Columns of PWM must add up to 1.0")

  chars <- c("A", "C", "G", "T", "X",
             "C->T", "C->A", "C->G",
             "T->A", "T->C", "T->G")

  letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
  npos <- ncol(pwm)


  if (ic.scale){
    if(yscale_change){
      if(max(ic)<1){ylim <- 1
      facs <- ic + 1 - max(ic)}
      if(max(ic)>1){ylim <- 2
      facs <- ic}
    }else{
      ylim <- 2
      facs <- ic
    }
    ylab <- "Information content"
  }else{
    ylim <- 1
    ylab <- "Probability"
    facs <- rep(1, npos)
  }

  wt <- c(rep(1, floor(npos/2)),mut_width,rep(1, floor(npos/2)))
  flanked_coord <- c((-floor(npos/2)):(-1), 0, 1:floor(npos/2))
  x.pos <- 0
  for (j in 1:npos){

    column <- pwm[,j]
    hts <- 0.99*column*facs[j]
    letterOrder <- order(hts)

    y.pos <- 0
    for (i in 1:length(chars)){
      letter <- chars[letterOrder[i]]
      ht <- hts[letterOrder[i]]
      if (ht>0) letters <- addLetter(letters,letter,x.pos,y.pos,ht,wt[j])
      y.pos <- y.pos + ht + start
    }
    x.pos <- x.pos + wt[j]

  }

#  bottomMargin = ifelse(xaxis, 2 + xaxis_fontsize/3.5, 3)
#  leftMargin = ifelse(yaxis, 2 + y_fontsize/3.5, 3)

  grid.newpage()
  vp <- viewport(x=logoport_x,y=logoport_y,width=logoport_width, height=logoport_height)
  pushViewport(vp)
  #pushViewport(plotViewport(c(bottomMargin,leftMargin,max(ylim)+1.5,max(ylim)+ 1.5)))
  # pushViewport(viewport(layout = grid.layout(2, 2),
  #              x = bottomMargin,
  #              y = leftMargin,
  #              width = max(xlim/2)+0.5,
  #              height = max(ylim/2)+0.5))
  pushViewport(dataViewport(0:ncol(pwm),0:ylim,name="vp1"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id, gp=gpar(fill=letters$fill,col="transparent"))
  grid.polygon(x=unit(letters$x,"native"), y=unit(letters$y,"native"),
               id=letters$id,
               gp=gpar(fill=letters$fill,col="transparent"))

  xlim <- cumsum(wt) - wt/2;
  # xlim <- c(wt[1]/2, wt[1] + wt[2]/2, wt[1]+wt[2]+wt[3]/2, wt[1]+wt[2]+wt[3], 5.5)
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  ylim_scale <- seq(0, ylim, length.out=6);
  ic_lim_scale <- seq(0, max(ic), length.out=6)

  for(n in 1:length(xlim)){
    grid.lines(x = unit(low_xlim[n], "native"),
               y = unit(c(0, ylim), "native"),
               gp=gpar(col="grey80"))
  }

  if(is.null(pop_name)){
    grid.text("Damage logo plot", y = unit(1, "npc") + unit(1.5, "lines"),
              gp = gpar(fontsize = 16))
  }else{
    grid.text(paste0("Logo plot for ", pop_name), y = unit(1, "npc") + unit(1.5, "lines"),
              gp = gpar(fontsize = 16))
  }

  if (xaxis){
    grid.xaxis(at=xlim,
               label=(c(paste0("\n left \n flank \n", flanked_coord[1:floor(npos/2)]),
                        "\n mutation",
                        paste0("\n right \n flank \n", -flanked_coord[(1:floor(npos/2))])
               )),
               gp=gpar(fontsize=xaxis_fontsize))
    grid.text("Position",y=unit(-3,"lines"),
              gp=gpar(fontsize=xlab_fontsize))
  }
  if (yaxis){
    if(yscale_change==TRUE){
      grid.yaxis(at = ylim_scale,
                 label = round(ic_lim_scale,1),
                 gp=gpar(fontsize=y_fontsize))
    }else{
      grid.yaxis(gp=gpar(fontsize=y_fontsize))
    }
    grid.text(ylab,x=unit(-3,"lines"),rot=90,
              gp=gpar(fontsize=y_fontsize))
  }

  upViewport(0)
  vp2 <- viewport(x=lineport_x,y=lineport_y,width=lineport_width, height=lineport_height, just=c("right","top"))
  pushViewport(vp2)
  par(plt = gridPLT(), new=TRUE)
  plot_graph(probs, max_pos)
  popViewport(0)
  par(ask=FALSE)
}


plot_graph <- function(probs, max_pos, col="red",
                       cex=unit(1, "npc"), pch=unit(16,"npc"),
                       xlab="read position", ylab="mutation probability",
                       main="",
                       cex.axis=unit(0.75, "npc"),
                       cex.main=unit(1, "npc")){
  if (length(probs) != max_pos){
    stop(cat('probability vector must be of length ', max_pos))
  }
  plot(1:max_pos, probs/max(probs), xlim = c(0, max_pos), ylim=c(0,1), xlab = xlab, ylab = ylab,
       type = "b", xaxt = "n", yaxt = "n", cex = cex, pch=pch, col=col, main=main,
       cex.main=cex.main)
  axis(side = 1, at = round(seq(0, max_pos, length.out=4),1), cex.axis = cex.axis, lwd.ticks = 2)
  ylimit <- c(0.0, 0.5, 1.0)*max(probs)
  axis(side = 2, at = c(0.0, 0.5, 1.0), labels = round(ylimit,1), cex.axis = cex.axis, lwd.ticks=2)
}


