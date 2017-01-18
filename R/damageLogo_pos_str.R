#' @title Builds damage logo plots along with mutational profile across read and the strand break composition for the clusters
#'
#' @description Damage Logo plots for each cluster representing the substitution frequency of the 6 substitutional
#' patterns (adjusting for strand bias) and the flanking bases arranged as per relative frequency on either side of
#' the substitution along with the probability of the cluster mutational profile along the read and the strand break composition.
#'
#' @param theta_pool The theta matrix obtained from running the grade of membership model that stores for each cluster, the
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
#' @param renyi_alpha The entropy scale for the Renyi entropy on the flanking bases and mutations.
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
#' @import grid
#' @import gridBase
#' @import dplyr
#' @import plyr
#'
#' @export



damageLogo_pos_str <- function(theta_pool,
                               sig_names = NULL,
                               ic.scale=TRUE,
                               max_pos = 20,
                               flanking_bases=2,
                               yscale_change = TRUE,
                               xaxis=TRUE,
                               yaxis=TRUE,
                               xaxis_fontsize=5,
                               xlab_fontsize=10,
                               y_fontsize=10,
                               mut_width=2,
                               start=0.0001,
                               renyi_alpha = 1,
                               pop_names=paste0("Cluster ",1:dim(theta_pool)[2]),
                               logoport_x = 0.24,
                               logoport_y= 0.50,
                               logoport_width= 0.30,
                               logoport_height= 0.40,
                               lineport_x = 0.70,
                               lineport_y=0.50,
                               lineport_width=0.20,
                               lineport_height=0.25,
                               pieport_x = 1.18,
                               pieport_y = 5,
                               pieport_width=0.70,
                               pieport_height=0.6,
                               barport_x = 0.58,
                               barport_y=0.65,
                               barport_width=0.25,
                               barport_height=0.25){

  signature_set <- rownames(theta_pool)

  indices_left <- grep("left", signature_set)
  indices_right <- grep("right", signature_set)

  breaks_theta <- theta_pool[c(indices_left, indices_right),]

  theta_new <- theta_pool[-c(indices_left, indices_right),]
  theta_new<- apply(theta_new, 2, function(x) return(x/sum(x)))

  indices_minus <- grep("_-_", signature_set)
  strand_theta <- data.frame("minus" = colSums(theta_new[indices_minus,]),
                             "plus" = colSums(theta_new[-indices_minus,]))
  strand_theta_vec <- data.frame(strand_theta[1,])

  signature_set_new <- rownames(theta_new)
  signature_set_split <- do.call(rbind, lapply(signature_set_new, function(x) {
    y = strsplit(as.character(x), split="")[[1]][-c((4+2*flanking_bases + 1):(4+2*flanking_bases + 2))]
    return(paste0(y, collapse=""))
  }))

  theta <- tbl_df(data.frame(theta_new)) %>% mutate(sig = signature_set_split) %>% group_by(sig) %>% dplyr::summarise_each(funs(sum)) %>% as.data.frame()
  rownames(theta) <-  theta[,1]
  theta <- theta[,-1]

  if(is.null(sig_names))
    sig_names <- rownames(theta)

  prob_mutation <- filter_signatures_only_location(t(theta_new), max_pos = max_pos, flanking_bases = flanking_bases)
  prob_mutation <- t(apply(prob_mutation, 1, function(x) {
    y <- x[!is.na(x)];
    return(y/sum(y))
  }))
  max_prob <- max(prob_mutation);

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

  ic <- damage.ic(prop_patterns_list, alpha=renyi_alpha)

  grob_list <- list()
  for(l in 1:length(prop_patterns_list)){
    damageLogo.pos.str.skeleton(pwm = prop_patterns_list[[l]],
                                probs = prob_mutation[l,],
                                breaks_theta_vec = breaks_theta[,l],
                                strand_theta_vec = strand_theta[l,],
                                ic = ic[,l],
                                max_pos = max_pos,
                                max_prob = max_prob,
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
                                lineport_height=lineport_height,
                                pieport_x = pieport_x,
                                pieport_y = pieport_y,
                                pieport_width=pieport_width,
                                pieport_height=pieport_height,
                                barport_x = barport_x,
                                barport_y = barport_y,
                                barport_width = barport_width,
                                barport_height = barport_height)

  }
}

damageLogo.pos.str.skeleton <- function(pwm,
                                        probs,
                                        breaks_theta_vec,
                                        strand_theta_vec,
                                        ic,
                                        max_pos,
                                        max_prob,
                                        ic.scale=TRUE,
                                        xaxis=TRUE,
                                        yaxis=TRUE,
                                        xaxis_fontsize=10,
                                        xlab_fontsize=15,
                                        y_fontsize=15,
                                        mut_width=2,
                                        start=0.0001,
                                        yscale_change=TRUE,
                                        pop_name = NULL,
                                        logoport_x=0.3,
                                        logoport_y= 0.5,
                                        logoport_width= 0.3,
                                        logoport_height= 0.8,
                                        lineport_x = 0.6,
                                        lineport_y=0.25,
                                        lineport_width=0.30,
                                        lineport_height=0.30,
                                        pieport_x = 1,
                                        pieport_y = 0.50,
                                        pieport_width=0.40,
                                        pieport_height=0.40,
                                        barport_x = 0.5,
                                        barport_y=0.80,
                                        barport_width=0.25,
                                        barport_height=0.25){

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
      #facs <- ic + 1 - max(ic)
      facs <- ic/max(ic)
      }
      if(max(ic)>1){ylim <- 2
      facs <- ic}
    }else{
      ylim <- ceiling(max(ic))
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
    hts <- as.numeric(0.99*column*facs[j])
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

  xlim <- cumsum(wt) - wt/2;
  # xlim <- c(wt[1]/2, wt[1] + wt[2]/2, wt[1]+wt[2]+wt[3]/2, wt[1]+wt[2]+wt[3], 5.5)
  low_xlim <- c(xlim - 0.5*wt, xlim[length(xlim)]+0.5*wt[length(xlim)])
  ylim_scale <- seq(0, ylim, length.out=6);
  if(ic.scale){
    ic_lim_scale <- seq(0, max(ic), length.out=6)
  }else{
    ic_lim_scale <- ylim_scale
  }
  if(ic.scale){
    if(yscale_change){
      if(ylim  > 1){
        letters$y <- letters$y*(ylim/max(ic));
      }
    }}


  #  bottomMargin = ifelse(xaxis, 2 + xaxis_fontsize/3.5, 3)
  #  leftMargin = ifelse(yaxis, 2 + y_fontsize/3.5, 3)

  grid.newpage()
  vp <- viewport(x=logoport_x, y=logoport_y, width=logoport_width, height=logoport_height)
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
    grid.text("Logo plot", y = unit(1, "npc") + unit(1.5, "lines"),
              gp = gpar(fontsize = 16))
  }else{
    grid.text(paste0(pop_name), x = unit(1.3, "npc"), y = unit(12, "lines"),
              gp = gpar(fontsize = 20, col="black"))
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
                 label = round(ic_lim_scale,2),
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
  plot_graph(probs, max_pos=max_pos, max_prob=max_prob)
  upViewport(0)
  vps <- baseViewports()
  #pushViewport(vps$figure)
  # vp3 <- viewport(x=pieport_x, y=pieport_y, width=pieport_width, height=pieport_height, just=c("right","bottom"))
 # pushViewport(vp3)
 # par(plt = gridPLT(), new=TRUE)
  vp3 <- viewport(width = pieport_width, height = pieport_height, x = pieport_x,
                 y = unit(pieport_y, "lines"), just = c("right","bottom"))
  p = plot_pie(breaks_theta_vec = breaks_theta_vec)
  print(p, vp = vp3)

  vp4 <- viewport(width = barport_width, height = barport_height, x = barport_x, y = barport_y)
  p = plot_bar(strand_theta_vec = strand_theta_vec)
  print(p, vp = vp4)
  #plot(1:4, col="red", pch=20)
  popViewport(0)
  par(ask=FALSE)
}

plot_pie <- function(breaks_theta_vec){

  names = names(breaks_theta_vec)
  bases <- c(as.character(sapply(names, function(x) return(strsplit(x, "[_]")[[1]][1]))))
  shuffle <- c(match(c("A", "G", "C", "T"), bases[1:4]), 4+ match(c("A", "G", "C", "T"), bases[5:8]))
  strand <- c(as.character(sapply(names, function(x) return(strsplit(x, "[_]")[[1]][2]))))
  strand <- plyr::revalue(factor(strand), c("left"="5' strand break", "right"="3' strand break"))
  sum1 = rep(tapply(breaks_theta_vec, strand, sum), each=4)
  breaks_theta_vec_2 <- breaks_theta_vec/sum1
  df <- data.frame("value" = bases[shuffle], "category"=strand[shuffle], "percentage" = breaks_theta_vec_2[shuffle])

  # get counts and melt it
  data.m = df
  names(data.m)[3] = "percentage"

  # calculate percentage:
  m1 = plyr::ddply(data.m, plyr::.(category), plyr::summarize, ratio=percentage/sum(percentage))

  #order data frame (needed to comply with percentage column):
  #m2 = data.m[order(data.m$category),]
  m2 = data.m

  # combine them:
  mydf = data.frame(m2,ratio=m1$ratio)

  # get positions of percentage labels:
  mydf = plyr::ddply(mydf, .(category), transform, position = 1 - cumsum(percentage) + 0.5*percentage)
  bases = factor(mydf$value, levels=c("A", "G", "C", "T"))
  # create bar plot
  pie = ggplot(mydf, aes(x = factor(1), y = percentage, fill = bases)) +
    geom_bar(stat = "identity", width = 1) +
    facet_wrap(~category, nrow=2) +
    labs(x = "") + labs(y="") + labs(fill="")
   # labs(y = "strand breaks composition", size=2)+
   # theme_nothing() + labs(x = NULL, y = NULL)
   # labs(title = "Strand breaks composition", size=5)

  # make a pie
  pie = pie + coord_polar(theta = "y") +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank(),
          panel.border = element_rect(colour = "white")) +
    theme(plot.margin = unit(c(0.01,0.01,0.01,0.01), "cm"))

  # add labels
  pie + geom_text(aes(label = sprintf("%1.2f%%", 100*ratio), y = position), size=2)

}


plot_bar <- function(strand_theta_vec){
  df <- data.frame("value" = as.numeric(t(strand_theta_vec)), "strand"=plyr::revalue(factor(names(strand_theta_vec)), c("plus" = "+",  "minus" = "-")))
  ggplot(df, aes(x=strand, y=value, fill=c("green", "lightpink")))  +
    geom_bar(stat='identity', width=0.8, position = position_dodge(width=0.9), colour = 'black') +
    xlab("") + ylim(0,1) + ylab("") + theme(legend.position="none") +
  #  geom_text(aes(x=strand, y=value+0.1, label=value)) +
    geom_text(aes(x=strand, y=value*0.5, label=as.character(strand)), colour="black", size=5)+
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid  = element_blank(),
          panel.border = element_rect(colour = "white")) +
    coord_flip()
}

