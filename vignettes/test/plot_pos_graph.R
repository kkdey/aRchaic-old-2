
plot_graph <- function(probs, col="red",
                       cex=1.25, pch=16,
                       xlab="", ylab="",
                       cex.axis=1.75){
  if (length(probs) != 23){
    stop('probability vector must be of length 23')
  }
  plot(c(2:20, 22.5, 30, 45), probs, xlim = c(0, 45), xlab = xlab, ylab = ylab,
       type = "b", xaxt = "n", yaxt = "n", ylim  = c(0,1), cex = cex, pch=pch, col=col)
  axis(side = 1, at = c(0, 10, 20, 30, 45), cex.axis = cex.axis, lwd.ticks = 2)
  axis(side = 2, at = c(0.0, 0.5, 1.0), labels = c(0.0, 0.5, 1.0), cex.axis = cex.axis, lwd.ticks=2)
}

probs = sort(runif(23), decreasing = TRUE)
plot_graph(probs)
