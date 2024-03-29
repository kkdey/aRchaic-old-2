% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/damageLogo_5.R
\name{damageLogo_five}
\alias{damageLogo_five}
\title{Builds damage logo plots along with mutational profile across read and the strand break composition for the clusters}
\usage{
damageLogo_five(theta_pool, sig_names = NULL, ic.scale = TRUE,
  max_pos = 20, flanking_bases = 1, yscale_change = TRUE, xaxis = TRUE,
  yaxis = TRUE, xaxis_fontsize = 5, xlab_fontsize = 10, y_fontsize = 10,
  mut_width = 2, start = 1e-04, renyi_alpha = 1,
  pop_names = paste0("Cluster ", 1:dim(theta_pool)[2]), logoport_x = 0.24,
  logoport_y = 0.5, logoport_width = 0.3, logoport_height = 0.4,
  lineport_x = 0.7, lineport_y = 0.5, lineport_width = 0.2,
  lineport_height = 0.25, stackbarport_x = 1, stackbarport_y = 5,
  stackbarport_width = 0.3, stackbarport_height = 0.7, barport_x = 0.58,
  barport_y = 0.65, barport_width = 0.25, barport_height = 0.25)
}
\arguments{
\item{theta_pool}{The theta matrix obtained from running the grade of membership model that stores for each cluster, the
probability distribution over all the mutational signatures.}

\item{sig_names}{The mutational signature names. Defaults to the rownames of the theta matrix above.}

\item{ic.scale}{A binary variable indicating whether the height of the bars for substitution and flanking bases should be
adjusted by the information criterion.}

\item{max_pos}{The maximum distance from the end of the read upto which mutations are considered.}

\item{flanking_bases}{The number of flanking bases of the mutational signature.}

\item{yscale_change}{A binary variable indicating whether the Y axis scale should be adjusted based on the size of the
logos, defaults to TRUE.}

\item{xaxis}{A binary indicating whether the X axis of the logo plot should be shown}

\item{yaxis}{A binary indicating whether the Y axis of the logo plot should be shown}

\item{xaxis_fontsize}{The fontsize of the X axis ticks.}

\item{xlab_fontsize}{The fontsize of the X axis labels.}

\item{y_fontsize}{The fontsize of the Y axis ticks.}

\item{mut_width}{Thw width of the bar for the mutation at the center.}

\item{start}{The starting point of the stacking of logos on the Y axis. Should be close to 0, defau;ts to 0.0001.}

\item{renyi_alpha}{The entropy scale for the Renyi entropy on the flanking bases and mutations.}

\item{pop_names}{The title of the plot. Defaults to the cluster labels.}

\item{logoport_x}{the X-axis position of the plot window for the logo plot}

\item{logoport_y}{the Y-axis position of the plot window for the logo plot}

\item{logoport_width}{the width of the plot window for the logo plot}

\item{logoport_height}{the width of the plot window for the logo plot}

\item{lineport_x}{the X-axis position of the plot window for the mutational profile line plot.}

\item{lineport_y}{the Y-axis position of the plot window for the mutational profile line plot.}

\item{lineport_width}{the width of the plot window for the mutational profile line plot.}

\item{lineport_height}{the width of the plot window for the mutational profile line plot.}
}
\value{
Returns logo plots for each cluster
}
\description{
Damage Logo plots for each cluster representing the substitution frequency of the 6 substitutional
patterns (adjusting for strand bias) and the flanking bases arranged as per relative frequency on either side of
the substitution along with the probability of the cluster mutational profile along the read and the strand break composition.
}

