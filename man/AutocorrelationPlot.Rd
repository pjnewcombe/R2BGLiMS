\name{AutocorrelationPlot}
\alias{AutocorrelationPlot}
\title{Autocorrelation plot for model selection. Model inclusion of each predictor is
represented across all iterations as a black and white "heatmap".}
\usage{
  AutocorrelationPlot(results, vars.to.include = NULL,
    plot.file = NULL, include.var.names = TRUE,
    var.dictionary = NULL, include.post.probs = TRUE)
}
\arguments{
  \item{plot.file}{Path to a PDF to print plot to.}

  \item{include.var.names}{Include variable names
  perpendicularly below the x-axis. Turn off if it looks
  too crowded.}

  \item{include.post.probs}{Include a barplot of posterior
  probabilities at the top.}

  \item{results}{Reversible Jump results object from
  running \code{\link{R2MHRJ}}.}

  \item{vars.to.include}{Optional character vector
  specifying a subset of variables to restrict output to.}

  \item{var.dictionary}{Alternative covariate names. This
  must be a character vector, named wtih the current
  covariate names. This can consist of expression
  functions, e.g. to achieve italics or superscripts when
  passed to a plot function.}
}
\value{
  NA
}
\description{
  Generates an autocorrelation plot from a Java MCMC
  results class object
}
\examples{
library(R2BGLiMS)

### --- Logistic regression with two model space prior components
# Load results from the biopsy logistic regression analysis (see R2MHRJ examples)
data(biopsyResults)
AutocorrelationPlot(biopsyResults)
}
\author{
  Paul Newcombe
}

