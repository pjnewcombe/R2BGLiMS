\name{ChainPlots}
\alias{ChainPlots}
\title{Coefficient history trace plots from a Reversible Jump results object}
\usage{
  ChainPlots(results, vars.to.include = NULL,
    var.dictionary = NULL, predictors.only = FALSE,
    plot.file = NULL, add.to.title = NULL,
    par.mfrow = NULL, add.panel.labels = FALSE)
}
\arguments{
  \item{predictors.only}{If true only chain plots for
  predictors are plotted (and intercept, likelihood etc.
  are excluded). Default is FALSE.}

  \item{plot.file}{PDF path}

  \item{add.to.title}{Optional character string to add to
  the start of each plot title}

  \item{add.panel.labels}{Add lower case letter labels to
  panel titles. Default FALSE. Note: Does not work when
  providing a var.dictionary}

  \item{par.mfrow}{Allow passing the mfrow layout vector to
  par. Default par(mfrow=c(5,2))}

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
  Generates a PDF of chain plots for all variables from a
  Reversible Jump results object
}
\examples{
library(R2BGLiMS)

### --- Logistic regression with two model space prior components
# Load results from the biopsy logistic regression analysis (see R2MHRJ examples)
data(biopsyResults)
ChainPlots(biopsyResults)
}
\author{
  Paul Newcombe
}

