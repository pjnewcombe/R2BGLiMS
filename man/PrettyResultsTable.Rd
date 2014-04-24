\name{PrettyResultsTable}
\alias{PrettyResultsTable}
\title{A pretty summary results table for reports - note that entries are formatted as character strings so,
for processing results use \code{\link{ResultsTable}}.}
\usage{
  PrettyResultsTable(results, round.digits.betas = 2,
    round.digits.postprob = 2, round.digits.bf = 1,
    measure = "OR", vars.to.include = NULL,
    var.dictionary = NULL, normalised.sds = NULL,
    post.prob.cut = NULL, abbreviated.names = FALSE,
    remove.col.white.spaces = FALSE)
}
\arguments{
  \item{round.digits.betas}{Number of decimal places to
  include for effect estimates. (Default is 2)}

  \item{round.digits.postprob}{Number of decimal places to
  include for posterior probabilities. (Default is 2)}

  \item{round.digits.bf}{Number of decimal places to
  include for Bayes factors. (Default is 1)}

  \item{measure}{Abbreviated measure to include in column
  titles, e.g. "HR". (Default is "OR")}

  \item{post.prob.cut}{Optionally can provide a posterior
  probability threshold for including predictors in the
  table for which model selection has been performed for.}

  \item{abbreviated.names}{Use abbreviated names for
  Posterior Probability and Bayes Factors columns? (Default
  is FALSE)}

  \item{remove.col.white.spaces}{Should the column titles
  have no white spaces? (Default is FALSE)}

  \item{results}{Reversible Jump results object from
  running \code{\link{R2MHRJ}}.}

  \item{vars.to.include}{Optional character vector
  specifying a subset of variables to restrict output to.}

  \item{var.dictionary}{Alternative covariate names. This
  must be a character vector, named wtih the current
  covariate names. This can consist of expression
  functions, e.g. to achieve italics or superscripts when
  passed to a plot function.}

  \item{normalised.sds}{If variables were normalised (e.g.
  to facilitate use of a common unknown prior on their
  effects), a named vector of their original sds can be
  provided such that estimates returned may be interpreted
  on the original scale.}
}
\value{
  Prints a nice summary results table
}
\description{
  Creates a `pretty' summary results table from a
  Reversible Jump results object
}
\examples{
library(R2BGLiMS)

### --- Logistic regression with two model space prior components
# Load results from the VA survival analysis (see R2MHRJ examples)
data(VAResults)
PrettyResultsTable(VAResults)
}
\author{
  Paul Newcombe
}

