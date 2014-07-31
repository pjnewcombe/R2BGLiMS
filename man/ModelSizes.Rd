\name{ModelSizes}
\alias{ModelSizes}
\title{Posterior model size/dimension summary}
\usage{
  ModelSizes(results)
}
\arguments{
  \item{results}{Reversible Jump results object from
  running \code{\link{R2BGLiMS}}.}
}
\value{
  Table of posterior probabilities and Bayes Factors for
  different model sizes
}
\description{
  Calulcates posterior probabiltiies and Bayes Factors for
  different model sizes, from a Reversible Jump Results
  object
}
\examples{
library(R2BGLiMS)

### --- Logistic regression with two model space prior components
# Load results from the biopsy logistic regression analysis (see R2MHRJ examples)
data(biopsyResults)
ModelSizes(biopsyResults)
}
\author{
  Paul Newcombe
}

