\name{GetBetaBinomParams}
\alias{GetBetaBinomParams}
\title{Derives the beta-binomial parameters for target mean and standard deviation}
\usage{
  GetBetaBinomParams(n, mu, sd)

  GetBetaBinomMuSd(n, a, b)
}
\arguments{
  \item{n}{The total number of covariates in the analysis}

  \item{mu}{Desired mean number of selected covariates}

  \item{sd}{Desired standard deviation around the selected
  covariates}

  \item{n}{The total number of covariates in the analysis}

  \item{a}{Beta-binomial `a' or alpha parameter}

  \item{b}{Beta-binomial `b' or beta parameter}
}
\value{
  A list containing the estimated beta-binomial parameters
  `a' and `b'.

  A list containing the beta-binomial mean and standard
  deviation.
}
\description{
  Derives the beta-binomial parameters to target a mean and
  standard deviation around the selected number of
  predictors

  Derives the beta-binomial mean and standard deviation
}
\author{
  Paul Newcombe

  Paul Newcombe
}

