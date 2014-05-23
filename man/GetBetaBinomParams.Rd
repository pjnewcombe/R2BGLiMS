\name{GetBetaBinomParams}
\alias{GetBetaBinomParams}
\title{Derives the beta-binomial parameters for target mean and standard deviation}
\usage{
  GetBetaBinomParams(n, mu, sd)
<<<<<<< HEAD
=======

  GetBetaBinomMuSd(n, a, b)
>>>>>>> 41f1c0a516fdd2dd4e0b463bb9aaf7829e877cb6
}
\arguments{
  \item{n}{The total number of covariates in the analysis}

  \item{mu}{Desired mean number of selected covariates}

  \item{sd}{Desired standard deviation around the selected
  covariates}
<<<<<<< HEAD
=======

  \item{n}{The total number of covariates in the analysis}

  \item{a}{Beta-binomial `a' or alpha parameter}

  \item{b}{Beta-binomial `b' or beta parameter}
>>>>>>> 41f1c0a516fdd2dd4e0b463bb9aaf7829e877cb6
}
\value{
  A list containing the estimated beta-binomial parameters
  `a' and `b'.
<<<<<<< HEAD
=======

  A list containing the beta-binomial mean and standard
  deviation.
>>>>>>> 41f1c0a516fdd2dd4e0b463bb9aaf7829e877cb6
}
\description{
  Derives the beta-binomial parameters to target a mean and
  standard deviation around the selected number of
  predictors
<<<<<<< HEAD
}
\author{
  Paul Newcombe
=======

  Derives the beta-binomial mean and standard deviation
}
\author{
  Paul Newcombe

  Paul Newcombe
>>>>>>> 41f1c0a516fdd2dd4e0b463bb9aaf7829e877cb6
}

