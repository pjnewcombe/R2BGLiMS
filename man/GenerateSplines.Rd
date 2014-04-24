\name{GenerateSplines}
\alias{GenerateSplines}
\title{Generates splines.}
\usage{
  GenerateSplines(data, which.vars = NULL, n.knots = 6)
}
\arguments{
  \item{data}{Data set to append splines to}

  \item{which.vars}{Which variables to generate splines
  for. If left at default NULL}

  \item{n.knots}{How many splines to generate for each
  variable (knots will be placed at quantiles). NOTE: if a
  covariate does not have enough unique quantiles it will
  be excluded from the spline modelling, and returned in
  it's original form.}
}
\value{
  spline.data A dataset containing the splines. Can be
  merged to the original data using cbind, and the names of
  the spline covariates can be obtained using colnames.
}
\description{
  Generates splines to append to original dataset, for
  spline modelling. Knots are created at quantiles, and
  binary indicators of being greater than or equal to the
  knot are returned. NOTE: if a covariate does not have
  enough unique quantiles it will be excluded from the
  spline modelling, and returned in it's original form.
}
\author{
  Paul Newcombe
}

