#' Derives the beta-binomial parameters to target a mean and standard deviation
#' around the selected number of predictors
#' @export
#' @title Derives the beta-binomial parameters for target mean and standard deviation
#' @name GetBetaBinomParams
#' @param n The total number of covariates in the analysis
#' @param mu Desired mean number of selected covariates
#' @param sd Desired standard deviation around the selected covariates
#' @return A list containing the estimated beta-binomial parameters `a' and `b'. 
#' @author Paul Newcombe
GetBetaBinomParams <- function(n, mu, sd) {
  # Convert mean and variance (2nd moment about the mean) to moments about the origin
  m1 <- mu
  m2 <- sd^2 + mu^2
  # Below is expression in terms of moments about the origin
  a <- (n*m1 - m2)/(n*(m2/m1 -m1 -1) + m1)
  b <- (n - m1)*(n - m2/m2) / ( n*(m2/m1 -m1 -1) + m1)
  return(params = list(a = a, b = b))
}

#' Derives the beta-binomial mean and standard deviation
#' @export
#' @title Derives the beta-binomial mean and standard deviation
#' @name GetBetaBinomMuSd
#' @param n The total number of covariates in the analysis
#' @param a Beta-binomial `a' or alpha parameter
#' @param b Beta-binomial `b' or beta parameter
#' @return A list containing the beta-binomial mean and standard deviation. 
#' @author Paul Newcombe
GetBetaBinomMuSd <- function(n, a, b) {
  pi <- a/(a+b)
  ro <- 1/(a+b+1)
  mu <- n*pi
  sd <- sqrt(n*pi*(1-pi)*(1+(n-1)*ro))
  return(params = list(mu = mu, sd = sd))
}
