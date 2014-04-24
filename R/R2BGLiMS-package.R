#' Provides an R interface to a Java package for fitting Bayesian GLMs. Approximate posterior
#' samples are drawn using an MCMC sampler with a (Reversible Jump) Metropolis-Hastings
#' acceptance ratio. Currently, Logistic and Weibull regression models are available.
#' - Allows models to be fitted with random intercepts.
#' - Allows model selection analyses via Reversible Jump.
#' 
#' \tabular{ll}{
#' Package: \tab R2BGLiMS\cr
#' Type: \tab Package\cr
#' Version: \tab 0.1\cr
#' Date: \tab 2013-11-21\cr
#' License: \tab GPL (>=2)\cr
#' LazyLoad: \tab yes\cr
#' }
#' 
#' @name R2BGLiMS-package
#' @docType package
#' @title R interface to the java package "BGLiMS" (Bayesian Generalised Linear Model Selection).
#' @author Paul Newcombe \email{paul.newcombe@@mrc-bsu.cam.ac.uk}
#' @keywords package
NULL