#' This function can be used to generate a posterior sample of effects under
#' the Normal-Inverse-Gamma conjugate model, for a particular combination of
#' covariates (i.e. conditional on a fixed model). For use with JAM specify
#' xTx and z instead of data and outcome.var.
#' Note to self: adapted from
#' Michael Jordan's lecture notes "Bayes Factors, g-priors, and Model 
#' Selection for Regression". Conugate expression is given under the 1/sigma prior.
#' For the InverseGamma(a, b) prior as below, simply add a and b to the respective 
#' hyperparameters in equation (6) of Michael Jordan's notes.
#' 
#' @export
#' @title Generate conjugate posterior sample of coefficients for a particular linear model.
#' @name NormInvGamPosteriorSample
#' @inheritParams R2BGLiMS
#' @param sigma2_invGamma_a Inverse-Gamma parameter a for the residual
#' variance. Not specifying means the value in default.arguments is used (type "data(DefaultArguments)").
#' @param sigma2_invGamma_b Inverse-Gamma parameter b for the residual
#' variance. Not specifying means the value in default.arguments is used (type "data(DefaultArguments)").
#' @param model Vector of covariate names to include in the model. Do not include confounders here - they
#' should be specified with the confounders argument.
#' @param tau Value to use for the g-prior sparsity parameter (tau*sigma^2 parameterisation).
#' @param n.samples Number of posterior samples to draw.
#' 
#' @return The posterior sample as a matrix. 
#' Rows are different posterior samples, and columns are parameters.
#' 
#' @author Paul Newcombe
NormInvGamPosteriorSample <- function(
  data=NULL,
  outcome.var=NULL,
  confounders=NULL,
  model=NULL,
  tau=NULL,
  n.samples=1000,
  xTx=NULL,
  z=NULL,
  sigma2_invGamma_a=NULL,
  sigma2_invGamma_b=NULL
) {
  library(MASS)
  # Get InverseGamma parameters from default arguments if not specified
  if (is.null(sigma2_invGamma_a)) {
    data(DefaultArguments)
    sigma2_invGamma_a <- default.arguments$GaussianResidualVarianceInvGammaPrior_a
  }
  if (is.null(sigma2_invGamma_b)) {
    data(DefaultArguments)
    sigma2_invGamma_b <- default.arguments$GaussianResidualVarianceInvGammaPrior_b
  }  
  
  # Deal with confounders
  if (!is.null(confounders)) {
    if (length(confounders %in% model)>0) {
      model <- model[!model%in%confounders]
    }
    cat("Obtaining confounder adjusted residuals...\n")
    data <- .ConfounderAdjustedResiduals(data, outcome.var, confounders)    
  }
  
  # Data setup
  if (!is.null(xTx)) {
    # Data setup for JAM
    L <- chol(xTx[[1]])
    X <- L[,model]
    y <- solve(t(L))%*%z
    n <- length(z)
  } else {
    X <- as.matrix(data[,model])
    y <- data[,outcome.var]
    n <- nrow(X)    
  }
  
  # Conjugate multivariate normal setup
  xTxInv <- solve(t(X)%*%X)
  beta.hat <- (xTxInv%*%t(X)) %*% y
  mvn.mean <- beta.hat*tau/(tau+1)
  mvn.sigma.multiplier <- xTxInv*tau/(1+tau)
  
  # Conjugate inverse-gamma setup
  s.squared <- t(y-X%*%beta.hat)%*%(y-X%*%beta.hat)
  InverseGamma_a <- sigma2_invGamma_a + n/2
  InverseGamma_b <- sigma2_invGamma_b + s.squared/2 + (t(beta.hat) %*% (t(X)%*%X) %*% beta.hat)[1,1]/(2*(tau+1))    
    
  # Draw from conjugate Normal-Inverse-Gamma using a gibbs sampler, i.e.
  # 1) Draw sigma from Inverse-Gamma
  # 2) Draw beta from MVN conditional on sigma
  posterior.samples <- NULL
  for (i in 1:n.samples) {
    sigma_sq <- 1/rgamma(1, InverseGamma_a, InverseGamma_b)
    posterior.samples <- rbind(
      posterior.samples, 
      c(mvrnorm(n = 1, mu=mvn.mean, Sigma=sigma_sq*mvn.sigma.multiplier),sqrt(sigma_sq)) )
  }
  colnames(posterior.samples) <- c(model,"sigma")
  
  ## Return
  return(posterior.samples)
}