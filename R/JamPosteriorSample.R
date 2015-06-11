#' Conjugate results are used to generate a posterior sample of parameters (beta and sigma) under the JAM framework, 
#' conditional on a particular model, i.e. selection of SNPs. 
#' 
#' @export
#' @title Generate JAM posterior sample
#' @name JamPosteriorSample
#' @param likelihood Type of model to fit. Current options are "Logistic" (for binary data), "Weibull" (for survival data), 
#' "Gaussian" (for continuous data), "GaussianConj" (linear regression exploiting conjugate results), "GaussianMarg" (for analysis of univariate associations from Gaussian linear 
#' regressions) and "GaussianMargConj" (for analysis under a marginal conjugate linear regression model).
#' @param xTx GaussianMarg and GaussianMargConj ONLY: List containing each block's plug-in estimate for X'X.
#' @param z GaussianMarg and GaussianMargConj ONLY: Vector of quantities calculated from the summary statistics.
#' @param sigma2_invGamma_a Guassian models ONLY: Inverse-Gamma parameter one for the residual
#' precision. For the conjugate model this parameter is intergrated out, and this may be provided as in the
#' Bottolo and Richardson 2010 notation. For an informative prior for the non-conjugate model
#' (taking into account Java parameterisation) choose N/2.
#' @param sigma2_invGamma_b Guassian models ONLY: Inverse-Gamma parameter one for the residual
#' precision. For the conjugate model this parameter is intergrated out, and this may be provided as in the
#' Bottolo and Richardson 2010 notation. For an informative prior for the non-conjugate model
#' (taking into account Java parameterisation) choose N/(2*variance estimate).
#' @param g.prior Gaussian conjugate models ONLY: Whether to use a g-prior for the beta's - i.e. a multivariate normal 
#' with correlation structure proportional to sigma^2*X'X^-1 or to use independence priors (default = FALSE).
#' @param tau Gaussian conjugate models ONLY: Value to use for sparsity parameter tau (tau*sigma^2 parameterisation).
#' Default residual.var.n. If modelling this parameter, this value is used to center the Zellner-Siow prior
#' and as an intial value.
#' @param model Character vector of covariate names to include in the model.
#' @param n.samples Number of posterior samples.
#' 
#' @return The posterior sample as a matrix. Rows are different posterior samples, and columns are the covairates
#' specified in model.
#' 
#' @author Paul Newcombe

JamPosteriorSample <- function(
  xTx=NULL,
  z=NULL,
  sigma2_invGamma_a=NULL,
  sigma2_invGamma_b=NULL,
  g.prior=TRUE,
  tau=NULL,
  model=NULL,
  n.samples=1000
) {
  # Calculations
  L <- chol(xTx[[1]])
  L_g <- L[,model]
  P <- length(z)
  y <- solve(t(L))%*%z
  beta_hat <- ((solve( t(L_g) %*% L_g ))%*%t(L_g))%*%y
  sSq <- t(y-L_g%*%beta_hat) %*% (y-L_g%*%beta_hat)
  
  # Hyper-parmaters
  InverseGamma_a <- sigma2_invGamma_a + P/2
  InverseGamma_b <- sigma2_invGamma_b + sSq/2 + (t(beta_hat) %*% (t(L_g)%*%L_g) %*% beta_hat)[1,1]/(2*(tau+1))
  mvn.mean <- tau*beta_hat/(1+tau)
  mvn.sigma <- solve(t(L_g)%*%L_g)*tau/(1+tau)
  
  beta.samples <- NULL
  for (i in 1:n.samples) {
    sigma_sq <- 1/rgamma(1, InverseGamma_a, InverseGamma_b)
    beta.samples <- rbind(
      beta.samples, 
      c(mvrnorm(n = 1, mu=mvn.mean, Sigma=sigma_sq*mvn.sigma),sqrt(sigma_sq)) )
  }
  colnames(beta.samples) <- c(model,"sigma")
  
  ## Return
  return(beta.samples)
}