#' Constructs a point estimate of multivariate SNP effects, based on the JAM estimator, from a vetorc of marginal
#' effects and a reference genotype matrix.
#' @export
#' @title JAM (Joint Analysis of Marginal statistics) multivariate point estimator.
#' @name JAM_PointEstimates
#' @inheritParams R2BGLiMS
#' @param just.get.z Return the estimated X'y outcome used in JAM (without actually calculating the corresponding point estimates). This is mainly
#' for internal use.
#' 
#' @return A vector of multivariate point estimates.
#' 
#' @seealso See also \code{\link{JAM}}.
#' 
#' @author Paul Newcombe

JAM_PointEstimates <- function(
  marginal.betas=NULL,
  X.ref=NULL,
  n=NULL,
  just.get.z=FALSE
) {
  
  # --- Setup sample sizes
  n.ref <- nrow(X.ref)
  if (is.null(n)) {
    n <- n.ref
  }
  
  ################################################
  # --- Construct the z = X'y outcome vector --- #
  ################################################
  # The element for each SNP is constructed from:
  # 1) Infer predicted y-values from the marginal.betas. Then mean-center for 0-intercept model
  # 2) Matrix multiply X.ref by the predicted y-values
  z <- rep(NA,length(marginal.betas))
  names(z) <- names(marginal.betas)
  for (v in 1:length(marginal.betas)) {
    y.pred <- X.ref[,v]*marginal.betas[v]
    y.pred.centered <- y.pred - mean(y.pred)
    z[v] <- X.ref[,v] %*% y.pred.centered # t(X.ref)%*%y
  }
  z <- z*n/n.ref
  
  #################################################
  # --- Construct multivariate beta estimates --- #
  #################################################
  
  if (just.get.z) {
    vec.return <- z
  } else {
    # 1) Get Cholesky decomposition L
    for (v in 1:ncol(X.ref)) {
      X.ref[,v] <- X.ref[,v] - mean(X.ref[,v]) # MUST mean-center since z is constructed under 0 intercept 
    }
    L <- chol(t(X.ref) %*% X.ref)
    
    # 2) Calculate MLE corresponding to the summary model, multipled through by L
    z_L <- solve(t(L)) %*% z
    multivariate.beta.hat <- t(solve(t(L) %*% L) %*% t(L) %*% z_L)
    vec.return <- multivariate.beta.hat
  }
  
  return(vec.return)
}
