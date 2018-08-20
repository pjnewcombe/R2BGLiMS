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
  just.get.z=FALSE,
  mafs.if.independent=NULL
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
  if (is.null(mafs.if.independent)) {
    for (v in 1:length(marginal.betas)) {
      y.pred <- X.ref[,v]*marginal.betas[v]
      y.pred.centered <- y.pred - mean(y.pred)
      z[v] <- X.ref[,v] %*% y.pred.centered # t(X.ref)%*%y
    }
    z <- z*n/n.ref
  } else {
    for (v in 1:length(marginal.betas)) {
      n0 <- n*(1 - mafs.if.independent[v])^2
      n1 <- n*2*mafs.if.independent[v]*(1 - mafs.if.independent[v])
      n2 <- n*(mafs.if.independent[v])^2
      y0 <- -marginal.betas[v]*(n1+2*n2)/n
      y1 <- y0 + marginal.betas[v]
      y2 <- y0 + 2*marginal.betas[v]
      z[v] <- n1*y1 + 2*n2*y2
    }
  }
  
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
