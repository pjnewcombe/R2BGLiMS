#' Constructs a point estimate of multivariate SNP effects, based on the JAM estimator, from a vetorc of marginal
#' effects and a reference genotype matrix.
#' @export
#' @title JAM (Joint Analysis of Marginal statistics) multivariate point estimator.
#' @name JAM_PointEstimate
#' @inheritParams R2BGLiMS
#' 
#' @return A vector of multivariate point estimates.
#' 
#' The function \code{\link{PrettyResultsTable}} can be used to print summary posterior results for all parameters. Other functions
#' for summarising results are listed under "see also".
#' 
#' @seealso See also \code{\link{JAM}}.
#' 
#' @author Paul Newcombe

JAM <- function(
  marginal.betas=NULL,
  X.ref=NULL
) {
  
  # Mean centre the columns of X, to avoid inferring the intercept
  for (v in 1:ncol(X)) {
    X[,v] <- X[,v] - mean(X[,v])
  }
  n <- nrow(X) # N is arbitrary for the point estimator, use the size of the reference matrix.
  
  # Construct the z_jam outcome vector. The element for each SNP is constructed from:
  # 1) Mean y-values within each of the three genotype groups - inferred from marginal effect
  # 2) Count of each genotype - inferred from reference data
  z_jam <- NULL
  for (v in 1:length(lmBetaHat)) {
    # Genotype group counts according to proportions in reference X
    n1 <- n * table(X[,v])[2]/sum(table(X[,v]))
    n2 <- n * table(X[,v])[3]/sum(table(X[,v]))
    # Genotype group y-means (for a mean-centred y)
    y0 <- -(n1 * lmBetaHat[v] + n2 * 2 *lmBetaHat[v])/n
    y1 <- y0 + lmBetaHat[v]
    y2 <- y0 + 2*lmBetaHat[v]
    z_jam <- c(z_jam, y1 * n1 + 2 * y2 * n2) # t(X)%*%y  
  }
  
  # Take Cholesky decomposition and construct estimate. Formula is just below eq (13) in paper.
  L <- chol(t(X) %*% X)
  z_L <- solve(t(L)) %*% z_jam
  multivariate.beta.hat <- t(solve(t(L) %*% L) %*% t(L) %*% z_L)  
  
  return(multivariate.beta.hat)
}
