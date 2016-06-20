#' Constructs a point estimate of multivariate SNP effects, based on the JAM estimator, from a vetorc of marginal
#' effects and a reference genotype matrix.
#' @export
#' @title JAM (Joint Analysis of Marginal statistics) multivariate point estimator.
#' @name JAM_PointEstimates
#' @inheritParams R2BGLiMS
#' 
#' @return A vector of multivariate point estimates.
#' 
#' @seealso See also \code{\link{JAM}}.
#' 
#' @author Paul Newcombe

JAM_PosteriorSampler <- function(
  marginal.betas=NULL,
  X.ref=NULL
) {
  
  # Get estimate of betahat
  
  n <- nrow(X.ref)
  
  # Construct the z_jam outcome vector. The element for each SNP is constructed from:
  # 1) Genotype counts - inferred from proportions in reference X.ref
  # 2) Mean y-values within each of the three genotype groups - inferred from marginal.betas
  z_jam <- rep(NA,length(marginal.betas))
  names(z_jam) <- names(marginal.betas)
  for (v in 1:length(marginal.betas)) {
    # Genotype group counts according to proportions in reference X.ref
    genotype.table <- table(round(X.ref[,v])) # Round in case of dosage data
    n1 <- n * genotype.table[2]/sum(genotype.table)
    n2 <- n * genotype.table[3]/sum(genotype.table)
    # Genotype group y-means (under a 0-intercept model)
    y0 <- -(n1 * marginal.betas[v] + n2 * 2 *marginal.betas[v])/n # Makes overall y mean 0
    y1 <- y0 + marginal.betas[v]
    y2 <- y0 + 2*marginal.betas[v]
    z_jam[v] <- y1 * n1 + 2 * y2 * n2 # t(X.ref)%*%y  
  }
  
  # MUST Mean centre the columns of X since z_jam is constructed under 0 intercept
  for (v in 1:ncol(X.ref)) {
    X.ref[,v] <- X.ref[,v] - mean(X.ref[,v])
  }
  
  # Take Cholesky decomposition and construct estimate. Formula is just below eq (13) in paper.
  L <- chol(t(X.ref) %*% X.ref)
  z_L <- solve(t(L)) %*% z_jam
  multivariate.beta.hat <- t(solve(t(L) %*% L) %*% t(L) %*% z_L)
  
  return(multivariate.beta.hat)
}
