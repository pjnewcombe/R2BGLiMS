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
  cor.ref=NULL,
  mafs.ref=NULL,
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
    if (!is.null(X.ref)) {
      # Calculate z according to X.ref
      for (v in 1:length(marginal.betas)) {
        y.pred <- X.ref[,v]*marginal.betas[v]
        y.pred.centered <- y.pred - mean(y.pred)
        z[v] <- X.ref[,v] %*% y.pred.centered # t(X.ref)%*%y
      }
      z <- z*n/n.ref
    } else if (!is.null(cor.ref)) {
      # Calculate z according to MAFs
      for (v in 1:length(marginal.betas)) {
        n0 <- n*(1 - mafs.ref[v])^2
        n1 <- n*2*mafs.ref[v]*(1 - mafs.ref[v])
        n2 <- n*(mafs.ref[v])^2
        y0 <- -marginal.betas[v]*(n1+2*n2)/n
        y1 <- y0 + marginal.betas[v]
        y2 <- y0 + 2*marginal.betas[v]
        z[v] <- n1*y1 + 2*n2*y2
      }
    }
  } else if (!is.null(mafs.if.independent)) {
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
    if (!is.null(X.ref)) {
      # 1) Mean-centre X.ref
      for (v in 1:ncol(X.ref)) {
        X.ref[,v] <- X.ref[,v] - mean(X.ref[,v]) # MUST mean-center since z is constructed under 0 intercept 
      }
      
      # 2) Calculate MLE corresponding to the summary model
      multivariate.beta.hat <- solve(t(X.ref) %*% X.ref) %*% z*n/n.ref
      vec.return <- multivariate.beta.hat
    } else if (!is.null(X.cor)) {
      ### --- Generate X'X, from correlation matrix and MAFs
      snp.sds <- sqrt(sapply(mafs.ref, function(p) 2*p*(1-p)))
      xTx <- cor.ref
      # Scale by SDs
      for (snp1 in colnames(xTx)) {
        for (snp2 in colnames(xTx)) {
          xTx[snp1,snp2] <- xTx[snp1,snp2]*(snp.sds[snp1]*snp.sds[snp2])
        }
      }
      # Multiply by N
      xTx <- xTx*n
      
      # 2) Calculate MLE corresponding to the summary model
      multivariate.beta.hat <- solve(xTx) %*% z
      vec.return <- multivariate.beta.hat
    }
  }
  
  return(vec.return)
}
