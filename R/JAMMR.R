#' JAM-MR does stochastic search (RJMCMC) to select instruments robustly associated with
#' the risk factor, while penalizing pleiotropic sets
#' of SNPs according to a heterogeneity loss function (similar to the
#' het-pen approach of Burgess et al., 2018). The statistical framework for 
#' combining Bayesian inference with loss functions is provided in Bissiri, Holmes and Walker (2016).
#' The stochastic search provides posterior model probbilities and causal
#' effect estimation is then performed by averaging the individual-model
#' IVW estimates weighted by the model probabilities.
#' @export
#' @title JAMMR
#' @name JAMMR
#' @param bx Univariate effect estimates between each SNP and the risk factor.
#' @param sx Standard errors for the SNP-risk factor genetic effects.
#' @param by Univariate effect estimates between each SNP and the outcome.
#' @param sy Standard errors for the SNP-outcome genetic effects.
#' @param N1 Sample size of the study from which the bx effects were obtained.
#' @param eafs  A vector of effect allele frequencies of length same as bx. Must 
#' be provided if genetic variants are assumed to be independent.
#' @param G.matrix A reference matrix for the genetic variants used in the analysis.
#' If present, will override eafs and force JAM-MR to model genetic 
#' variants as correlated. If absent, JAM-MR will assume that genetic
#' variants are independent. One of eafs ad G.matrix must be specified.
#' @param trait.var An estimate of the variance in risk factor measurements
#' (can be obtained from the G-X GWAS and will be equal to 1
#' if the GWAS estimates are reported based on standardized data).
#' If not provided, it is internally estimated by JAM-MR.
#' @param iter The number of reversible-jump MCMC iterations to perform. 
#' @param w Tuning parameter(s) for pleiotropy penalization. If scalar, a single
#' JAM-MR implementation will be run. If vector, one implementation
#' will be run for each value and the minimum-causal-standard -error run
#' will be selected. Values should typically be multiples of N1.
#' Larger values encourage stronger penalization and smaller models.
#' @param initial.model The model to be used in the first iteration of the stochastic search.
#' @param n.models The (maximum) number of highest-posterior-probability models to be reported
#' by JAM-MR. If set to zero, no models will be reported.
#' @param jam.seed The seed to use for initializing the RJMCMC algorithm, if any.
#' If multiple JAM-MR runs are implemented, these will be seeded using
#' jam.seed, jam.seed + 1, jam.seed + 2, ...
#' 
#' @return A JAMMR results object, which is a list with the following components:
#' \enumerate{
#' \item causal: The estimated causal effect.  
#' \item se: The corresponding standard error.  
#' \item se.adj: Adjusted standard error (equal to se * 1.3178).  
#' \item model.matrix: A matrix containing the highest posterior probability models 
#' visited by JAM-MR during the variable selection procedure.
#' See the function TopModels for more details.  
#' \item snp.probs: Posterior inclusion probabilities for each SNP.  
#' \item w: The value of the tuning parameter selected.
#' \item all.causal: Causal effect estimates from all JAM-MR runs (for various w).
#' \item all.se: Causal standard errors from all JAM-MR runs (for various w).
#' \item all.w: All w values used.
#' \item all.probs: SNP inclusion probabilities for all JAM-MR runs (for various w).
#' }
#' 
#' sensitivity plot
## of causal effects and standard errors for various w.

#' @seealso See \code{\link{JAMMR_ManhattanPlot}} for a visual 
#' summary of individual SNP selection probabilities and \code{\link{JAMMR_SensitivityPlot}}
#' for producing sensitivity plots of causal effects and standard errors for various values of w.
#' 
#' @author Apostolos Gkatzionis
#' 
#' @example Examples/JAMMR_Examples.R

JAMMR <- function (
  bx,
  sx,
  by,
  sy,
  N1,
  eafs = NULL,
  G.matrix = NULL,
  trait.var = NULL,
  iter = 1e6, w,
  initial.model = rep(1, length(bx)),
  n.models = 10,
  jam.model.priors = c(1, length(bx)),
  loss.function = "variance",
  random.effects = TRUE,
  jam.seed = NULL) {
  
  ## Decide if JAM should run with independent or correlated SNPs.
  are.snps.independent <- is.null(G.matrix)
  
  ## If the trait variance is not provided, estimate it.
  if (is.null(trait.var)) {
    
    if (are.snps.independent) {
      residual.var <- sx^2 * 2 * N1 * eafs * (1 - eafs)
      trait.var <- median(2 * bx^2 * eafs * (1 - eafs) + residual.var)
    } else {
      snp.variances <- diag(t(G.matrix) %*% G.matrix) / nrow(G.matrix)
      residual.var <- sx^2 * snp.variances * N1
      trait.var <- median(bx^2 * snp.variances + residual.var)
    }
    
  }
  
  ## JAM needs to have SNP names, so generate some.
  P <-  length(bx)   ## Number of SNPs.
  snp.names <-  paste("SNP", c(1:P), sep = "")
  names(bx) <- snp.names
  if (are.snps.independent) {
    names(eafs) <- snp.names
  } else {
    colnames(G.matrix) <- snp.names
  }
  
  ## Decide whether to run a single JAM-MR implementation or many.
  nw <- length(w)
  
  if (nw == 1) {
    ## If a single run is needed, do it.
    
    paste("---------- w = ", w, " ----------")
    
    ## Now run the algorithm...
    if (are.snps.independent) {
      
      ## ... when genetic variants are independent.
      jam.results <- JAM( marginal.betas = bx, n = N1, mafs.if.independent = eafs, initial.model = initial.model,
                          model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                          n.iter = iter, seed = jam.seed, trait.variance = trait.var,
                          mrloss.w = w, mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
    } else {
      
      ## ... when genetic variants are correlated.
      jam.results <- JAM( marginal.betas = bx, n = N1, X.ref = G.matrix, initial.model = initial.model,
                          model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                          n.iter = iter, seed = jam.seed, trait.variance = trait.var,
                          mrloss.w = w, mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
    }
    
    ## Generate a table with all models visited and their posterior probabilities.
    models.visited <- as.matrix(TopModels(jam.results, n.top.models = iter, remove.empty.cols = FALSE))
    n.to.report <- min(nrow(models.visited), n.models)
    
    ## This checks if the null model was visited by JAM-MR.
    if (nrow(models.visited) == 1) {
      
      ## If the null model was the *only* model visited, there is little we can do.
      n.per.model <- sum(models.visited[1:P])
      if (n.per.model == 0) {
        theta <- NA
        sigma <- NA
        snp.probs <- rep(0, P)
        models.visited.backup <- models.visited
        warning("JAM-MR did not select any genetic variants. Genetic associations with the risk factor are probably very weak.")
      }
      
    } else {
      
      ## If the null model was one of many, discard it and use the rest for causal effect estimation.
      n.per.model <- rowSums(models.visited[, - (P + 1)])
      if ( min(n.per.model) == 0) {
        models.visited.backup <- models.visited
        zero.prob <- models.visited[which(n.per.model == 0), P + 1]
        models.visited <- models.visited[- which(n.per.model == 0), ]
        models.visited[, P + 1] <- models.visited[, P + 1] / (1 - zero.prob)
        warning("JAM-MR assigned non-zero probability to the null model.")
      }
      
    }
    
    ## This does causal effect estimation for "usual" JAMMR outputs.
    if (nrow(models.visited) > 1) {
      
      ## Compute model probabilities and inclusion probabilities per SNP.
      model.probs <- models.visited[, P + 1]
      snp.probs <- rep(0, P)
      for (ii in 1:P) snp.probs[ii] <- sum(model.probs * models.visited[, ii])
      
      ## For each model, compute IVW estimates and standard errors.
      alleffects <- rep(0, nrow(models.visited))
      allerrors <- rep(0, nrow(models.visited))
      for (i in 1:nrow(models.visited)) {
        current.model.snps <- unname(which(models.visited[i, - (P + 1)] == 1)) 
        current.thetas <- ( by / unname(bx) )[current.model.snps]
        current.sigmas <- ( sy / abs(unname(bx)) )[current.model.snps]
        
        # OLD: Relied on MendelianRandomization package
        #current.model.ivw <- mr_ivw( mr_input(bx = unname(bx[current.model.snps]), bxse = sx[current.model.snps], by = by[current.model.snps], byse = sy[current.model.snps]) )
        #alleffects[i] <- current.model.ivw$Estimate
        #allerrors[i] <- current.model.ivw$StdError
        
        # NEW: No longer depends on MendelianRandomization package
        current.ivw <- sum(current.thetas / current.sigmas^2) / sum(1 / current.sigmas^2)
        if (random.effects) current.phi <- max(1, sum( (current.thetas - current.ivw)^2 / (current.sigmas^2) ) / (length(current.model.snps) - 1) ) else current.phi <- 1
        alleffects[i] <- current.ivw
        allerrors[i] <- sqrt(current.phi) * sqrt( 1 / sum(1 / current.sigmas^2) ) 
      }
      
      ## Combine into an overall estimator and its standard error.
      theta <- sum(model.probs * alleffects)
      sigma <- sqrt( sum( model.probs * allerrors^2 ) + sum( model.probs * alleffects * (alleffects - theta) ) )
      
    } else { 
      
      ## This is the case when JAMMR assigns 100% probability to a single model. Done separately for numerical stability.
      single.model.snps <- models.visited[, - (P + 1)]
      snp.probs <- single.model.snps
      
      # OLD: Relied on MendelianRandomization package
      #single.model.ivw <- mr_ivw( mr_input(bx = unname(bx[single.model.snps]), bxse = sx[single.model.snps], by = by[single.model.snps], byse = sy[single.model.snps]) )
      
      # NEW: No longer relies on MendelianRandomization package
      single.model.thetas <- ( by / unname(bx) )[single.model.snps]
      single.model.sigmas <- ( sy / abs(unname(bx)) )[single.model.snps]
      single.model.ivw <- sum(single.model.thetas / single.model.sigmas^2) / sum(1 / single.model.sigmas^2)
      if (random.effects) single.model.phi <- max(1, sum( (single.model.thetas - single.model.ivw)^2 / (single.model.sigmas^2) ) / (length(single.model.snps) - 1) ) else single.model.phi <- 1
      theta <- single.model.ivw
      sigma <- sqrt(single.model.phi) * sqrt( 1 / sum(1 / single.model.sigmas^2) ) 
      
    }
    
    ## This adds the null model back into the top model matrix, if it is to be returned.
    if ( min(n.per.model) == 0) {
      model.matrix.to.be.returned <- models.visited.backup[1:n.to.report, ]
    } else {
      model.matrix.to.be.returned <- models.visited[1:n.to.report, ]
    }
    
    ## Return a list of results.
    return(list("causal" = theta, "se" = sigma, "se.adj" = sigma * 1.3178, "model.matrix" = model.matrix.to.be.returned,
                "snp.probs" = snp.probs, "w" = w, "all.causal" = theta, "all.se" = sigma, "all.w" = w, "all.probs" = snp.probs))
    
  } else {
    ## If multiple w values need to be run, start looping.
    
    ## Store multiple-w results here.
    all.thetas <- rep(0, nw)
    all.sigmas <- rep(0, nw)
    all.probs <- matrix(0, nw, P)
    
    ## This will be used to store the top-model matrix for the best JAM-MR run.
    current.min.stderr <- Inf
    current.top.models <- 0
    
    ## Get it going.
    for (jj in 1:nw) {
      
      paste("---------- w = ", w[jj], " ----------")
      
      ## Now run the algorithm...
      if (are.snps.independent) {
        
        ## ... when genetic variants are independent.
        jam.results <- JAM( marginal.betas = bx, n = N1, mafs.if.independent = eafs, initial.model = initial.model,
                            model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                            n.iter = iter, seed = jam.seed + jj - 1, trait.variance = trait.var,
                            mrloss.w = w[jj], mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
      } else {
        
        ## ... when genetic variants are correlated.
        jam.results <- JAM( marginal.betas = bx, n = N1, X.ref = G.matrix, initial.model = initial.model,
                            model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                            n.iter = iter, seed = jam.seed + jj - 1, trait.variance = trait.var,
                            mrloss.w = w[jj], mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
      }
      
      ## Generate a table with all models visited and their posterior probabilities.
      models.visited <- as.matrix(TopModels(jam.results, n.top.models = iter, remove.empty.cols = FALSE))
      n.to.report <- min(nrow(models.visited), n.models)
      
      ## This checks if the null model was visited by JAM-MR.
      if (nrow(models.visited) == 1) {
        
        ## If the null model was the *only* model visited, there is little we can do.
        n.per.model <- sum(models.visited[1:P])
        if (n.per.model == 0) {
          theta <- NA
          sigma <- NA
          snp.probs <- rep(0, P)
          models.visited.backup <- models.visited
          warning("JAM-MR did not select any genetic variants. Genetic associations with the risk factor are probably very weak.")
        }
        
      } else {
        
        ## If the null model was one of many, discard it and use the rest for causal effect estimation.
        n.per.model <- rowSums(models.visited[, - (P + 1)])
        if ( min(n.per.model) == 0) {
          models.visited.backup <- models.visited
          zero.prob <- models.visited[which(n.per.model == 0), P + 1]
          models.visited <- models.visited[- which(n.per.model == 0), ]
          models.visited[, P + 1] <- models.visited[, P + 1] / (1 - zero.prob)
          warning("JAM-MR assigned non-zero probability to the null model.")
        }
        
      }
      
      ## This does causal effect estimation for "usual" JAMMR outputs.
      if (nrow(models.visited) > 1) {
        
        ## Compute model probabilities and inclusion probabilities per SNP.
        model.probs <- models.visited[, P + 1]
        snp.probs <- rep(0, P)
        for (ii in 1:P) snp.probs[ii] <- sum(model.probs * models.visited[, ii])
        
        ## For each model, compute IVW estimates and standard errors.
        alleffects <- rep(0, nrow(models.visited))
        allerrors <- rep(0, nrow(models.visited))
        for (i in 1:nrow(models.visited)) {
          
          current.model.snps <- unname(which(models.visited[i, - (P + 1)] == 1)) 
          current.thetas <- ( by / unname(bx) )[current.model.snps]
          current.sigmas <- ( sy / abs(unname(bx)) )[current.model.snps]
          
          # OLD: Relied on MendelianRandomization package
          #current.model.ivw <- mr_ivw( mr_input(bx = unname(bx[current.model.snps]), bxse = sx[current.model.snps], by = by[current.model.snps], byse = sy[current.model.snps]) )
          #alleffects[i] <- current.model.ivw$Estimate
          #allerrors[i] <- current.model.ivw$StdError
          
          # NEW: No longer relies on MendelianRandomization package
          current.ivw <- sum(current.thetas / current.sigmas^2) / sum(1 / current.sigmas^2)
          if (random.effects) current.phi <- max(1, sum( (current.thetas - current.ivw)^2 / (current.sigmas^2) ) / (length(current.model.snps) - 1) ) else current.phi <- 1
          alleffects[i] <- current.ivw
          allerrors[i] <- sqrt(current.phi) * sqrt( 1 / sum(1 / current.sigmas^2) ) 
        }
        
        ## Combine into an overall estimator and its standard error.
        all.thetas[jj] <- sum(model.probs * alleffects)
        all.sigmas[jj] <- sqrt( sum( model.probs * allerrors^2 ) + sum( model.probs * alleffects * (alleffects - all.thetas[jj]) ) )
        all.probs[jj, ] <- snp.probs
        
      } else { 
        
        ## This is the case when JAMMR assigns 100% probability to a single model. Done separately for numerical stability.
        single.model.snps <- models.visited[, - (P + 1)]
        all.probs[jj, ] <- single.model.snps
        
        # OLD: Relied on MendelianRandomization package
        #single.model.ivw <- mr_ivw( mr_input(bx = unname(bx[single.model.snps == 1]), bxse = sx[single.model.snps == 1], by = by[single.model.snps == 1], byse = sy[single.model.snps == 1]) )
        #all.probs[jj, ] <- single.model.snps
        #all.thetas[jj] <- single.model.ivw$Estimate
        #all.sigmas[jj] <- single.model.ivw$StdError
        
        # NEW: No longer relies on MendelianRandomization package
        single.model.thetas <- ( by / unname(bx) )[single.model.snps]
        single.model.sigmas <- ( sy / abs(unname(bx)) )[single.model.snps]
        single.model.ivw <- sum(single.model.thetas / single.model.sigmas^2) / sum(1 / single.model.sigmas^2)
        if (random.effects) single.model.phi <- max(1, sum( (single.model.thetas - single.model.ivw)^2 / (single.model.sigmas^2) ) / (length(single.model.snps) - 1) ) else single.model.phi <- 1
        all.thetas[jj] <- single.model.ivw
        all.sigmas[jj] <- sqrt(single.model.phi) * sqrt( 1 / sum(1 / single.model.sigmas^2) ) 
        
      }
      
      ## If we have obtained a smaller standard error, store the top models matrix.
      if (current.min.stderr > all.sigmas[jj]) {
        current.min.stderr <- all.sigmas[jj]
        if ( min(n.per.model) == 0) {
          current.top.models <- models.visited.backup[1:n.to.report, ]
        } else {
          current.top.models <- models.visited[1:n.to.report, ]
        }
      }
      
    }
    
    ## Select the best w value.
    best.w <- w[which.min(all.sigmas)]
    theta <- all.thetas[which.min(all.sigmas)]
    sigma <- min(all.sigmas, na.rm = TRUE)
    snp.probs <- all.probs[which.min(all.sigmas), ]
    
    ## Return a list of results.
    return(list("causal" = theta, "se" = sigma, "se.adj" = sigma * 1.3178, "model.matrix" = current.top.models,
                "snp.probs" = snp.probs, "w" = best.w, "all.causal" = all.thetas, "all.se" = all.sigmas, 
                "all.w" = w, "all.probs" = all.probs))
    
  }
  
}
