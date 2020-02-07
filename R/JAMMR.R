#' This function implements the JAM-MR algorithm for Mendelian randomization.
#' Given a set of genetic variants and corresponding summary statistics for
#' associations with a risk factor and a disease outcome of interest, JAM-MR
#' performs Bayesian stochastic search through a reversible-jump MCMC algorithm 
#' to identify the most suitable variants for inclusion in a subsequent Mendelian
#' randomization analysis. The algorithm prioritizes genetic variants with strong 
#' associations with the risk factor while downweighting variants with heterogeneous 
#' ratio estimates through the use of the general Bayesian framework and a heterogeneity-
#' penalizing loss function. The stochastic search returns a list of models with 
#' associated posterior probabilties. For each of them, JAM-MR obtains a model-specific
#' causal effect estimate by fitting a normal or (if mretn = TRUE) a truncated normal
#' multiplictive random-effects model. Finally, the algorithm averages across the 
#' model-specific estimates to return an aggregate causal effect estimate and its 
#' standard error.
#' @export
#' @title JAMMR - Bayesian variable selection for Mendelian randomization.
#' @name JAMMR
#' @param bx Univariate effect estimates between each genetic variant and the risk factor.
#' @param sx Standard errors for the varant-risk factor genetic effects.
#' @param by Univariate effect estimates between each genetic variant and the outcome.
#' @param sy Standard errors for the variant-outcome genetic effects.
#' @param N1 Sample size of the study from which the bx effects were obtained.
#' @param eafs A vector of effect allele frequencies of length same as bx.
#' Must be provided if genetic variants are assumed to be independent.
#' @param G.matrix A reference matrix for the genetic variants used in the analysis.
#' If present, will override eafs and force JAM-MR to model genetic variants as 
#' correlated. If absent, JAM-MR will assume that genetic variants are independent. 
#' One of eafs ad G.matrix must be specified. Note: the algorithm takes longer to
#' run with coorrelated instruments, therefore when genetic variants are truly  
#' independent it is better to specify eafs and not G.matrix.
#' @param G.cor The genetic correlation matrix. If an analysis with correlated SNPs 
#' is implemented but the matrix of raw genetic data (G.matrix) is unavailable or is too 
#' big to store in memory, JAM-MR can be run by providing both eafs and G.cor instead.
#' @param trait.var An estimate of the variance in risk factor measurements. Can be 
#' obtained from the G-X GWAS and will be equal to 1 if the GWAS estimates are reported
#' based on standardized data. If not provided, it is internally estimated by JAM-MR.
#' @param iter The number of reversible-jump MCMC iterations to perform. 
#' @param w Tuning parameter(s) for pleiotropy penalization. If scalar, a single
#' JAM-MR implementation will be run. If vector, one implementation will be run 
#' for each value and the minimum-causal-standard-error run will be selected. 
#' Values should typically be multiples of N1. Larger values encourage stronger 
#' penalization and smaller models.
#' @param n.grid If w is not specified, it is estimated by a grid search; n.grid is 
#' the number of grid-points to visit during the grid search. Defaults to 26. 
#' Redundant if w is specified.
#' @param grid.limits Lower and upper limits for the w values to consider during the 
#' grid search. The algorithm will create a grid of size n.grid by logarithmic
#' interpolation between the limits, plus the value w = 0. The default is a grid  
#' search between w = 0.01 N1 and w = 100 N1. Redundant if w is specified.
#' @param initial.model The model to be used in the first iteration of the stochastic 
#' search. By default the algorithm starts from a model with all genetic variants 
#' included, when running on independent variants, and a model with only the smallest
#' p-value genetic variant included, when running on correlated variants.
#' @param n.models The (maximum) number of highest posterior probability models to be 
#' reported by JAM-MR. If set to zero, no models will be reported.
#' @param mretn Logical. If TRUE, the algorithm will fit a multiplicative random-effects
#' truncated normal model in order to compute model-specific standard errors when 
#' averaging across models. If FALSE, random-effects IVW estimates will be used instead.  
#' Currently this is only implemented for independent genetic variants; with correlated 
#' variants the algorithm will always use IVW as model-specific estimates.
#' @param jam.model.priors Prior parameters for JAM's beta-binomial model space prior.
#' @param loss.function The loss function to use. Weighted if "steve", unweighted if "variance".
#' @param jam.seed The seed to use for initializing the RJMCMC algorithm, if any. If 
#' a grid search for w is implemented, the implementations for different grid points  
#' will be seeded using the values jam.seed, jam.seed + 1, jam.seed + 2, ...
#' 
#' 
#' @return A list with the following arguments:
#' \itemize{
#'   \item causal - Estimated causal effect.
#'   \item se - Causal standard error.
#'   \item model.matrix - If n.models > 0, a matrix containing the highest posterior
#'   probability models visited by JAM-MR during the variable selection procedure. See
#'   the function TopModels for more details.
#'   \item snp.probs - Posterior inclusion probabilities for each genetic variant.
#'   \item w - The value of the tuning parameter selected.
#' }
#' In addition, if a grid search for w is implemented:
#' \itemize{
#'   \item all.causal - Causal effect estimates for each w value.
#'   \item all.se - Causal standard errors for each w value.
#'   \item all.probs - Posterior inclusion probabilities per genetic variant for each w value.
#'   \item all.w - All w values used as part of the grid search.
#' }
#'   
#' @seealso See \code{\link{JAMMR_ManhattanPlot}} for a visual summary of covariate
#' selection probabilities and \code{\link{JAMMR_SensitivityPlot}} for a visual 
#' comparison of results obtained using different w values. Also see \code{\link{JAM}} 
#' for an implementation of the fine-mapping algorithm on which JAM-MR relies.
#' 
#' @author Apostolos Gkatzionis and Paul Newcombe
#'
#' @example Examples/JAMMR_Examples.R

JAMMR <- function (bx, sx, by, sy, N1, eafs = NULL, G.matrix = NULL, G.cor = NULL, trait.var = NULL, iter = 1e6, 
                   w = NULL, n.grid = 26, grid.limits = NULL, initial.model = NULL, n.models = 10, 
                   mretn = TRUE, jam.model.priors = c(1, length(bx)), loss.function = "steve", jam.seed = NULL) {
  
  ## ---------- SETUP ---------- ##
  
  ## Set up the grid of w values, if not user-specified.
  if (is.null(w)) {
    
    if (n.grid == 1) {
      w <- N1
    } else {
      if (is.null(grid.limits)) {
        if (n.grid < 10) {
          grid.limits <- c(0.1 * N1, 10 * N1)
        } else {
          grid.limits <- c(0.01 * N1, 100 * N1)
        }
      }
      w <- c(0, exp(seq(from = log(grid.limits[1]), to = log(grid.limits[2]), length = n.grid - 1)))
    }
    
  }
  
  ## Decide whether to run a single JAM-MR implementation or many.
  nw <- length(w)
  
  ## Compute univariate causal effect estimates and standard errors.
  theta.univ <- by / bx
  sigma.univ <- sy / abs(bx)
  
  ## Decide if JAM should run with independent or correlated SNPs.
  are.snps.independent <- is.null(G.matrix) & is.null(G.cor)
  
  
  
  ## The cases of independent and correlated SNPs are treated separately.
  if (are.snps.independent) {
    
    
    ## -------------------------------------------------- ##
    ## ---------- JAM-MR WITH INDEPENDENT SNPS ---------- ##
    ## -------------------------------------------------- ##
    
    
    ## If the trait variance is not provided, estimate it.
    if (is.null(trait.var)) {
      residual.var <- sx^2 * 2 * N1 * eafs * (1 - eafs)
      trait.var <- median(2 * bx^2 * eafs * (1 - eafs) + residual.var)
    }
    
    ## If the initial model is not provided, use the full model.
    if (is.null(initial.model)) initial.model <- rep(1, length(bx))
    
    ## JAM needs to have SNP names, so generate some.
    P <-  length(bx)   ## Number of SNPs.
    snp.names <-  paste("SNP", c(1:P), sep = "")
    names(bx) <- snp.names
    names(eafs) <- snp.names
    
    ## If nw = 1 run the algorithm with no grid search (a single w) otherwise run it multiple times.
    if (nw == 1) {
      
      
      ## ---------- NW = 1 ---------- ##
      
      
      ## If a single run is needed, do it.
      print(paste("---------- w = ", w, " ----------"))
      
      ## Run JAM.
      jam.results <- JAM( marginal.betas = bx, n = N1, mafs.if.independent = eafs, initial.model = initial.model,
                          model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                          n.iter = iter, seed = jam.seed, trait.variance = trait.var, mrloss.w = w, 
                          mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
      
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
          if (zero.prob > 0.2) warning("Genetic associations with the risk factor are probably fairly weak.")
          models.visited <- models.visited[- which(n.per.model == 0), ]
          models.visited[, P + 1] <- models.visited[, P + 1] / (1 - zero.prob)
        }
        
      }
      
      
      ## For numerical stability we treat separately the case where JAM-MR only visits one model.
      if (nrow(models.visited) > 1) {
        
        ## ----- JAM-MR visits multiple models ----- ##
        
        ## Compute model probabilities and inclusion probabilities per SNP.
        model.probs <- models.visited[, P + 1]
        snp.probs <- rep(0, P)
        for (ii in 1:P) snp.probs[ii] <- sum(model.probs * models.visited[, ii])
        
        ## For each model, compute IVW estimates and standard errors.
        alleffects <- rep(0, nrow(models.visited))
        allerrors <- rep(0, nrow(models.visited))
        
        for (i in 1:nrow(models.visited)) {
          
          ## Extract the current model's SNPs, thetas and sigmas.
          current.model.snps <- unname(which(models.visited[i, - (P + 1)] == 1)) 
          current.thetas <- theta.univ[current.model.snps]
          current.sigmas <- sigma.univ[current.model.snps]
          
          ## If the multiplicative random-effects truncated normal model is needed, fit it.
          if (mretn == TRUE) {
            
            ## Compute truncation points.
            left.cutoff <- max(-Inf, theta.univ[theta.univ < min(current.thetas)])
            right.cutoff <- min(Inf, theta.univ[theta.univ > max(current.thetas)])
            
            ## Fit the model by maximum likelihood estimation.
            suppressWarnings({
              try({
                mt.fit <- constrOptim(theta = c(mean(current.thetas), 1.1), f = .mt.lik, grad = .mt.grad, ui = matrix(c(0, 1), 1, 2), ci = 1, 
                                      method = "BFGS", l = left.cutoff, u = right.cutoff, thetas = current.thetas, sigmas = current.sigmas)
                mt.fit.hessian <- .mt.hess(mt.fit$par, l = left.cutoff, u = right.cutoff, thetas = current.thetas, sigmas = current.sigmas)
                alleffects[i] <- mt.fit$par[1]
                allerrors[i] <- sqrt( solve(mt.fit.hessian)[1, 1] )
              }, silent = TRUE)
            })
            
          } else {
            
            ## Otherwise, use standard Mendelian randomization.
            current.ivw <- sum(current.thetas / current.sigmas^2) / sum(1 / current.sigmas^2)
            current.phi <- max(1, sum( (current.thetas - current.ivw)^2 / (current.sigmas^2) ) / (length(current.model.snps) - 1) )
            alleffects[i] <- current.ivw
            allerrors[i] <- sqrt(current.phi) * sqrt( 1 / sum(1 / current.sigmas^2) ) 
            
          }
          
        }
        
        ## If some (but not all) models returned "NaN", make the necesary adjustments.
        if (sum(model.probs[allerrors != "NaN"]) > 0.1) {
          model.probs <- model.probs[allerrors != "NaN"] / sum(model.probs[allerrors != "NaN"])
          alleffects <- alleffects[allerrors != "NaN"]
          allerrors <- allerrors[allerrors != "NaN"]
        }
        
        ## Combine into an overall estimator and its standard error.
        theta <- sum(model.probs * alleffects)
        sigma <- sqrt( sum( model.probs * allerrors^2 ) + sum( model.probs * alleffects * (alleffects - theta) ) )
        
      } else { 
        
        ## ----- JAM-MR visits only a single model ----- ##
        
        ## Extract the current model's SNPs, thetas, sigmas.
        single.model.snps <- unname(which(models.visited[, - (P + 1)] == 1))
        snp.probs <- models.visited[, - (P + 1)]
        single.model.thetas <- theta.univ[single.model.snps]
        single.model.sigmas <- sigma.univ[single.model.snps]
        
        ## If the multiplicative random-effects truncated normal model is needed, fit it.
        if (mretn == TRUE) {
          
          ## Compute truncation points.
          left.cutoff <- max(-Inf, theta.univ[theta.univ < min(single.model.thetas)])
          right.cutoff <- min(Inf, theta.univ[theta.univ > max(single.model.thetas)])
          
          ## Fit the model by maximum likelihood estimation.
          suppressWarnings({
            try({
              mt.fit <- constrOptim(theta = c(mean(single.model.thetas), 1.1), f = .mt.lik, grad = .mt.grad, ui = matrix(c(0, 1), 1, 2), ci = 1, 
                                    method = "BFGS", l = left.cutoff, u = right.cutoff, thetas = single.model.thetas, sigmas = single.model.sigmas)
              mt.fit.hessian <- .mt.hess(mt.fit$par, l = left.cutoff, u = right.cutoff, thetas = single.model.thetas, sigmas = single.model.sigmas)
              theta <- mt.fit$par[1]
              sigma <- sqrt( solve(mt.fit.hessian)[1, 1] )
            }, silent = TRUE)
          })
          
          ## For a single model and a single w there is no adjustment we can do, so use standard MR if it fails.
          if(sigma == "NaN") {
            
            single.model.ivw <- sum(single.model.thetas / single.model.sigmas^2) / sum(1 / single.model.sigmas^2)
            single.model.phi <- max(1, sum( (single.model.thetas - single.model.ivw)^2 / (single.model.sigmas^2) ) / (length(single.model.snps) - 1) )
            theta <- single.model.ivw
            sigma <- sqrt(single.model.phi) * sqrt( 1 / sum(1 / single.model.sigmas^2) ) 
            warning("The truncated normal fit was unstable -  standard IVW has been used instead.")
            
          } 
          
        } else {
          
          ## Otherwise, use standard Mendelian randomization.
          single.model.ivw <- sum(single.model.thetas / single.model.sigmas^2) / sum(1 / single.model.sigmas^2)
          single.model.phi <- max(1, sum( (single.model.thetas - single.model.ivw)^2 / (single.model.sigmas^2) ) / (length(single.model.snps) - 1) )
          theta <- single.model.ivw
          sigma <- sqrt(single.model.phi) * sqrt( 1 / sum(1 / single.model.sigmas^2) ) 
          
        }
        
      }
      
      ## This adds the null model back into the top model matrix, if it is to be returned.
      if ( min(n.per.model) == 0) {
        model.matrix.to.be.returned <- models.visited.backup[1:n.to.report, ]
      } else {
        model.matrix.to.be.returned <- models.visited[1:n.to.report, ]
      }
      
      ## Return a list of results. This ends the "nw = 1" case.
      return(list("causal" = theta, "se" = sigma, "model.matrix" = model.matrix.to.be.returned,
                  "snp.probs" = snp.probs, "w" = w))
      
      
    } else {
      
      
      ## ---------- NW > 1 ---------- ##
      
      
      ## If multiple w values need to be run, loop over them.
      
      ## Store per-w results here.
      all.thetas <- rep(0, nw)
      all.sigmas <- rep(0, nw)
      all.probs <- matrix(0, nw, P)
      
      ## This will be used to store the top-model matrix for the best JAM-MR run.
      current.min.stderr <- Inf
      current.top.models <- 0
      
      ## This stores warning signs when JAM visits the null model.
      warning.signs <- rep(0, nw)
      
      ## Get the loop going.
      for (jj in 1:nw) {
        
        print(paste("---------- w = ", w[jj], " ----------"))
        
        ## Run JAM.
        jam.results <- JAM( marginal.betas = bx, n = N1, mafs.if.independent = eafs, initial.model = initial.model,
                            model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                            n.iter = iter, seed = jam.seed + jj - 1, trait.variance = trait.var, mrloss.w = w[jj], 
                            mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
        
        
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
            warning.signs[jj] <- 2
          }
          
        } else {
          
          ## If the null model was one of many, discard it and use the rest for causal effect estimation.
          n.per.model <- rowSums(models.visited[, - (P + 1)])
          if ( min(n.per.model) == 0) {
            models.visited.backup <- models.visited
            zero.prob <- models.visited[which(n.per.model == 0), P + 1]
            models.visited <- models.visited[- which(n.per.model == 0), ]
            models.visited[, P + 1] <- models.visited[, P + 1] / (1 - zero.prob)
            if (zero.prob > 0.2) warning.signs[jj] <- 1
          }
          
        }
        
        
        ## For numerical stability we treat separately the case where JAM-MR only visits one model.
        if (nrow(models.visited) > 1) {
          
          ## ----- JAM-MR visits multiple models ----- ##
          
          ## Compute model probabilities and inclusion probabilities per SNP.
          model.probs <- models.visited[, P + 1]
          snp.probs <- rep(0, P)
          for (ii in 1:P) snp.probs[ii] <- sum(model.probs * models.visited[, ii])
          
          ## For each model, compute IVW estimates and standard errors.
          alleffects <- rep(0, nrow(models.visited))
          allerrors <- rep(0, nrow(models.visited))
          
          for (i in 1:nrow(models.visited)) {
            
            ## Extract the current model's SNPs, thetas and sigmas.
            current.model.snps <- unname(which(models.visited[i, - (P + 1)] == 1)) 
            current.thetas <- theta.univ[current.model.snps]
            current.sigmas <- sigma.univ[current.model.snps]
            
            ## If the multiplicative random-effects truncated normal model is needed, fit it.
            if (mretn == TRUE) {
              
              ## Compute truncation points.
              left.cutoff <- max(-Inf, theta.univ[theta.univ < min(current.thetas)])
              right.cutoff <- min(Inf, theta.univ[theta.univ > max(current.thetas)])
              
              ## Fit the model by maximum likelihood estimation.
              suppressWarnings({
                try({
                  mt.fit <- constrOptim(theta = c(mean(current.thetas), 1.1), f = .mt.lik, grad = .mt.grad, ui = matrix(c(0, 1), 1, 2), ci = 1, 
                                        method = "BFGS", l = left.cutoff, u = right.cutoff, thetas = current.thetas, sigmas = current.sigmas)
                  mt.fit.hessian <- .mt.hess(mt.fit$par, l = left.cutoff, u = right.cutoff, thetas = current.thetas, sigmas = current.sigmas)
                  alleffects[i] <- mt.fit$par[1]
                  allerrors[i] <- sqrt( solve(mt.fit.hessian)[1, 1] )
                }, silent = TRUE)
              })
              
            } else {
              
              ## Otherwise, use standard Mendelian randomization.
              current.ivw <- sum(current.thetas / current.sigmas^2) / sum(1 / current.sigmas^2)
              current.phi <- max(1, sum( (current.thetas - current.ivw)^2 / (current.sigmas^2) ) / (length(current.model.snps) - 1) )
              alleffects[i] <- current.ivw
              allerrors[i] <- sqrt(current.phi) * sqrt( 1 / sum(1 / current.sigmas^2) ) 
              
            }
            
          }
          
          ## If some (but not many) models returned "NaN", make the necesary adjustments.
          if (sum(model.probs[allerrors != "NaN"]) > 0.5) {
            model.probs <- model.probs[allerrors != "NaN"] / sum(model.probs[allerrors != "NaN"])
            alleffects <- alleffects[allerrors != "NaN"]
            allerrors <- allerrors[allerrors != "NaN"]
          }
          
          ## Combine into an overall estimator and its standard error.
          all.thetas[jj] <- sum(model.probs * alleffects)
          all.sigmas[jj] <- sqrt( sum( model.probs * allerrors^2 ) + sum( model.probs * alleffects * (alleffects - all.thetas[jj]) ) )
          all.probs[jj, ] <- snp.probs
          
          
        } else { 
          
          ## ----- JAM-MR visits only a single model ----- ##
          
          ## Extract the current model's SNPs, thetas and sigmas.
          single.model.snps <- unname(which(models.visited[, - (P + 1)] == 1)) 
          all.probs[jj, ] <- models.visited[, - (P + 1)]
          single.model.thetas <- theta.univ[single.model.snps]
          single.model.sigmas <- sigma.univ[single.model.snps]
          
          ## If the multiplicative random-effects truncated normal model is needed, fit it.
          if (mretn == TRUE) {
            
            ## Compute truncation points.
            left.cutoff <- max(-Inf, theta.univ[theta.univ < min(single.model.thetas)])
            right.cutoff <- min(Inf, theta.univ[theta.univ > max(single.model.thetas)])
            
            ## Fit the model by maximum likelihood estimation.
            suppressWarnings({
              try({
                mt.fit <- constrOptim(theta = c(mean(single.model.thetas), 1.1), f = .mt.lik, grad = .mt.grad, ui = matrix(c(0, 1), 1, 2), ci = 1, 
                                      method = "BFGS", l = left.cutoff, u = right.cutoff, thetas = single.model.thetas, sigmas = single.model.sigmas)
                mt.fit.hessian <- .mt.hess(mt.fit$par, l = left.cutoff, u = right.cutoff, thetas = single.model.thetas, sigmas = single.model.sigmas)
                all.thetas[jj] <- mt.fit$par[1]
                all.sigmas[jj] <- sqrt( solve(mt.fit.hessian)[1, 1] )
              }, silent = TRUE)
            })
            
          } else {
            
            ## Otherwise, use standard Mendelian randomization.
            single.model.ivw <- sum(single.model.thetas / single.model.sigmas^2) / sum(1 / single.model.sigmas^2)
            single.model.phi <- max(1, sum( (single.model.thetas - single.model.ivw)^2 / (single.model.sigmas^2) ) / (length(single.model.snps) - 1) )
            all.thetas[jj] <- single.model.ivw
            all.sigmas[jj] <- sqrt(single.model.phi) * sqrt( 1 / sum(1 / single.model.sigmas^2) ) 
            
          }
          
        }
        
        ## ----- Round up the per-w iteration ----- ##
        
        ## If we have obtained a smaller standard error, store the top models matrix.
        if (!(all.sigmas[jj] %in% c(NA, NaN))) {
          if (current.min.stderr > all.sigmas[jj]) {
            current.min.stderr <- all.sigmas[jj]
            if ( min(n.per.model) == 0) {
              current.top.models <- models.visited.backup[1:n.to.report, ]
            } else {
              current.top.models <- models.visited[1:n.to.report, ]
            }
          }
        }
        
        
      }
      
      ## ----- Round up the nw > 1 case ----- ##
      
      ## There is still a small chance that things will go wrong.
      if (all(is.na(all.sigmas))) {
        
        ## If things go wrong, warn.
        warning("JAM-MR failed to obtain a reasonable causal standard error. This may be due to the specified range of w values being too narrow.")
        best.w <- w[1]
        theta <- all.thetas[1]
        sigma <- NA
        snp.probs <- all.probs[1, ]
        
      } else {
        
        ## Otherwise, select the best w value and estimates.
        best.w <- w[which.min(all.sigmas)]
        theta <- all.thetas[which.min(all.sigmas)]
        sigma <- min(all.sigmas, na.rm = TRUE)
        snp.probs <- all.probs[which.min(all.sigmas), ]
        
        ## If the best w corresponds to a JAM-MR run that visited the null model, warn.
        if (warning.signs[which.min(all.sigmas)] == 2) {
          warning("JAM-MR did not select any genetic variants. Genetic associations with the risk factor are probably very weak.")
        }
        if (warning.signs[which.min(all.sigmas)] == 1) {
          warning("Genetic associations with the risk factor are probably fairly weak.")
        }

      }
      
      ## Return a list of results. This ends the nw > 1 case.
      return(list("causal" = theta, "se" = sigma, "model.matrix" = current.top.models, "snp.probs" = snp.probs, 
                  "w" = best.w, "all.causal" = all.thetas, "all.se" = all.sigmas, "all.w" = w, "all.probs" = all.probs))
      
    }
    
    
  } else {
    
    
    ## ------------------------------------------------- ##
    ## ---------- JAM-MR WITH CORRELATED SNPS ---------- ##
    ## ------------------------------------------------- ##
    
    
    ## Is a full genetic matrix provided?
    is.G.provided <- !(is.null(G.matrix))
    
    ## If the trait variance is not provided, estimate it.
    if (is.null(trait.var)) {
      
      if (is.G.provided) {
        
        G.centered <- t(t(G.matrix) - colMeans(G.matrix))
        snp.variances <- diag(t(G.centered) %*% G.centered) / nrow(G.centered)
        residual.var <- sx^2 * snp.variances * N1
        trait.var <- median(bx^2 * snp.variances + residual.var)
        rm(G.centered)   ## To save memory.
        
      } else {
        
        residual.var <- sx^2 * 2 * N1 * eafs * (1 - eafs)
        trait.var <- median(2 * bx^2 * eafs * (1 - eafs) + residual.var)
        
      }
    }
    
    ## If the initial model is not provided, use a 1-SNP model.
    if (is.null(initial.model)) {
      initial.model <- rep(0, length(bx))
      initial.model[which.max(abs(bx / sx))] <- 1
    }
      
    ## JAM needs to have SNP names, so generate some.
    P <-  length(bx)   ## Number of SNPs.
    snp.names <-  paste("SNP", c(1:P), sep = "")
    names(bx) <- snp.names
    if (is.G.provided) {
      colnames(G.matrix) <- snp.names
    } else {
      colnames(G.cor) <- snp.names
      rownames(G.cor) <- snp.names
      names(eafs) <- snp.names
    }
    
    ## Compute genetic correlations from the reference matrix.
    if (is.G.provided) Rho <- cor(G.matrix) else Rho <- G.cor
    
    
    ## If nw = 1 run the algorithm with no grid search (a single w) otherwise run it multiple times.
    if (nw == 1) {
      
      
      ## ---------- NW = 1 ---------- ##
      
      
      ## If a single run is needed, do it.
      print(paste("---------- w = ", w, " ----------"))
      
      ## Run JAM.
      if (is.G.provided) {
        jam.results <- JAM( marginal.betas = bx, n = N1, X.ref = G.matrix, initial.model = initial.model,
                            model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                            n.iter = iter, seed = jam.seed, trait.variance = trait.var, mrloss.w = w, 
                            mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
      } else {
        jam.results <- JAM( marginal.betas = bx, n = N1, cor.ref = G.cor, mafs.ref = eafs, initial.model = initial.model,
                            model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                            n.iter = iter, seed = jam.seed, trait.variance = trait.var, mrloss.w = w, 
                            mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
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
          if (zero.prob > 0.2) warning("Genetic associations with the risk factor are probably fairly weak.")
        }
        
      }
      
      
      ## For numerical stability we treat separately the case where JAM-MR only visits one model.
      if (nrow(models.visited) > 1) {
        
        ## ----- JAM-MR visits multiple models ----- ##
        
        ## Compute model probabilities and inclusion probabilities per SNP.
        model.probs <- models.visited[, P + 1]
        snp.probs <- rep(0, P)
        for (ii in 1:P) snp.probs[ii] <- sum(model.probs * models.visited[, ii])
        
        ## For each model, compute IVW estimates and standard errors.
        alleffects <- rep(0, nrow(models.visited))
        allerrors <- rep(0, nrow(models.visited))
        
        for (i in 1:nrow(models.visited)) {
          
          ## Extract the current model's SNPs, thetas and sigmas.
          current.model.snps <- unname(which(models.visited[i, - (P + 1)] == 1)) 
          current.bx <- unname(bx[current.model.snps])
          current.by <- by[current.model.snps]
          current.sy <- sy[current.model.snps]
          
          ## Truncation and random-effects are not currently supported for correlated data, so do IVW instead.
          current.invOmega <- solve( (current.sy %*% t(current.sy)) * Rho[current.model.snps, current.model.snps] )
          allerrors[i] <- sqrt( as.numeric( solve( t(current.bx) %*% current.invOmega %*% current.bx ) ) )
          alleffects[i] <- as.numeric( t(current.bx) %*% current.invOmega %*% current.by ) * allerrors[i]^2
          
        }
        
        ## Combine into an overall estimator and its standard error.
        theta <- sum(model.probs * alleffects)
        sigma <- sqrt( sum( model.probs * allerrors^2 ) + sum( model.probs * alleffects * (alleffects - theta) ) )
        
      } else { 
        
        ## ----- JAM-MR visits only a single model ----- ##
        
        ## Extract the current model's SNPs, thetas, sigmas.
        single.model.snps <- unname(which(models.visited[, - (P + 1)] == 1))
        snp.probs <- models.visited[, - (P + 1)]
        current.bx <- unname(bx[single.model.snps])
        current.by <- by[single.model.snps]
        current.sy <- sy[single.model.snps]
        
        ## Truncation and random-effects are not currently supported for correlated data, so do IVW instead.
        current.invOmega <- solve( (current.sy %*% t(current.sy)) * Rho[single.model.snps, single.model.snps] )
        sigma <- sqrt( as.numeric( solve( t(current.bx) %*% current.invOmega %*% current.bx ) ) )
        theta <- as.numeric( t(current.bx) %*% current.invOmega %*% current.by ) * sigma^2
        
      }
      
      ## This adds the null model back into the top model matrix, if it is to be returned.
      if ( min(n.per.model) == 0) {
        model.matrix.to.be.returned <- models.visited.backup[1:n.to.report, ]
      } else {
        model.matrix.to.be.returned <- models.visited[1:n.to.report, ]
      }
      
      ## Return a list of results. This ends the "nw = 1" case.
      return(list("causal" = theta, "se" = sigma, "model.matrix" = model.matrix.to.be.returned,
                  "snp.probs" = snp.probs, "w" = w))
      
      
    } else {
      
      
      ## ---------- NW > 1 ---------- ##
      
      
      ## If multiple w values need to be run, loop over them.
      
      ## Store per-w results here.
      all.thetas <- rep(0, nw)
      all.sigmas <- rep(0, nw)
      all.probs <- matrix(0, nw, P)
      
      ## This will be used to store the top-model matrix for the best JAM-MR run.
      current.min.stderr <- Inf
      current.top.models <- 0
      
      ## This stores warning signs when JAM visits the null model.
      warning.signs <- rep(0, nw)
      
      ## Get the loop going.
      for (jj in 1:nw) {
        
        print(paste("---------- w = ", w[jj], " ----------"))
        
        ## Run JAM.
        if (is.G.provided) {
          jam.results <- JAM( marginal.betas = bx, n = N1, X.ref = G.matrix, initial.model = initial.model,
                              model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                              n.iter = iter, seed = jam.seed + jj - 1, trait.variance = trait.var, mrloss.w = w[jj], 
                              mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
        } else {
          jam.results <- JAM( marginal.betas = bx, n = N1, cor.ref  = G.cor, mafs.ref = eafs, initial.model = initial.model,
                              model.space.prior = list("a" = jam.model.priors[1], "b" = jam.model.priors[2], "Variables" = snp.names), 
                              n.iter = iter, seed = jam.seed + jj - 1, trait.variance = trait.var, mrloss.w = w[jj], 
                              mrloss.function = loss.function, mrloss.marginal.by = by, mrloss.marginal.sy = sy )
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
            warning.signs[jj] <- 2
          }
          
        } else {
          
          ## If the null model was one of many, discard it and use the rest for causal effect estimation.
          n.per.model <- rowSums(models.visited[, - (P + 1)])
          if ( min(n.per.model) == 0) {
            models.visited.backup <- models.visited
            zero.prob <- models.visited[which(n.per.model == 0), P + 1]
            models.visited <- models.visited[- which(n.per.model == 0), ]
            models.visited[, P + 1] <- models.visited[, P + 1] / (1 - zero.prob)
            if (zero.prob > 0.2) warning.signs[jj] <- 1
          }
          
        }
        
        
        ## For numerical stability we treat separately the case where JAM-MR only visits one model.
        if (nrow(models.visited) > 1) {
          
          ## ----- JAM-MR visits multiple models ----- ##
          
          ## Compute model probabilities and inclusion probabilities per SNP.
          model.probs <- models.visited[, P + 1]
          snp.probs <- rep(0, P)
          for (ii in 1:P) snp.probs[ii] <- sum(model.probs * models.visited[, ii])
          
          ## For each model, compute IVW estimates and standard errors.
          alleffects <- rep(0, nrow(models.visited))
          allerrors <- rep(0, nrow(models.visited))
          
          for (i in 1:nrow(models.visited)) {
            
            ## Extract the current model's SNPs, thetas and sigmas.
            current.model.snps <- unname(which(models.visited[i, - (P + 1)] == 1)) 
            current.bx <- unname(bx[current.model.snps])
            current.by <- by[current.model.snps]
            current.sy <- sy[current.model.snps]
            
            ## Truncation and random-effects are not currently supported for correlated data, so do IVW instead.
            current.invOmega <- solve( (current.sy %*% t(current.sy)) * Rho[current.model.snps, current.model.snps] )
            allerrors[i] <- sqrt( as.numeric( solve( t(current.bx) %*% current.invOmega %*% current.bx ) ) )
            alleffects[i] <- as.numeric( t(current.bx) %*% current.invOmega %*% current.by ) * allerrors[i]^2
            
          }
          
          ## If some (but not many) models returned "NaN", make the necesary adjustments.
          if (sum(model.probs[allerrors != "NaN"]) > 0.5) {
            model.probs <- model.probs[allerrors != "NaN"] / sum(model.probs[allerrors != "NaN"])
            alleffects <- alleffects[allerrors != "NaN"]
            allerrors <- allerrors[allerrors != "NaN"]
          }
          
          ## Combine into an overall estimator and its standard error.
          all.thetas[jj] <- sum(model.probs * alleffects)
          all.sigmas[jj] <- sqrt( sum( model.probs * allerrors^2 ) + sum( model.probs * alleffects * (alleffects - all.thetas[jj]) ) )
          all.probs[jj, ] <- snp.probs
          
          
        } else { 
          
          ## ----- JAM-MR visits only a single model ----- ##
          
          ## Extract the current model's SNPs, thetas and sigmas.
          single.model.snps <- unname(which(models.visited[, - (P + 1)] == 1))
          all.probs[jj, ] <- models.visited[, - (P + 1)]
          current.bx <- unname(bx[single.model.snps])
          current.by <- by[single.model.snps]
          current.sy <- sy[single.model.snps]
          
          ## Truncation and random-effects are not currently supported for correlated data, so do IVW instead.
          current.invOmega <- solve( (current.sy %*% t(current.sy)) * Rho[single.model.snps, single.model.snps] )
          all.sigmas[jj] <- sqrt( as.numeric( solve( t(current.bx) %*% current.invOmega %*% current.bx ) ) )
          all.thetas[jj] <- as.numeric( t(current.bx) %*% current.invOmega %*% current.by ) * all.sigmas[jj]^2
          
        }
        
        ## ----- Round up the per-w iteration ----- ##
        
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
      
      ## ----- Round up the nw > 1 case ----- ##
      
      ## Select the best w value.
      best.w <- w[which.min(all.sigmas)]
      theta <- all.thetas[which.min(all.sigmas)]
      sigma <- min(all.sigmas, na.rm = TRUE)
      snp.probs <- all.probs[which.min(all.sigmas), ]
      
      ## If the best w corresponds to a JAM-MR run that visited the null model, warn.
      if (warning.signs[which.min(all.sigmas)] == 2) {
        warning("JAM-MR did not select any genetic variants. Genetic associations with the risk factor are probably very weak.")
      }
      if (warning.signs[which.min(all.sigmas)] == 1) {
        warning("Genetic associations with the risk factor are probably fairly weak.")
      }
      
      ## Return a list of results. This ends the nw > 1 case.
      return(list("causal" = theta, "se" = sigma, "model.matrix" = current.top.models, "snp.probs" = snp.probs, 
                  "w" = best.w, "all.causal" = all.thetas, "all.se" = all.sigmas, "all.w" = w, "all.probs" = all.probs))
      
    }
    
  }
  
  ## Goodbye!
  
}



## Auxiliary functions for fitting the multiplicative random-effects truncated normal model.

## Compute the multiplicative random-effects truncated normal likelihood.
.mt.lik <- function (pars, l, u, thetas, sigmas) {
  dPhi <- pnorm( (u - pars[1]) / (sigmas * sqrt(pars[2])) ) - pnorm( (l - pars[1]) / (sigmas * sqrt(pars[2])) )
  if (min(dPhi) <= 0) lik <- Inf else  lik <- 1/2 * length(thetas) * log(pars[2]) + 1/2 * sum( log(sigmas^2) ) + 1/2 * 1 / pars[2] * sum( (thetas - pars[1])^2 / sigmas^2 ) + sum( log( dPhi ) )
  lik
}

## Compute the corresponding Gradient.
.mt.grad <- function (pars, l, u, thetas, sigmas) {
  Phi1 <- pnorm( (u - pars[1]) / (sigmas * sqrt(pars[2])) )
  Phi2 <- pnorm( (l - pars[1]) / (sigmas * sqrt(pars[2])) )
  phi1 <- dnorm( (u - pars[1]) / (sigmas * sqrt(pars[2])) )
  phi2 <- dnorm( (l - pars[1]) / (sigmas * sqrt(pars[2])) )
  if (u == Inf) arg.u <- rep(0, length(sigmas)) else arg.u <- (u - pars[1]) / (sigmas * sqrt(pars[2]))
  if (l == - Inf) arg.l <- rep(0, length(sigmas)) else arg.l <- (l - pars[1]) / (sigmas * sqrt(pars[2]))
  g1 <- - sum( (phi1 - phi2) / (sigmas * sqrt(pars[2]) * (Phi1 - Phi2)) ) - sum( (thetas - pars[1]) / (sigmas^2 * pars[2]) )
  g2 <- 1/2 * length(thetas) / pars[2] - 1/2 * 1 / (pars[2])^2 * sum( (thetas - pars[1])^2 / sigmas^2) - 1/2 *  1/pars[2] * sum( (arg.u * phi1 - arg.l * phi2) / (Phi1 - Phi2) )
  c(g1, g2)
}

## Compute the corresponding Hessian.
.mt.hess <- function (pars, l, u, thetas, sigmas) {
  Phi1 <- pnorm( (u - pars[1]) / (sigmas * sqrt(pars[2])) )
  Phi2 <- pnorm( (l - pars[1]) / (sigmas * sqrt(pars[2])) )
  phi1 <- dnorm( (u - pars[1]) / (sigmas * sqrt(pars[2])) )
  phi2 <- dnorm( (l - pars[1]) / (sigmas * sqrt(pars[2])) )
  if (u == Inf) arg.u <- rep(0, length(sigmas)) else arg.u <- (u - pars[1]) / (sigmas * sqrt(pars[2]))
  if (l == - Inf) arg.l <- rep(0, length(sigmas)) else arg.l <- (l - pars[1]) / (sigmas * sqrt(pars[2]))
  H11 <- sum(1 / (sigmas^2 * pars[2])) - sum( 1 / (sigmas^2 * pars[2]) * ( (arg.u * phi1 - arg.l * phi2) * (Phi1 - Phi2) + (phi1 - phi2)^2 ) / (Phi1 - Phi2)^2 )
  H12 <- sum( (thetas - pars[1]) / (sigmas^2 * pars[2]^2) ) - sum( 1 / (2 * sigmas * pars[2]^(3/2)) * ( ((arg.u^2 - 1) * phi1 - (arg.l^2 - 1) * phi2) * (Phi1 - Phi2) + (phi1 - phi2) * (arg.u * phi1 - arg.l * phi2) ) / (Phi1 - Phi2)^2 )
  H22 <- - length(thetas) / (2 * pars[2]^2) + 1 / pars[2]^3 * sum( (thetas - pars[1])^2 / sigmas^2 ) - 1 / (4 * pars[2]^2) * sum( ( (arg.u^3 * phi1 - arg.l^3 * phi2) / (Phi1 - Phi2) ) - 3 * ( (arg.u * phi1 - arg.l * phi2) / (Phi1 - Phi2) ) + ( (arg.u * phi1 - arg.l * phi2)^2 / (Phi1 - Phi2)^2 ) )
  matrix(c(H11, H12, H12, H22), 2, 2)
}

