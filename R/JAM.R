#' @include R2BGLiMS.R
NULL

#' Wrapper for JAM (Joint Analysis of Marginal statistics). See the Vignette.
#' 
#' @export
#' @title JAM (Joint Analysis of Marginal statistics)
#' @name JAM
#' @inheritParams R2BGLiMS
#' @param full.mcmc.sampling By default JAM only samples models, since the parameters can be analytically integrated out.
#' If you would like to force JAM to perform full Reversible Jump MCMC of models and parameters (effects and residual), then
#' set this to TRUE. The posterior summaries can be seen using \code{\link{PrettyResultsTable}}. Note that this option is
#' forced to false if inference via model enumeration is requested by setting enumerate.up.to.dim>0.
#' @param n.cases If the marginal.betas contain log-Odds Ratios, please specify the number of cases
#' with this option, so that JAM can calculate the case proportion in order to invoke an approximate 
#' transformation between the linear and logistic scales.
#' 
#' @return A Reversible Jump results object is returned. This is a list of two elements: "args" which records various modelling
#' arguments used in the analysis, and "results" - a matrix containing all saved posterior samples from the analysis. Columns
#' correspond to parameters and each row contains values from a particular itertation, with 0 indicating exclusion from the model.
#' 
#' The function \code{\link{PrettyResultsTable}} can be used to print summary posterior results for all parameters. Other functions
#' for summarising results are listed under "see also".
#' 
#' @seealso Summary results are stored in the slot posterior.summary.table. See also
#' \code{\link{PrettyResultsTable}} and \code{\link{ManhattanPlot}}. For posterior model space
#' summaries see \code{\link{TopModels}}. For convergence checks
#' plots see \code{\link{ChainPlots}} and \code{\link{AutocorrelationPlot}}.
#' 
#' @author Paul Newcombe
#' 
#' @example Examples/R2BGLiMS_Examples.R

JAM <- function(
  marginal.betas=NULL,
  n=NULL,
  X.ref=NULL,
  model.space.priors=NULL,
  g.prior=TRUE,
  tau=NULL,
  enumerate.up.to.dim=0,
  full.mcmc.sampling=FALSE,
  n.mil=1,
  seed=1,
  n.cases=NULL,
  extra.arguments=NULL
) {
  
  #######################################
  ### --- Set tau if not provided --- ###
  #######################################
  
  if (is.null(tau)) {
    cat("\nSetting tau to max(n,P^2)\n")
    tau <- max(n, length(marginal.betas)^2)
  }
  
  #####################################################
  ### --- Whether to perform full MCMC sampling --- ###
  #####################################################
  
  if (full.mcmc.sampling) {
    which.blgims.jam.method <- "JAM_MCMC"
  } else {
    which.blgims.jam.method <- "JAM"    
  }
  
  #######################################
  ### --- Logistic transformation --- ###
  #######################################
  
  if (!is.null(n.cases)) {
    cat("\nLog-Odds Ratios were provided.\n")
    cat("\nInvoking the logistic -> linear transformation.")
    phi <- n.cases/n
    marginal.betas <- marginal.betas*phi*(1-phi)
    if (is.null(extra.arguments)) {
      extra.arguments=list(
        "GaussianResidualVarianceInvGammaPrior_a" = n,
        "GaussianResidualVarianceInvGammaPrior_b" = (n-1)*phi*(1-phi)
      )      
    } else {
      extra.arguments[["GaussianResidualVarianceInvGammaPrior_a"]] = n
      extra.arguments[["GaussianResidualVarianceInvGammaPrior_b"]] = (n-1)*phi*(1-phi)
    }
  }
  
  ############################
  ### --- Run R2BGLiMS --- ###
  ############################
  
  if (enumerate.up.to.dim > 0) {
    which.blgims.jam.method <- "JAM" # Full MCMC sampling is disabled for model enumeration
    ### --- Enumeration
    n.mil <- 0.000001 # Set to minimum number of iterations
    if (!is.list(X.ref)|length(X.ref)==1) {
      ######################################
      ### --- Enumeration - 1 region --- ###
      ######################################
      results <- R2BGLiMS(
        likelihood=which.blgims.jam.method,
        marginal.betas=marginal.betas,
        n=n,
        X.ref=X.ref,
        model.space.priors=model.space.priors,
        g.prior=g.prior,
        tau=tau,
        enumerate.up.to.dim=enumerate.up.to.dim,
        n.mil=n.mil,
        seed=seed,
        extra.arguments=extra.arguments
      )
    } else {
      ##############################################
      ### --- Enumeration - Multiple regions --- ###
      ##############################################
      model.space.priors.ld.block <- model.space.priors
      for (ld.block in 1:length(X.ref)) {
        cat("\n--------------------------\n")
        cat("\n--------------------------\n")
        cat("\nEnumeration for block",ld.block,"\n")        
        cat("\n--------------------------\n")
        cat("\n--------------------------\n")
        vars.ld.block <- colnames(X.ref[[ld.block]])
        model.space.priors.ld.block$Variables <- vars.ld.block
        results.r2bglims.g <- R2BGLiMS(
          likelihood=which.blgims.jam.method,
          marginal.betas=marginal.betas[vars.ld.block],
          n=n,
          X.ref=X.ref[[ld.block]],
          model.space.priors=model.space.priors.ld.block,
          g.prior=g.prior,
          tau=tau,
          enumerate.up.to.dim=enumerate.up.to.dim,
          n.mil=n.mil,
          seed=seed,
          extra.arguments=extra.arguments
        )
        if (ld.block==1) {
          results <- results.r2bglims.g
          results@posterior.summary.table <- results@posterior.summary.table[vars.ld.block,] # Only display SNPs
          results@n.covariate.blocks.for.jam <- length(X.ref)
          results@enumerated.posterior.inference <- list(results@enumerated.posterior.inference) # Setup as a list
        } else {
          results@posterior.summary.table <- rbind( # Append extra rows to posterior summary table
            results@posterior.summary.table,
            results.r2bglims.g@posterior.summary.table[vars.ld.block,]
          )
          results@model.space.priors[[ld.block]] <- model.space.priors.ld.block # Append extra model space components
          results@enumerated.posterior.inference[[ld.block]] <- results.r2bglims.g@enumerated.posterior.inference # Add enumeration to the list
        }
      }
    }
  } else {
    ### --- RJMCMC    
    results <- R2BGLiMS(
      likelihood=which.blgims.jam.method,
      marginal.betas=marginal.betas,
      n=n,
      X.ref=X.ref,
      model.space.priors=model.space.priors,
      g.prior=g.prior,
      tau=tau,
      enumerate.up.to.dim=enumerate.up.to.dim,
      n.mil=n.mil,
      seed=seed,
      extra.arguments=extra.arguments      
    )
  }

  return(results)
}
