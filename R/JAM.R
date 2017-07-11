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
#' @param n The size of the dataset in which the summary statistics were calculated. This must be specified if use.da.v2=TRUE.
#' @param n.cases If the marginal.betas contain log-Odds Ratios, please specify the number of cases
#' with this option, so that JAM can calculate the case proportion in order to invoke an approximate 
#' transformation between the linear and logistic scales.
#' @param use.da.v2 Whether to use Daniel Ahfock's new formulation of the marginal JAM model likelihood. NB: Requires specification of trait.variance.ref and n.
#' @param trait.variance.ref Reference estimate of the trait variance. Must specify if use.da.v2=TRUE.
#' 
#' @return An R2BGLiMS_Results class object is returned. See the slot 'posterior.summary.table' for a posterior summary of all parameters. 
#' See slot 'mcmc.output' for a matrix containing the raw MCMC output from the saved posterior samples (0 indicates a covariate is excluded 
#' from the model in a particular sample. Some functions for summarising results are listed under "see also".
#' 
#' @seealso Summary results are stored in the slot posterior.summary.table. See \code{\link{ManhattanPlot}} for a visual 
#' summary of covariate selection probabilities. For posterior model space summaries see \code{\link{TopModels}}. For 
#' convergence checks see \code{\link{ChainPlots}} and \code{\link{AutocorrelationPlot}}.
#' 
#' @author Paul Newcombe
#' 
#' @example Examples/R2BGLiMS_Examples.R

JAM <- function(
  marginal.betas=NULL,
  n=NULL,
  X.ref=NULL,
  model.space.priors=NULL,
  beta.priors=NULL,
  beta.prior.partitions=NULL,
  g.prior=TRUE,
  tau=NULL,
  enumerate.up.to.dim=0,
  full.mcmc.sampling=FALSE,
  n.iter=1e6,
  n.mil.iter=NULL,
  thinning.interval=NULL,
  seed=NULL,
  n.cases=NULL,
  extra.arguments=NULL,
  initial.model=NULL,
  save.path=NULL,
  max.model.dim=-1,
  use.da.v2=FALSE,
  trait.variance.ref=NULL,
  extra.java.arguments=NULL
) {
  
  ##################################################################
  ### --- Error messages (moved from R2BGLiMS on 2016-01-29) --- ###
  ##################################################################
  
  # --- X.ref
  
  # Force to list
  if (is.data.frame(X.ref)) {
    X.ref <- matrix(X.ref) # convert to matrix
  }
  if (!is.list(X.ref)) {
    X.ref <- list(X.ref) # convert to list if there is a single block
  }
  
  # Check Rank
  # for (ld.block in 1:length(X.ref)) {
  #  qr.decomp <- qr(X.ref[[ld.block]])
  #  if (qr.decomp$rank < ncol(X.ref[[ld.block]])) stop (
  #    paste("The reference matrix for block",ld.block,"/",length(X.ref),"is not full rank.
  #            Ideally a larger reference sample should be used, 
  #            or you could try pruning correlated SNPs.")
  #  )
  #}
    
  # Check format
  if (sum(unlist(lapply(X.ref, function(x) !is.numeric(x) )))>0) {stop("Reference genotype matrices must be numeric, coded as risk allele countsin the 0 to 2 range")}
  if (max(unlist(X.ref))>2 | min(unlist(X.ref)) < 0) {stop("Reference genotype matrices must be coded coded as risk allele counts in the 0 to 2 range")}
  for (ld.block in 1:length(X.ref)) {
    if (is.null(colnames(X.ref[[ld.block]]))) stop ("All columns of the X reference matrice(s) must be named, corresponding to SNP effects in marginal.betas")
  }
  if (use.da.v2 & is.null(trait.variance.ref)) {stop("If using Daniel Ahfock's reformulation, must specify an estimate of the trait variance.")}
  if (full.mcmc.sampling & use.da.v2){stop("Full MCMC sampling of models and parameters not yet implemented for Daniel Ahfock's formulation.")}

  # --- Marginal betas
  if (is.null(marginal.betas)) { stop("For analysis with JAM you must provide a vector of marginal summary statistics") }
  if (is.null(names(marginal.betas))) stop ("All effects in marginal.betas must be named, corresponding to columns in X.ref")
  if (sum(names(marginal.betas) %in% unlist(lapply(X.ref, colnames))) < length(marginal.betas)) {stop("Reference genotype matrices do not include all SNPs in the marginal.betas vector")}
  
  # -- N
  if (is.null(n)) { stop("You must specificy the number of individuals the summary statistics were calculated in.") }

  ######################################################
  ### --- Set yTy.ref for Daniel Ahfock's method --- ###
  ######################################################
  
  if (use.da.v2) {
    yTy.ref <- trait.variance.ref*n
  } else {
    yTy.ref <- NULL
  }
  
  #######################################################################
  ### --- Take subset of X.refs correpsonding to elements of beta --- ###
  #######################################################################
  
  for (ld.block in 1:length(X.ref)) {
    original.n.snps <- ncol(X.ref[[ld.block]])
    X.ref[[ld.block]] <- X.ref[[ld.block]][,colnames(X.ref[[ld.block]]) %in% names(marginal.betas)]
    final.n.snps <- ncol(X.ref[[ld.block]])
    if (final.n.snps!=original.n.snps) {
      cat((original.n.snps-final.n.snps),"extra SNPs removed from the reference matrix for LD block",ld.block,"\n")
    }
  }
  
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
    if (use.da.v2) {
      which.blgims.jam.method <- "JAMv2"
    } else {
      which.blgims.jam.method <- "JAM"
    }
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
    if (use.da.v2) { # Force not using JAM_MCMC
      which.blgims.jam.method <- "JAMv2"
    } else {
      which.blgims.jam.method <- "JAM"
    }
    ### --- Enumeration
    n.iter <- 1 # Set to minimum number of iterations
    if (!is.list(X.ref)|length(X.ref)==1) {
      ######################################
      ### --- Enumeration - 1 region --- ###
      ######################################
      results <- R2BGLiMS(
        likelihood=which.blgims.jam.method,
        marginal.betas=marginal.betas,
        n=n,
        X.ref=X.ref,
        yTy.ref=yTy.ref,
        n.for.jam=n,
        model.space.priors=model.space.priors,
        g.prior=g.prior,
        tau=tau,
        enumerate.up.to.dim=enumerate.up.to.dim,
        n.iter=n.iter,
        n.mil.iter=n.mil.iter,
        thinning.interval=thinning.interval,
        seed=seed,
        extra.arguments=extra.arguments,
        initial.model=initial.model,
        save.path=save.path,
        max.model.dim=max.model.dim,
        extra.java.arguments=extra.java.arguments
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
          yTy.ref=yTy.ref,
          n.for.jam=n,
          model.space.priors=model.space.priors.ld.block,
          g.prior=g.prior,
          tau=tau,
          enumerate.up.to.dim=enumerate.up.to.dim,
          n.iter=n.iter,
          n.mil.iter=n.mil.iter,
          thinning.interval=thinning.interval,
          seed=seed,
          extra.arguments=extra.arguments,
          initial.model=initial.model,
          save.path=save.path,
          max.model.dim=max.model.dim,
          extra.java.arguments=extra.java.arguments
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
      yTy.ref=yTy.ref,
      n.for.jam=n,
      model.space.priors=model.space.priors,
      beta.priors=beta.priors,
      beta.prior.partitions=beta.prior.partitions,
      g.prior=g.prior,
      tau=tau,
      enumerate.up.to.dim=enumerate.up.to.dim,
      n.iter=n.iter,
      n.mil.iter=n.mil.iter,
      thinning.interval=thinning.interval,
      seed=seed,
      extra.arguments=extra.arguments,
      initial.model=initial.model,
      save.path=save.path,
      max.model.dim=max.model.dim,
      extra.java.arguments=extra.java.arguments
    )
  }

  return(results)
}
