#' @include R2BGLiMS.R
NULL

#' Wrapper for JAM (Joint Analysis of Marginal statistics). See the Vignette.
#' 
#' @export
#' @title JAM (Joint Analysis of Marginal statistics)
#' @name JAM
#' @inheritParams R2BGLiMS
#' @param trait.variance NB EXPERIMENTAL: Invokes a slightly different summary statistic likelihood, which requires an estimate of
#' the trait variance to be provided. (Not yet implemented to work with mJAM or enumeration)
#' @param full.mcmc.sampling By default JAM only samples models, since the parameters can be analytically integrated out.
#' If you would like to force JAM to perform full Reversible Jump MCMC of models and parameters (effects and residual), then
#' set this to TRUE. The posterior summaries can be seen using \code{\link{PrettyResultsTable}}. Note that this option is
#' forced to false if inference via model enumeration is requested by setting enumerate.up.to.dim>0.
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
#' @example Examples/JAM_Examples.R

JAM <- function(
  marginal.betas=NULL,
  n=NULL,
  ns.each.ethnicity=NULL,
  X.ref=NULL,
  cor.ref=NULL,
  mafs.ref=NULL,
  model.space.priors=NULL,
  beta.priors=NULL,
  beta.prior.partitions=NULL,
  g.prior=TRUE,
  tau=NULL,
  xtx.ridge.term=0,
  enumerate.up.to.dim=0,
  full.mcmc.sampling=FALSE,
  n.iter=1e6,
  n.mil.iter=NULL,
  thinning.interval=NULL,
  seed=NULL,
  extra.arguments=NULL,
  initial.model=NULL,
  save.path=NULL,
  debug=FALSE,
  max.model.dim=-1,
  burnin.fraction = 0.5,
  trait.variance = NULL,
  mrloss.w=0,
  mrloss.function="variance",
  mrloss.marginal.by=NULL,
  mrloss.marginal.sy=NULL,
  mafs.if.independent=NULL,
  extra.java.arguments=NULL
) {
  
  ##################################################################
  ### --- Error messages (moved from R2BGLiMS on 2016-01-29) --- ###
  ##################################################################
  
  # --- X.ref checks (if mafs.if.independent not provided)
  
  if (!is.null(mafs.if.independent)) {
    if (enumerate.up.to.dim > 0) { stop("When specifiying indepdent SNPs not yet possible to do enumeration.")}
  }

  # --- Marginal betas
  if (is.null(marginal.betas)) { stop("For analysis with JAM you must provide a vector of marginal summary statistics") }

  # -- N
  if (!is.null(ns.each.ethnicity)) {
    if (!is.vector(ns.each.ethnicity)) stop("ns.each.ethnicity must be a vector")
    n = sum(ns.each.ethnicity)
  } else {
    if (is.null(n)) { stop("You must specificy the number of individuals the summary statistics were calculated in.") }
  }

  #####################################################
  ### --- Whether to perform full MCMC sampling --- ###
  #####################################################
  
  if (full.mcmc.sampling) {
    which.blgims.jam.method <- "JAM_MCMC"
  } else {
    which.blgims.jam.method <- "JAM"
  }
  
  ############################
  ### --- Run R2BGLiMS --- ###
  ############################
  
  if (enumerate.up.to.dim > 0) {
    if (
      (is.list(X.ref) & length(X.ref>1)) |
      (is.list(cor.ref) & length(cor.ref>1))
    ) {stop("For enumeration please just supply one region at a time")}
    # Force not using JAM_MCMC
    which.blgims.jam.method <- "JAM"
    ### --- Enumeration
    n.iter <- 1 # Set to minimum number of iterations
    ######################################
    ### --- Enumeration - 1 region --- ###
    ######################################
    results <- R2BGLiMS(
      likelihood=which.blgims.jam.method,
      marginal.betas=marginal.betas,
      n=n,
      X.ref=X.ref,
      cor.ref=cor.ref,
      mafs.ref=mafs.ref,
      ns.each.ethnicity=ns.each.ethnicity,
      model.space.priors=model.space.priors,
      g.prior=g.prior,
      tau=tau,
      xtx.ridge.term=xtx.ridge.term,
      enumerate.up.to.dim=enumerate.up.to.dim,
      n.iter=n.iter,
      n.mil.iter=n.mil.iter,
      thinning.interval=thinning.interval,
      seed=seed,
      extra.arguments=extra.arguments,
      initial.model=initial.model,
      save.path=save.path,
      debug=debug,
      max.model.dim=max.model.dim,
      burnin.fraction = burnin.fraction,
      trait.variance = trait.variance,
      mrloss.w = mrloss.w,
      mrloss.function = mrloss.function,
      mrloss.marginal.by = mrloss.marginal.by,
      mrloss.marginal.sy = mrloss.marginal.sy,
      mafs.if.independent = mafs.if.independent,
      extra.java.arguments=extra.java.arguments
    )
  } else {
    ### --- RJMCMC    
    results <- R2BGLiMS(
      likelihood=which.blgims.jam.method,
      marginal.betas=marginal.betas,
      n=n,
      X.ref=X.ref,
      cor.ref=cor.ref,
      mafs.ref=mafs.ref,
      ns.each.ethnicity=ns.each.ethnicity,
      model.space.priors=model.space.priors,
      beta.priors=beta.priors,
      beta.prior.partitions=beta.prior.partitions,
      g.prior=g.prior,
      tau=tau,
      xtx.ridge.term=xtx.ridge.term,
      enumerate.up.to.dim=enumerate.up.to.dim,
      n.iter=n.iter,
      n.mil.iter=n.mil.iter,
      thinning.interval=thinning.interval,
      seed=seed,
      extra.arguments=extra.arguments,
      initial.model=initial.model,
      save.path=save.path,
      debug=debug,
      max.model.dim=max.model.dim,
      burnin.fraction = burnin.fraction,
      trait.variance = trait.variance,
      mrloss.w = mrloss.w,
      mrloss.function = mrloss.function,
      mrloss.marginal.by = mrloss.marginal.by,
      mrloss.marginal.sy = mrloss.marginal.sy,
      mafs.if.independent = mafs.if.independent,
      extra.java.arguments=extra.java.arguments
    )
  }

  return(results)
}
