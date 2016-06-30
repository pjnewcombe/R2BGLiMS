#' Infer credible sets of predictors from a Reversible Jump Results object.
#' @export
#' @title Infer credible sets of predictors
#' @name CredibleSet
#' @inheritParams ManhattanPlot
#' @param credible.percentile.threshold Percentile to determine credible set for (default 0.95 for a 95\% Credible Set).
#' @return List of variables included in the credible set.
#' @author Paul Newcombe
CredibleSet <- function(
  results,
  credible.percentile.threshold=0.95) {
  
  ######################
  ### --- Errors --- ###
  ######################
  
  ### Error checks
  
	if (results@enumerate.up.to.dim==0) {
	  #########################
	  # --- MCMC was done --- #
	  #########################
	  vars.to.include.in.table <- unlist(lapply(results@model.space.priors, function(x) x$Variables))      
	  all.models <- apply(
	    results@mcmc.output[,vars.to.include.in.table],
	    MAR=1,
	    function(r) paste( as.integer(r!=0), collapse="_")
	  )
	  
	  # --- Tabulate models and threshold at credible set threshold
	  models.table <- sort(table(all.models),d=T)/nrow(results@mcmc.output) # Tabulates to unique models
	  cumulative.model.probs <- sapply(c(1:length(models.table)), function(i) sum(models.table[1:i]) )
	  last.model.ind <- min(which(cumulative.model.probs> (credible.percentile.threshold) ))
	  models.table <- models.table[1:last.model.ind]

	  # --- Add covariate names for resulting credible set
	  models.tab.str <- strsplit(names(models.table), split="_")
	  credible.set <- vars.to.include.in.table[sort(unique(unlist(lapply(models.tab.str, function(i) which(i=="1")))))]
	  
	} else {
	  ################################
	  # --- Enmueration was done --- #
	  ################################
	  credible.set <- list()
	  for (ld.block in 1:results@n.covariate.blocks.for.jam) {
	    if (results@n.covariate.blocks.for.jam==1) {
	      enumerated.model.probs <- results@enumerated.posterior.inference$model.probs
	    } else {
	      enumerated.model.probs <- results@enumerated.posterior.inference[[ld.block]]$model.probs
	    }
	    
	    enumerated.model.probs <- sort(enumerated.model.probs, dec=T)
	    cumulative.model.probs <- sapply(c(1:length(enumerated.model.probs)), function(i) sum(enumerated.model.probs[1:i]) )
	    last.model.ind <- min(which(cumulative.model.probs> (credible.percentile.threshold) ))
	    credible.set[[ld.block]] <- unique(unlist(strsplit(names(enumerated.model.probs[1:last.model.ind]), "_AND_")))
	  }
	}
  
  # --- Return
	return(credible.set)
}
