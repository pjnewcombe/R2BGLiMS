#' Generates a table of the top models from a Reversible Jump Results object.
#' @export
#' @title Table of top models
#' @name TopModels
#' @inheritParams ManhattanPlot
#' @param n.top.models number of top models to include in the table (default 20)
#' @return Ordered table of best models, with corresponding posterior probabilities. If JAM
#' was run with multiple covariate blocks, a list of model tables is returned.
#' @author Paul Newcombe
#' @example Examples/TopModels_Examples.R 
TopModels <- function(
  results,
  n.top.models=20,
  remove.empty.cols=TRUE) {
  
  ######################
  ### --- Errors --- ###
  ######################
  
  if (n.top.models < 2) {stop("n.top.models must be atleast 2")}

	### Get list of all models differently depending on MCMC vs enumeration
	if (results@enumerate.up.to.dim==0) {
	  #########################
	  # --- MCMC was done --- #
	  #########################
	  vars.to.include <- unlist(lapply(results@model.space.priors, function(x) x$Variables))
	  all.models <- apply(
	    results@mcmc.output[,vars.to.include],
	    MAR=1,
	    function(r) paste( as.integer(r!=0), collapse="_")
	  )
	  ### --- Make table
	  models.table <- sort(table(all.models),d=T)/nrow(results@mcmc.output) # Tabulates to unique models
	  n.top.models <- min(n.top.models, length(models.table))
	  models.table <- models.table[1:n.top.models]
	  models.tab.str <- strsplit(names(models.table), split="_")
	  models.tab <- NULL
	  for (m in 1:length(models.tab.str)) {
	    models.tab <- rbind(models.tab, as.integer(models.tab.str[[m]]) )
	  }
	  colnames(models.tab) <- vars.to.include
	  models.tab <- cbind(models.tab, "Post Prob"=models.table)
	  rownames(models.tab) <- NULL
	  # --- Remove empty columns
	  if (remove.empty.cols) {
	    empty.cols <- NULL
	    for (v in colnames(models.tab)[colnames(models.tab)!="Post Prob"]) {
	      if ( sum(models.tab[,v])==0) {
	        empty.cols <- c(empty.cols,which(colnames(models.tab)==v))
	      }
	    }
	    if (!is.null(empty.cols)) {
	      models.tab <-models.tab[,-empty.cols]        
	    }
	  }
	} else {
	  ################################
	  # --- Enmueration was done --- #
	  ################################
	  models.tabs.list <- list()
	  for (ld.block in 1:results@n.covariate.blocks.for.jam) {
	    if (results@n.covariate.blocks.for.jam==1) {
	      enumerated.model.probs <- results@enumerated.posterior.inference$model.probs
	    } else {
	      enumerated.model.probs <- results@enumerated.posterior.inference[[ld.block]]$model.probs
	    }
	    all.models.avail <- names(enumerated.model.probs)
	    all.unique.vars <- unique(unlist(strsplit(all.models.avail, "_AND_")))[-1]
	    models.tab <- matrix(0,length(all.models.avail), length(all.unique.vars))
	    colnames(models.tab) <- all.unique.vars
	    all.models.avail <- paste("_",all.models.avail,"_",sep="") # To guard against overlapping names
	    for (v in all.unique.vars) {
	      models.tab[grep(paste("_",v,"_",sep=""), all.models.avail),v] <- 1
	    }
	    models.tab <- cbind(models.tab, "Post Prob" = enumerated.model.probs)
	    rownames(models.tab) <- NULL
	    models.tab <- models.tab[order(models.tab[,"Post Prob"], decreasing = TRUE),]
	    n.top.models <- min(n.top.models, nrow(models.tab))
	    models.tab <- models.tab[1:n.top.models,]
	    # --- Remove empty columns
	    if (remove.empty.cols) {
	      empty.cols <- NULL
	      for (v in colnames(models.tab)[colnames(models.tab)!="Post Prob"]) {
	        if ( sum(models.tab[,v])==0) {
	          empty.cols <- c(empty.cols,which(colnames(models.tab)==v))
	        }
	      }
	      if (!is.null(empty.cols)) {
	        models.tab <-models.tab[,-empty.cols]        
	      }
	    }
	    # --- Add to the list incase there are multiple LD blocks
	    models.tabs.list[[ld.block]] <- models.tab
	    if (results@n.covariate.blocks.for.jam > 1) {
	      cat("\n Table for LD block",ld.block,"finished\n")
	    }
	  }
	  if (results@n.covariate.blocks.for.jam>1) {
	    models.tab <- models.tabs.list # Only return a list if there are multiple LD blocks
	  }
	}
  
  # --- Return
	return(models.tab)		
}	
