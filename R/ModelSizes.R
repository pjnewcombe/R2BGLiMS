#' @include HiddenFunctions.R
NULL

#' Calulcates posterior probabiltiies and Bayes Factors for different model sizes,
#' from a Reversible Jump Results object
#' @export
#' @title Posterior model size/dimension summary
#' @name ModelSizes
#' @inheritParams ResultsTable
#' @return Table of posterior probabilities and Bayes Factors for different model sizes 
#' @author Paul Newcombe
#' @example Examples/ModelSizes_Examples.R 
ModelSizes <- function(results) {

	if (results$args$nRjComp==1) {
		# Setup
		vars <- c((results$args$startRJ+1):results$args$V)
		prior.probs <- rep(1,length(vars))
		post.probs <- rep(1,length(vars))
		prior.probs.geq <- rep(1,length(vars))
		post.probs.geq <- rep(1,length(vars))
		bf.geq <- rep(1,length(vars))

		# Prior probabilities of different numbers of variables present
		for (v in c(1:length(vars))) {
		  if (results$args$ModelSpacePriorFamily=="Poisson") {
		    prior.probs[v] <- dpois(v, results$args$poisson.mu)
		  } else if (results$args$ModelSpacePriorFamily=="BetaBinomial") {
		    prior.probs[v] <- dbinom(v, length(vars), results$args$beta.binom.a/(results$args$beta.binom.a+results$args$beta.binom.b))
		  }
		}
		# Posterior probabilities of different numbers of variables present
		vars.present <- apply(results$results[,vars], MAR=1, function(r) sum(r!=0))
		for (v in c(1:length(vars))) {
		  post.probs[v] <- sum(vars.present==v)/nrow(results$results)
		}	
		### Prior probs, posterior probs and Bayes Factors for GEQ
		for (v in 1:length(vars) ) {
			prior.probs.geq[v] <- sum(prior.probs[v:length(vars)])
			post.probs.geq[v] <- sum(post.probs[v:length(vars)])
			bf.geq[v] <- .BayesFactor( prior.probs.geq[v], post.probs.geq[v] )
		}
	
		### Results table
		results.table <- cbind(round(prior.probs.geq,4),post.probs.geq,bf.geq)
		rownames(results.table) <- paste(">=",c(1:length(vars))," variants",sep = "")
		colnames(results.table) <- c("Prior prob","Post prob","Bayes Factor")		
	} else {
	  # Setup - use lists, each element corresponding to a model space component
		prior.probs <- list(NA)
		post.probs <- list(NA)
		prior.probs.geq <- list(NA)
		post.probs.geq <- list(NA)
		bf.geq <- list(NA)		
		results.table <- list(NA)		
		if (results$args$ModelSpacePriorFamily=="Poisson") {
		  mu.normalised <- results$args$poisson.mu      
		}
		
		for (c in 1:results$args$nRjComp) {
      # Setup for model space component
			vars <- c((results$modSpaceSplits[c]+1):results$modSpaceSplits[c+1])	
			prior.probs[[c]] <- rep(1,length(vars))
			post.probs[[c]] <- rep(1,length(vars))
			prior.probs.geq[[c]] <- rep(1,length(vars))
			post.probs.geq[[c]] <- rep(1,length(vars))
			bf.geq[[c]] <- rep(1,length(vars))

			# Prior probs
			for (v in c(1:length(vars)) ) {
        if (results$args$ModelSpacePriorFamily=="Poisson") {
          prior.probs[[c]][v] <- dpois(v,results$args$poisson.mu[c])
        } else if (results$args$ModelSpacePriorFamily=="BetaBinomial") {
          prior.probs[[c]][v] <- dbinom(v, length(vars), results$args$beta.binom.a[c]/(results$args$beta.binom.a[c]+results$args$beta.binom.b[c]))          
        }
			}
      # Posterior probs
			vars.present <- apply(results$results[,vars], MAR=1, function(r) sum(r!=0))
			for (v in c(1:length(vars))) {
				post.probs[[c]][v] <- sum(vars.present==v)/nrow(results$results)
			}	
			
			### Prior probs, posterior probs and Bayes Factors for GEQ
			for (v in 1:length(vars) ) {
				prior.probs.geq[[c]][v] <- sum(prior.probs[[c]][v:length(vars)])
				post.probs.geq[[c]][v] <- sum(post.probs[[c]][v:length(vars)])
				bf.geq[[c]][v] <- .BayesFactor( prior.probs.geq[[c]][v], post.probs.geq[[c]][v] )
			}
			
			### Results table
			results.table[[c]] <- cbind(round(prior.probs.geq[[c]],4),post.probs.geq[[c]],bf.geq[[c]])
			rownames(results.table[[c]]) <- paste(">=",c(1:length(vars))," variants",sep = "")
			colnames(results.table[[c]]) <- c("Prior prob","Post prob","Bayes Factor")		
			
		}
	}
	return(results.table)	
}	
