#' @include HiddenFunctions.R
NULL

#' Creates a summary results table from a Java MCMC results class object
#' @export
#' @title Summary results table
#' @name ResultsTable
#' @param results Reversible Jump results object from running \code{\link{R2MHRJ}}.
#' @param vars.to.include Optional character vector specifying a subset of variables
#' to restrict output to.
#' @param normalised.sds If variables were normalised (e.g. to facilitate use of a common
#' unknown prior on their effects), a named vector of their original sds can be provided
#' such that estimates returned may be interpreted on the original scale.
#' @param var.dictionary Alternative covariate names. This must be a character vector,
#' named wtih the current covariate names. This can consist of expression functions, e.g. to
#' achieve italics or superscripts when passed to a plot function.
#' @return Prints a summary results table 
#' @author Paul Newcombe
#' @example Examples/ResultsTable_Examples.R 
ResultsTable <- function(
  results,
  vars.to.include=NULL,
  var.dictionary=NULL) {
  
  # Update names if variable dictionary provided
  if (!is.null(var.dictionary)) {
    colnames(results$results)[colnames(results$results)%in%names(var.dictionary)] <- var.dictionary[colnames(results$results)[colnames(results$results)%in%names(var.dictionary)]]
  }
  
	# Data results
	results.table <- matrix(NA,ncol(results$results),8)
	colnames(results.table) = c(
    "PostProb",
    "Median",
    "CrI_Lower",
    "CrI_Upper",
    "Median_Present",
    "CrI_Lower_Present",
    "CrI_Upper_Present",
    "BF")
  rownames(results.table) <- colnames(results$results)
  
	# Model space means
  # Needs work!!!
	prior.probs <- rep(NA,ncol(results$results))
	names(prior.probs) <- rownames(results.table)
  vars.fix <- NULL
  if (results$args$startRJ>0) {
    vars.fix <- colnames(results$results)[(grep("alpha", colnames(results$results))+1):(grep("alpha", colnames(results$results))+results$args$startRJ)]
  }
	if (results$args$nRjComp==1) {
	  first.beta.rj <- grep("alpha", colnames(results$results))+results$args$startRJ+1
	  last.beta.rj <- first.beta.rj + results$args$V - 1
    if (results$args$ModelSpacePriorFamily=="Poisson") {
      mu.normalised <- .ModelSpaceSpecProb(
        (results$args$V-results$args$startRJ), results$args$poisson.mu)
      prior.probs[first.beta.rj:last.beta.rj] <- mu.normalised
    } else if (results$args$ModelSpacePriorFamily=="BetaBinomial") {
      prior.probs[first.beta.rj:last.beta.rj] <- results$args$beta.binom.a/(results$args$beta.binom.a+results$args$beta.binom.b)
    }
	} else {
	  results$args$modSpaceSplits <- grep("alpha", colnames(results$results))+results$args$modSpaceSplits+1
    if (results$args$ModelSpacePriorFamily=="Poisson") {
      mu.normalised <- results$args$poisson.mu
      for (c in 1:results$args$nRjComp) {
        mu.normalised[c] <- .ModelSpaceSpecProb(
          (results$args$modSpaceSplits[c+1]-results$args$modSpaceSplits[c]), results$args$poisson.mu[c])
        prior.probs[results$args$modSpaceSplits[c]:(results$args$modSpaceSplits[c+1]-1)] <- mu.normalised[c]
      }      
    } else if (results$args$ModelSpacePriorFamily=="BetaBinomial") {
      for (c in 1:results$args$nRjComp) {
        prior.probs[results$args$modSpaceSplits[c]:(results$args$modSpaceSplits[c+1]-1)] <- results$args$beta.binom.a[c]/(results$args$beta.binom.a[c]+results$args$beta.binom.b[c])
      }
    }
	}
  
	for (v in rownames(results.table) ) {
    if (!v %in% c(
      "LogWeibullScale",
      "alpha",
      paste("LogBetaPriorSd",c(1:results$args$nBetaHyperPriorComp),sep=""),
      "LogLikelihood",
      vars.fix)) {
      results.table[v,"PostProb"] <- length( results$results[,v][results$results[,v]!=0] ) / nrow(results$results)
      results.table[v,"BF"] <- .BayesFactor( prior.probs[v], results.table[v,"PostProb"])      
    }
		if (results$args$Likelihood %in% c("Weibull", "Logistic") ) {
      # Exonentiate log-HRs or log-ORs
      results.table[v,c("CrI_Lower", "Median", "CrI_Upper")] <- exp( quantile(results$results[,v],c(0.025, 0.5, 0.975)) )
      results.table[v,c("CrI_Lower_Present", "Median_Present", "CrI_Upper_Present")] <- exp( quantile(results$results[,v][results$results[,v]!=0],c(0.025, 0.5, 0.975) ) )
    }
	}
	
	return(results.table)	
}	
