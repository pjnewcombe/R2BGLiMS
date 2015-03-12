#' @include HiddenFunctions.R
NULL

#' Calulcates approximate posterior probabiltiies from a Reversible Jump Results object
#' for an anlaysis which has been run with enumeration turned on for the beginning.
#' @export
#' @title Enumertated approximate posterior probabilities
#' @name EnumeratedApproxPostProbs
#' @param results R2BGLiMS results object
#' @param a Beta-binomial hyperparamter a - optional. If not specified the value used in
#' the original analysis is used.
#' @param b Beta-binomial hyperparamter b - optional. If not specified the value used in
#' the original analysis is used.
#' @param max.dim Maxmimum dimension to truncate space too (can be less than or equal to
#' what was passed to R2BGLiMS). If left NULL the value passed to R2BGLiMS will be used.
#' @return Table of posterior probabilities and Bayes Factors for different model sizes 
#' @author Paul Newcombe
EnumeratedApproxPostProbs <- function(results,a=NULL,b=NULL,max.dim=NULL) {
  
  # Setup hyper-parmaters and options
  if (is.null(a)) {
    a <- results$args$model.space.priors[[1]]$a
  }
  if (is.null(b)) {
    b <- results$args$model.space.priors[[1]]$b
  }
  if (is.null(max.dim)) {
    max.dim <- results$args$allModelScoresUpToDim
  }

  # Setup approx probs and model dimensions
  approx.probs <- results$model.scores
  approx.probs$Prob <- NULL
  model.dims <- unlist(lapply(strsplit(split="AND",as.character(results$model.scores$Model)),length))
  model.dims[1] <- 0
  approx.probs <- approx.probs[which(model.dims<=max.dim),]
  model.dims <- model.dims[which(model.dims<=max.dim)]
  
  # Null model
  prior.prob.0 <- .BetaBinomialProbability(k=0,n=length(snp.list),a=a,b=b)
  approx.probs[approx.probs$Model=="Null","Prob"] <- exp(approx.probs[approx.probs$Model=="Null","PosteriorScore"])*prior.prob.0
  
  # Single SNP model
  prior.prob.1 <- .BetaBinomialProbability(k=1,n=length(snp.list),a=a,b=b)
  approx.probs[model.dims==1,"Prob"] <- exp(approx.probs[model.dims==1,"PosteriorScore"])*prior.prob.1
  
  # Dual SNP models
  if (max.dim>=2) {
    prior.prob.2 <- .BetaBinomialProbability(k=2,n=length(snp.list),a=a,b=b)
    approx.probs[model.dims==2,"Prob"] <- exp(approx.probs[model.dims==2,"PosteriorScore"])*prior.prob.2      
  }
  
  # Triple SNP models
  if (max.dim>=3) {
    prior.prob.3 <- .BetaBinomialProbability(k=3,n=length(snp.list),a=a,b=b)
    approx.probs[model.dims==3,"Prob"] <- exp(approx.probs[model.dims==3,"PosteriorScore"])*prior.prob.3      
  }

  # Quadruple SNP models
  if (max.dim>=4) {
    prior.prob.4 <- .BetaBinomialProbability(k=4,n=length(snp.list),a=a,b=b)
    approx.probs[model.dims==4,"Prob"] <- exp(approx.probs[model.dims==4,"PosteriorScore"])*prior.prob.4
  }
  
  # Normalise
  approx.post.probs <- approx.probs$Prob/sum(approx.probs$Prob)  
  names(approx.post.probs) <- approx.probs$Model
  return(approx.post.probs)
}	
