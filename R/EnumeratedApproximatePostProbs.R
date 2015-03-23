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
#' @return A list containing all model, marginal and dimension approximate posterior probabilities 
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
  vars <- results$args$model.space.priors[[1]]$Variables
  P <- length(vars)

  # Setup approx probs and model dimensions
  approx.probs <- results$model.scores
  approx.probs$Prob <- NULL
  model.dims <- unlist(lapply(strsplit(split="AND",as.character(results$model.scores$Model)),length))
  model.dims[1] <- 0
  approx.probs <- approx.probs[which(model.dims<=max.dim),]
  model.dims <- model.dims[which(model.dims<=max.dim)]
  
  # Null model
  prior.prob.0 <- .BetaBinomialProbability(k=0,n=P,a=a,b=b)
  approx.probs[approx.probs$Model=="Null","Prob"] <- exp(approx.probs[approx.probs$Model=="Null","PosteriorScore"])*prior.prob.0
  
  # Single SNP model
  prior.prob.1 <- .BetaBinomialProbability(k=1,n=P,a=a,b=b)
  approx.probs[model.dims==1,"Prob"] <- exp(approx.probs[model.dims==1,"PosteriorScore"])*prior.prob.1
  
  # Dual SNP models
  if (max.dim>=2) {
    prior.prob.2 <- .BetaBinomialProbability(k=2,n=P,a=a,b=b)
    approx.probs[model.dims==2,"Prob"] <- exp(approx.probs[model.dims==2,"PosteriorScore"])*prior.prob.2      
  }
  
  # Triple SNP models
  if (max.dim>=3) {
    prior.prob.3 <- .BetaBinomialProbability(k=3,n=P,a=a,b=b)
    approx.probs[model.dims==3,"Prob"] <- exp(approx.probs[model.dims==3,"PosteriorScore"])*prior.prob.3      
  }

  # Quadruple SNP models
  if (max.dim>=4) {
    prior.prob.4 <- .BetaBinomialProbability(k=4,n=P,a=a,b=b)
    approx.probs[model.dims==4,"Prob"] <- exp(approx.probs[model.dims==4,"PosteriorScore"])*prior.prob.4
  }
  
  # Normalise model probs
  model.probs <- approx.probs$Prob/sum(approx.probs$Prob)  
  names(model.probs) <- approx.probs$Model
  
  # Infer marginal probs
  marg.probs <- rep(0,P)
  names(marg.probs) <- vars
  for (v in vars) {
    marg.probs[v] <- sum(model.probs[grep(v,names(model.probs))])
  }
  
  # Infer dimension probabilities
  dim.probs <- c(1:max.dim)
  for (i in 1:max.dim) {
    dim.probs[i] = sum(model.probs[model.dims==i])
  }
  
  return(list("model.probs"=model.probs, "marg.probs"=marg.probs, "dim.probs"=dim.probs))
}	
