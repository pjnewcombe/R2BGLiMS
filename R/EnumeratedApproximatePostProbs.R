#' @include HiddenFunctions.R
NULL

#' Calulcates approximate posterior probabiltiies from a Reversible Jump Results object
#' for an anlaysis which has been run with enumeration turned on for the beginning.
#' @export
#' @title Enumertated approximate posterior probabilities
#' @name EnumeratedApproxPostProbs
#' @param results R2BGLiMS results object
#' @param max.dim Maxmimum dimension to truncate space too (can be less than or equal to
#' what was passed to R2BGLiMS). If left NULL the value passed to R2BGLiMS will be used.
#' @return A list containing all model, marginal and dimension approximate posterior probabilities 
#' @author Paul Newcombe
EnumeratedApproxPostProbs <- function(results,max.dim=NULL) {
  
  # Determine model space prior
  if ("a" %in% names(results$args$model.space.priors[[1]])) {
    model.space.prior <- "beta.binom"
    a <- results$args$model.space.priors[[1]]$a
    b <- results$args$model.space.priors[[1]]$b
  } else if ("Rate" %in% names(results$args$model.space.priors[[1]])) {
    model.space.prior <- "poisson"
    poisson.lambda <- length(results$args$model.space.priors[[1]]$Variables)*results$args$model.space.priors[[1]]$Rate
  }
  
  # Setup hyper-parmaters and options
  if (is.null(max.dim)) {
    max.dim <- results$args$enumerateUpToDim
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
  if (model.space.prior=="beta.binom") {
    prior.prob.0 <- .BetaBinomialProbability(k=0,n=P,a=a,b=b)    
  } else if (model.space.prior=="poisson") {
    prior.prob.0 <- dpois(0, poisson.lambda)
  }
  approx.probs[approx.probs$Model=="Null","Prob"] <- approx.probs[approx.probs$Model=="Null","PosteriorScore"] + log(prior.prob.0)
  
  # Single SNP model
  if (model.space.prior=="beta.binom") {
    prior.prob.1 <- .BetaBinomialProbability(k=1,n=P,a=a,b=b)
  } else if (model.space.prior=="poisson") {
    prior.prob.1 <- dpois(1, poisson.lambda)
  }  
  approx.probs[model.dims==1,"Prob"] <- approx.probs[model.dims==1,"PosteriorScore"] + log(prior.prob.1)
  
  # Dual SNP models
  if (max.dim>=2) {
    if (model.space.prior=="beta.binom") {
      prior.prob.2 <- .BetaBinomialProbability(k=2,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.2 <- dpois(2, poisson.lambda)
    }    
    approx.probs[model.dims==2,"Prob"] <- approx.probs[model.dims==2,"PosteriorScore"] + log(prior.prob.2)
  }
  
  # Triple SNP models
  if (max.dim>=3) {
    if (model.space.prior=="beta.binom") {
      prior.prob.3 <- .BetaBinomialProbability(k=3,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.3 <- dpois(3, poisson.lambda)
    }
    approx.probs[model.dims==3,"Prob"] <- approx.probs[model.dims==3,"PosteriorScore"] + log(prior.prob.3)
  }

  # Quadruple SNP models
  if (max.dim>=4) {
    if (model.space.prior=="beta.binom") {
      prior.prob.4 <- .BetaBinomialProbability(k=4,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.4 <- dpois(4, poisson.lambda)
    }
    approx.probs[model.dims==4,"Prob"] <- approx.probs[model.dims==4,"PosteriorScore"] + log(prior.prob.4)
  }
  
  # Normalise model probs
  probs.norm <- approx.probs$Prob - mean(approx.probs$Prob)
  model.probs <- exp(probs.norm) /sum(exp(probs.norm))
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
