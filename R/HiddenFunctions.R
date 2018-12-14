#' Calculates a Bayes Facotor, give a prior and posterior probabiltiies. Check2
#' 
#' @param prior.prob prior probability
#' @param post.prob posterior probability
#' @return Bayes Factor
#' @author Paul Newcombe
.BayesFactor <- function(prior.prob, post.prob) {
  priorOdds <- prior.prob/(1-prior.prob)
  postOdds <- post.prob/(1-post.prob)
  bayesFactor <- postOdds/priorOdds
  return(bayesFactor)
}

#' Calculates a Beta-Binomial prior probability for a SPECIFIC model. From Bottolo et al.
#' This is not the prior on a particular dimension (that would required the binomial co-efficient).
#' 
#' @param k dimension of specific model
#' @param n total number of covariates
#' @param a beta-binomial hyper parameter a
#' @param b beta-binomial hyper parameter b
#' @return Probability
#' @author Paul Newcombe
.BetaBinomialProbabilitySpecificModel <- function(k, n, a, b) {
  beta( (k+a), (n-k+b) )/beta(a,b)
}

#' Calculates prior probabability of causality for a particular variable, when a Poisson prior is used for model space
#' @param V total number of variables
#' @param mu model space mean
#' @param n.var number of specific variables (default is 1)
#' @return prior probability for a single variable
#' @author Paul Newcombe
.ModelSpaceSpecProb <- function(V, mu, n.var=1) {
  ### Normalising constant for truncated Poisson
  poiDenom <- 0
  for (v in 0:V) {	
    poiDenom <- poiDenom + (mu^v)*exp(-mu)/gamma(v+1)
  }
  ## m - 1 places to fill (saince 1st place taken by SNP in question), from M-1 SNPs (SNP in question, that is in 1st place, is not an option)
  priorProb <- 0
  for (v in 1:V) {	# Add across all ways can select specific nvar
    p.v.selected <- (mu^v)*exp(-mu)/(gamma(v+1)*poiDenom) # prob of v variables included under the poisson
    # Union of rule for P(A U B U C etc)
    p.atleast.one.given.v <- 0
    i.max <- min(v,n.var)
    for (i in 1:i.max) {
      p.atleast.one.given.v <- p.atleast.one.given.v + (choose(n.var,i)*choose( (V-i) , (v-i) )/choose( V , v ))*((-1)^(i+1))
    }
    priorProb <- priorProb + p.v.selected*p.atleast.one.given.v
  }
  
  return(priorProb)
}

#' Replaces the outcome variable with residuals from a linear regression on the confounders. I.e. removes effect
#' of confounders for the conjugate models.
#' @param data matrix or data.frame containing outcome and confounders
#' @param outcome.var name of outcome.variable
#' @param confounders list of confounders
#' @return data original dataset with the confounders removed and the outcome replaces with adjusted residuals
#' @author Paul Newcombe
.ConfounderAdjustedResiduals <- function(data, outcome.var, confounders) {
  
  # Get confounder adjusted residuals
  formula <- formula(paste(outcome.var,"~",paste(confounders,collapse="+")))
  lm.res <- lm(formula, data)
  
  # Replace outcome variable with confounder adjusted residuals
  data[,outcome.var] <- lm.res$residuals
  
  # Remove confounders from data matrix
  keep.cols <- colnames(data)[!colnames(data)%in%confounders]
  data <- data[,keep.cols]
  
  # Return dataset
  return(data)
}

#' Calculates prior probababilities for x, and >=x causal variables, when a truncated Poisson prior is used for model space
#' @param V total number of variables
#' @param mu model space mean
#' @return prior matrix whereby 1st row is prior probabilities for =x causal variables, and 2nd row for >=x variables (x = 1,...V)
#' @author Paul Newcombe
.ModelSpaceProbs <- function(V, mu) {
  # Normalising constant for truncated Poisson
  poiDenom <- 0
  for (v in 0:V) {	
    poiDenom <- poiDenom + (mu^v)*exp(-mu)/gamma(v+1)
  }
  
  # =x
  numPresPriorProb <- c(1:(V+1))
  for (v in 1:V) {
    # Prior prob of  >= m SNPs for m = 1,..M
    numPresPriorProb[v] <- (mu^v)*exp(-mu)/(gamma(v+1)*poiDenom)
  }
  numPresPriorProb[V+1] <- exp(-mu)/(gamma(1)*poiDenom)
  
  # >=x
  geqPriorProb <- c(1:(V+1))
  for (v in 1:V) {
    geqPriorProb[v] <- sum(numPresPriorProb[v:V]) 
  }	
  geqPriorProb[V+1] <- sum(numPresPriorProb[0:V])
  
  # combining
  return.mat <- rbind(numPresPriorProb,geqPriorProb)
  colnames(return.mat) <- paste(c(c(1:V),0),sep="")
  rownames(return.mat) <- c("=",">=")
  
  return(return.mat)
}

#' Calculates FDR thresholds according to Algorithm 18.3 on P.689 of elemnts of statistical learning by Hastie et al.
#' This is a conservative caclulation - where resolution is not sufficient
#' 
#' @param obs.probs Posterior probabilities from analysis of actual data
#' @param permuted.probs Posterior probabilties from analysis of permuted outcome analyses
#' @param target.fdrs Vector of FDRs to estimate Posterior Probability thresholds for (defaults to 1%, 5% and 10%)
#' @param n.cuts.order Order of magnitude for vector length of possible thresholds to explore
#' @return Matrix of tagret FDRs, their estimated posterior probability thresholds, and the estimated FDR at each
#' threshold
#' @author Paul Newcombe
.GetFdrThresholds <- function(
  obs.probs,
  permuted.probs,
  target.fdrs=c(0.01,0.05,0.1),
  n.cuts.order=5
) {
  if (is.matrix(permuted.probs)) {
    n.permute <- ncol(permuted.probs)
  } else {
    n.permute <- 1
  }
  cuts <- seq(from=min(obs.probs), to=max(obs.probs), length.out=10^n.cuts.order)
  cuts.fdr <- sapply(cuts, function(x) sum(permuted.probs>=x)/(n.permute*sum(obs.probs>=x)) )
  
  fdr.thresholds <- matrix(NA,length(target.fdrs),3,
                           dimnames=list(paste(target.fdrs),c("FDR", "PostProbThreshold", "FDR_hat")))
  for (fdr in target.fdrs) {
    fdr.diff <- fdr - cuts.fdr
    cut.index <- which(fdr.diff==min(fdr.diff[fdr.diff>=0]))
    if (length(cut.index)>1) {
      cut.index <- cut.index[1]
    }
    fdr.thresholds[paste(fdr),"FDR"] <- fdr
    fdr.thresholds[paste(fdr),"PostProbThreshold"] <- cuts[cut.index]
    fdr.thresholds[paste(fdr),"FDR_hat"] <- cuts.fdr[cut.index]
  }
  
  return(fdr.thresholds)
}

#' Calulcates approximate posterior probabiltiies from a Reversible Jump Results object
#' for an anlaysis which has been run with enumeration turned on for the beginning.
#' @export
#' @title Enumertated approximate posterior probabilities
#' @name .ApproxPostProbsByModelEnumeration
#' @param results R2BGLiMS results object
#' @param enumerate.up.to.dim Maxmimum dimension to truncate space too (can be less than or equal to
#' what was passed to R2BGLiMS). If left NULL the value passed to R2BGLiMS will be used.
#' @return A list containing all model, marginal and dimension approximate posterior probabilities 
#' @author Paul Newcombe
.ApproxPostProbsByModelEnumeration <- function(enumerated.model.likelihood.scores, model.space.priors, enumerate.up.to.dim) {
  
  # Determine model space prior
  if ("a" %in% names(model.space.priors[[1]])) {
    model.space.prior <- "beta.binom"
    a <- model.space.priors[[1]]$a
    b <- model.space.priors[[1]]$b
  } else if ("Rate" %in% names(model.space.priors[[1]])) {
    model.space.prior <- "poisson"
    poisson.lambda <- length(model.space.priors[[1]]$Variables)*model.space.priors[[1]]$Rate
  }
  
  # Setup hyper-parmaters and options
  vars <- model.space.priors[[1]]$Variables
  P <- length(vars)
  
  # Setup approx probs and model dimensions
  approx.probs <- enumerated.model.likelihood.scores
  approx.probs$Prob <- NULL
  model.dims <- unlist(lapply(strsplit(split="_AND_",as.character(enumerated.model.likelihood.scores$Model)),length))
  model.dims[1] <- 0
  approx.probs <- approx.probs[which(model.dims<=enumerate.up.to.dim),]
  model.dims <- model.dims[which(model.dims<=enumerate.up.to.dim)]
  
  # Null model
  if (model.space.prior=="beta.binom") {
    prior.prob.0 <- .BetaBinomialProbabilitySpecificModel(k=0,n=P,a=a,b=b)    
  } else if (model.space.prior=="poisson") {
    prior.prob.0 <- dpois(0, poisson.lambda)
  }
  approx.probs[approx.probs$Model=="Null","Prob"] <- approx.probs[approx.probs$Model=="Null","PosteriorScore"] + log(prior.prob.0)
  
  # Single SNP model
  if (model.space.prior=="beta.binom") {
    prior.prob.1 <- .BetaBinomialProbabilitySpecificModel(k=1,n=P,a=a,b=b)
  } else if (model.space.prior=="poisson") {
    prior.prob.1 <- dpois(1, poisson.lambda)
  }  
  approx.probs[model.dims==1,"Prob"] <- approx.probs[model.dims==1,"PosteriorScore"] + log(prior.prob.1)
  
  # Dual SNP models
  if (enumerate.up.to.dim>=2) {
    if (model.space.prior=="beta.binom") {
      prior.prob.2 <- .BetaBinomialProbabilitySpecificModel(k=2,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.2 <- dpois(2, poisson.lambda)
    }    
    approx.probs[model.dims==2,"Prob"] <- approx.probs[model.dims==2,"PosteriorScore"] + log(prior.prob.2)
  }
  
  # Triple SNP models
  if (enumerate.up.to.dim>=3) {
    if (model.space.prior=="beta.binom") {
      prior.prob.3 <- .BetaBinomialProbabilitySpecificModel(k=3,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.3 <- dpois(3, poisson.lambda)
    }
    approx.probs[model.dims==3,"Prob"] <- approx.probs[model.dims==3,"PosteriorScore"] + log(prior.prob.3)
  }
  
  # Quadruple SNP models
  if (enumerate.up.to.dim>=4) {
    if (model.space.prior=="beta.binom") {
      prior.prob.4 <- .BetaBinomialProbabilitySpecificModel(k=4,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.4 <- dpois(4, poisson.lambda)
    }
    approx.probs[model.dims==4,"Prob"] <- approx.probs[model.dims==4,"PosteriorScore"] + log(prior.prob.4)
  }
  
  # Quintuple SNP models
  if (enumerate.up.to.dim>=5) {
    if (model.space.prior=="beta.binom") {
      prior.prob.5 <- .BetaBinomialProbabilitySpecificModel(k=5,n=P,a=a,b=b)
    } else if (model.space.prior=="poisson") {
      prior.prob.5 <- dpois(5, poisson.lambda)
    }
    approx.probs[model.dims==5,"Prob"] <- approx.probs[model.dims==5,"PosteriorScore"] + log(prior.prob.5)
  }
  
  # Normalise model probs
  probs.norm <- approx.probs$Prob - mean(approx.probs$Prob)
  model.probs <- exp(probs.norm) /sum(exp(probs.norm))
  names(model.probs) <- approx.probs$Model
  
  # Infer marginal probs DO NOT USE GREP TO PICK OUT MODELS - can have intersecting names
  models.snps.list <- strsplit(x=names(model.probs),split= "_AND_" )
  marg.probs <- rep(0,P)
  names(marg.probs) <- vars
  for (v in vars) {
    marg.probs[v] <- sum(model.probs[unlist( lapply(models.snps.list, function(x) v %in% x) )])
  }
  
  # Infer dimension probabilities
  dim.probs <- c(1:enumerate.up.to.dim)
  for (i in 1:enumerate.up.to.dim) {
    dim.probs[i] = sum(model.probs[model.dims==i])
  }
  
  return(list("model.probs"=model.probs, "marg.probs"=marg.probs, "dim.probs"=dim.probs))
}

#' Reads in a .txt data file, formatted for my Java program, and returns as a list
#' @export
#' @title Read formatted data file
#' @name .ReadData
#' @param data.file data file to be read
#' @return dataRead a list containing 
#' covariates: matrix of covariate values
#' disease: binary outcome vector
#' n: number of subjects
#' R: number of clusters
#' startRJ: 1st variable included in RJ 
#' V: number of variables 
#' var.names: names 
#' @author Paul Newcombe
.ReadData <- function(data.file) {
  
  dataRead <- list(NA)	
  
  # Always have markers, between study var and Loglike
  dataRead$model <- scan(data.file, nlines = 1)
  dataRead$V <- scan(data.file, skip=1, nlines = 1)
  dataRead$var.names <- read.table(data.file,skip=2,nrows=1)
  dataRead$var.names <- as.character(unlist(dataRead$var.names))
  dataRead$startRJ <- scan(data.file, skip = 3, nlines = 1)
  dataRead$n <- scan(data.file, skip = 4, nlines = 1)
  dataRead$R <- scan(data.file, skip = 5, nlines = 1)
  dataRead$covariates <- matrix(scan (data.file, skip = 6, nlines=dataRead$n), ncol=dataRead$V, byrow=TRUE)
  colnames(dataRead$covariates) <- dataRead$var.names
  if (dataRead$R>0) {
    dataRead$randInts <- scan(data.file, skip = (6+dataRead$n), nlines=1) # Skips arguments and covariate rows
    dataRead$R <- max(dataRead$randInts)
    # Random intercept matrix helps with linear predictor
    dataRead$randIntMat <- matrix(0,dataRead$n,dataRead$R)
    for (i in 1:dataRead$n ) {
      dataRead$randIntMat[i,dataRead$randInts[i]] <- 1
    }
  }
  
  # Disease
  dataRead$disease <- scan(data.file, skip = (6+as.integer(dataRead$R>0)+dataRead$n), nlines=1)		# Skips arguments and covariate rows
  if (dataRead$model=="Weibull") {
    dataRead$times <- scan(data.file, skip = (7+as.integer(dataRead$R>0)+dataRead$n), nlines=1)		# Skips arguments and covariate rows    
  }
  
  # Return list
  return(dataRead)
}

#' This function is used by JAMPred_Step2AdjustmentAndPredictions 
.ReweightSnpEffectsAccordingToBlockCorrelations <- function(iteration, snp.blocks, X.ref) {
  
  # --- 1) Determine which blocks have atleast 1 SNP selected
  non.null.blocks <- unlist(lapply(snp.blocks, function(block) sum(iteration[block]!=0)>0 ))
  
  # --- 2) Calculate relative weightings of blocks by applying JAM
  block.weights <- rep(0, length(snp.blocks)) # If all blocks are null, each is weighted 0
  if (sum(non.null.blocks) == 1) {
    # Only one non-null block it is given weight 1
    block.weights[non.null.blocks] <- 1
  } else if (sum(non.null.blocks) > 1) {
    # If multiple non-null blocks apply JAM to get relative weights of the block scores
    # --- Construct score for each block in reference data
    Block.scores.ref <- Reduce(
      "cbind",
      lapply(snp.blocks[non.null.blocks], function(x) X.ref[,x] %*% iteration[x]))
    # --- Re-apply JAM to block scores (all marginal effects are 1)
    block.weights[non.null.blocks] <- as.vector(JAM_PointEstimates(marginal.betas = rep(1, sum(non.null.blocks)), X.ref = Block.scores.ref ))
  }
  
  # --- 3) Re-weight SNP effects according to block weights
  reweighted.snp.effects <- unlist(Reduce(
    "c",
    sapply(c(1:length(non.null.blocks)), function(b) {iteration[snp.blocks[[b]]]*block.weights[b] })))
  
  # --- 4) Return re-weighted SNP effects
  return(reweighted.snp.effects)
}
