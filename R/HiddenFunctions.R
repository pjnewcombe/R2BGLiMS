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
