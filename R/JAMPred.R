#' This function runs JAMPred by analysing many SNP blocks independently in parallel.
#' @export
#' @title JAMPred
#' @name JAMPred
#' @inheritParams JAMPred_SplitIntoPositiveDefiniteBlocks
#' @param marginal.betas Vector of named one-at-a-time SNP effects. NB: This must be a named vector; it is where.
#' the SNP names are derived from. These could be from linear regressions, or log-Odds Ratios. 
#' NB: If providing log-ORs, you must also provide marginal.logor.ses and p.cases.training so
#' that a linear transformaiton can be applied.
#' @param n.training The sample size in which the marginal.betas were calculated.
#' @param marginal.logor.ses IF marginal log-ORs were provided, please provide a named vector of standard errors here.
#' These are necessary to convert them to a linear scale. Do not set if providing linear regression estimates.
#' @param p.cases.training IF marginal log-ORs were provided, please provide the proportion of cases in the sample in
#' which they were calculated
#' @param ref.geno Reference genotype matrix which will be used to calculate SNP-SNP correlations. Individual's genotype 
#' must be coded as a numeric risk allele count 0/1/2. Non-integer values reflecting imputation uncertaimnty may be given. 
#' NB: The risk allele coding MUST correspond to that used in marginal.betas and column names must match the names of
#' marginal.betas elements.
#' @param total.snps.genome.wide Total number of SNPs being analysed across the genome. Used to determine the level of
#' prior sparsity.
#' @param n.cores Number of CPU cores to parralelise the SNP blocks over. If unspecified the number available will be detected. NB: 
#' On Windows systems this will be forced to 1, since mclapply does not yet support forking.
#' @param n.mil Number of million iterations to run. In the paper we found 0.2 was sufficient to reach adequeate
#' convergence.
#' @param beta.binom.b.lambda Spartisty tuning parameter. Suggest a range of values is tried, e.g. 0.01, 0.1, 1, 10.
#' See the paper for more detail. Default is 1.
#' @param beta.binom.a First hyper-parameter of the beta-binomial prior on model sparsity. See the paper for more detail. Default 1 can usually
#' be used.
#' @param snps.blocks OPTIONAL: A partitioning of the SNPs into blocks, created using the function 
#' \code{\link[R2BGLiMS]{JAMPred_SplitIntoPositiveDefiniteBlocks}}. This is automatically generated 
#' if not provided.
#' @param seed An integer specifying the RJMCMC seed. If not set a random number will be used every time this is run.
#' @param thinning.factor Determines every ith iteration to use from the RJMCMC sample. The more thinning that is applied
#' the quicker the analysis will be, for the same number of iterations. Leaving at the default should be fine.
#' @param effect.sd.uniform.prior Upper and lower uniform hyper-parameters for the prior on the standard deviation
#' of SNP effects. Default is what we used in the paper, a Uniform(0.05, 2) distribution.
#' @param residual.var.invgamma.prior Hyper-parameters for the inversegamma prior on the residual variance.
#' Default is what we used in the paper, a Inverse-Gamma(0.01, 0.01) distribution.
#' @param save.path Optional path for JAM to store and save its temporary files in. Can help with debugging (default NULL).
#' @param debug An option passed to JAM requesting more verbose output (default FALSE).
#' 
#' @return A JAMPred results object, which is a list including as elements step1.posterior.mean.snp.weights 
#' (which do not adjust for long range LD) and step2.posterior.mean.snp.weights (which do adjust for long 
#' range LD). These SNP weights can be used to generate predictions from individual level genotype data. 
#' 
#' @seealso See \code{\link{JAMPred_SplitIntoPositiveDefiniteBlocks}}; this is the function used to create an initial 
#' partitioning of SNPs into blocksm which satisfy the property of all being positvie definite with respect to the 
#' reference genotype matrix (also required by JAMPred). By default this function is used to generate the snps.blocks
#' list input, if not provided by the user.
#' 
#' @author Paul Newcombe
#' 
#' @example Examples/JAMPred_Examples.R

JAMPred <- function(
  marginal.betas = NULL,
  n.training = NULL,
  marginal.logor.ses = NULL,
  p.cases.training = NULL,
  ref.geno = NULL,
  total.snps.genome.wide = NULL,
  n.cores = NULL,
  n.mil = 0.2,
  beta.binom.b.lambda = 1,
  beta.binom.a = 1,
  snps.blocks = NULL,
  initial.block.size = 100,
  initial.snps.blocks = NULL,
  seed = NULL,
  thinning.factor = round(1e3/n.mil),
  effect.sd.uniform.prior = c(0.05,2),
  residual.var.invgamma.prior = c(0.01,0.01),
  save.path=NULL,
  debug=FALSE
) {
  
  library(parallel)
  
  if(is.null(names(marginal.betas))) stop("marginal.betas must be a NAMED vector. It is where the SNP names are
                                          derived from.")
  
  # --- Deriving the names of SNPs from marginal.betas
  snps <- names(marginal.betas)
  
  # --- Map logistic to linear effects IF p.cases.training and ses are provided
  if ( !is.null(marginal.logor.ses) | !is.null(p.cases.training) ) {
    
    if (is.null(marginal.logor.ses) | is.null(p.cases.training) ) stop("If providing log-ORs must provide both their SEs and the proportion of cases.")
    cat("\nLog-ORs provided. Applying linear transformation.\n")
    marginal.betas <- JAM_LogisticToLinearEffects(
      log.ors = marginal.betas,
      log.or.ses = marginal.logor.ses[snps],
      mafs = apply(ref.geno[,snps],MAR=2,mean)/2,
      n = n.training,
      p.cases = p.cases.training
    )    
  }
  
  # --- Split SNPs into blocks
  if (is.null(snps.blocks)) {
    cat("\n.Generating a SNP partitioning of step 1 blocks, since one was not provided.\n")
    snps.blocks <- JAMPred_SplitIntoPositiveDefiniteBlocks(
      snps=snps,
      initial.block.size = initial.block.size,
      initial.snps.blocks = initial.snps.blocks,
      ref.geno = ref.geno[,snps]
    )
  }
  
  # --- Derive parallel block indices
  if (.Platform$OS.type=="windows") {
    n.cores <- 1
    cat("Forcing to the use of 1 core since mclapply can not fork jobs in Windows.\n")
  } else {
    if (is.null(n.cores)) {
      n.cores <- detectCores(logical=FALSE)
      cat("Detected",n.cores,"available to use.\n")
    }
  }
  parallel.block.indices <- JAMPred_ParallelBlockIndices(
    n.cores = n.cores,
    n.blocks.to.analyse = length(snps.blocks)
  )
  
  # --- JAM Step 1: Run JAM MCMC treating blocks as independent
  n.cores <- length(parallel.block.indices) # Determin number of cores to use
  step1.res <- unlist( mclapply(c(1:n.cores), 
                          function(i)
                            JAM(
                              marginal.betas = marginal.betas[unlist(snps.blocks[parallel.block.indices[[i]]])],
                              X.ref=lapply(snps.blocks[parallel.block.indices[[i]]], function(b) as.matrix(ref.geno[,b])),
                              n=n.training,
                              model.space.priors = list("a"=beta.binom.a,"b"=beta.binom.b.lambda*total.snps.genome.wide,"Variables"=unlist(snps.blocks[parallel.block.indices[[i]]])),
                              beta.prior.partitions = list(list("UniformA"=effect.sd.uniform.prior[1], "UniformB"=effect.sd.uniform.prior[2], "Variables"=unlist(snps.blocks[parallel.block.indices[[i]]]))),
                              full.mcmc.sampling = TRUE,
                              n.mil = n.mil,
                              thinning.interval = round(n.mil*thinning.factor),
                              seed = seed,
                              extra.arguments = 
                                list("GaussianResidualVarianceInvGammaPrior_a"=residual.var.invgamma.prior[1],
                                     "GaussianResidualVarianceInvGammaPrior_b"=residual.var.invgamma.prior[2]),
                              save.path=save.path,
                              debug=debug
                            ), 
                          mc.cores=n.cores), recursive=F)
  
  # --- Step 1
  step1.posterior.mean.snp.weights <- unlist(
    lapply(
      step1.res,
      function(i) i@posterior.summary.table[i@model.space.priors[[1]]$Variables,"Mean"]
  ))
  
  # --- Step 2
  step2.posterior.mean.snp.weights <- .ReweightSnpEffectsAccordingToBlockCorrelations(
    step1.posterior.mean.snp.weights, snps.blocks, ref.geno
  )

  # --- Combine results
  step1.res <- list(
    "rjmcmc.res"=step1.res,
    "snps.blocks"=snps.blocks,
    "parallel.block.indices"=parallel.block.indices,
    "step1.posterior.mean.snp.weights"=step1.posterior.mean.snp.weights,
    "step2.posterior.mean.snp.weights"=step2.posterior.mean.snp.weights,
    "snps"=snps
  )
  
  return(step1.res)
}
