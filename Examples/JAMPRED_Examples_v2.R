# --- Load libraries
library(R2BGLiMS)
data("JAMPRED_Example")

# In practice, we recommend parallelisation of different lambdas and chromosomes as different HPC jobs
# We then recommend parallelisation within a chromosome on the same compute node, using R's "parallel" package.
# This is because not every block requires it's own compute core; JAMPred is designed to analyse a number 
# of blocks jointly and in parallel using the same compute node.

# We recommend parallelisation of different lambdas and chromosomes as different HPC jobs. We then recommend
# using R's parallel package to help parallelise analysis of blocks within a chromosome.
# The synthetic example JAMPred datset, included in R2BGLiMS, consists of data from 2 chromosomes, with 200 SNPs 
# per chromosome. We split each chromosome into 4 blocks of 50 SNPs., to demonstrate how JAMPred is designed to analyse a number 
# of blocks jointly

# First, we analyse

##################################
# --- Binary outcome example --- #
##################################

# NB: We recommend running each lambda sparsity and chromosome as
# a seperate job on an HPC, rather than one after the other like in
# the nested loops below.

jam.pred.res <- list()
for (beta.binom.b.lambda in c(0.01, 0.1, 1, 10)) {
  jam.pred.res[[paste(beta.binom.b.lambda)]] <- list()
  for (chromosome in c(1:2)) {
    # Ideally, the data corresponding to a complete chromosome should be passed to every
    # JAMPred analysis, in order to account for LD across (potentially) the complete 
    # chromosome. JAMPred will internally parallelise the analysis using R's parallel library.
    jam.pred.res[[paste(beta.binom.b.lambda)]][[chromosome]] <- JAMPred(
      marginal.betas = marginal.logors[chromosome.snps[[chromosome]]],
      n.training = n.training,
      marginal.logor.ses = marginal.logor.ses, # Only necessary if passing log-ORs for a binary trait
      p.cases.training = n.cases.training/n.training, # Only necessary for a binary trait
      ref.geno = data.validation[,chromosome.snps[[chromosome]]],
      total.snps.genome.wide = 500000, # Total SNPs genomewide being analysed
      n.mil = 0.2,
      beta.binom.b.lambda = beta.binom.b.lambda,
      initial.block.size = 100, # The 200 SNPs on each chromosome are split into 4 blocks
      n.cores = 2 
      # n.cores=2 leads to two parallel analyses of 2 blocks each. Significantly more cores
      # would be available on an HPC, per node, perhaps 16 or 32, so this could be increased 
      # to spread the analysis of many more blocks over more cores. NB: It is not
      # necessary to analyse every block on it's own core - memory permitting JAMPred
      # is designed to jointly analyse a number of blocks on the same CPU core, as we demonstrate
      # here (two blocks are analysed per core). For the data on a particular chromosome, we We recommend simply specifying the number 
      # of cores here, and letting JAMPred divide
      # 100 SNP blocks across the cores.
    )
  }
}

save.image("/Users/paul/Downloads/JAMPredExample.RData")

# Generate predictions for each lambda
out.of.sample.predictions <- list()
for (beta.binom.b.lambda in c(0.01, 0.1, 1, 10)) {
  out.of.sample.predictions[[paste(beta.binom.b.lambda)]] <- 0
  for (chromosome in c(1:2)) {
    out.of.sample.predictions[[paste(beta.binom.b.lambda)]] <-
      out.of.sample.predictions[[paste(beta.binom.b.lambda)]] +
      data.validation[,jam.pred.res[[paste(beta.binom.b.lambda)]][[chromosome]]$snps] %*% 
      t(t(jam.pred.res[[paste(beta.binom.b.lambda)]][[chromosome]]$step2.posterior.mean.snp.weights))
  }
}

# Predictive r2
sapply(out.of.sample.predictions, function(preds) cor(preds,data.validation[,"d"])^2)
