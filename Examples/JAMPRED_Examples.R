# --- Load libraries
library(R2BGLiMS)
data("JAMPred_Example")

###########################################################
# --- Simple demonstation of syntax for binary traits --- #
###########################################################
# NB: See third example for modelling different chromosomes and sparsities

# Run JAMPred
snps <- chromosome.snps[[1]] # Only use chromosome 1 data
jampred.res.bin <- JAMPred(
  marginal.betas = marginal.logors[snps],
  n.training = n.training,
  marginal.logor.ses = marginal.logor.ses, # Only necessary for a binary trait
  p.cases.training = n.cases.training/n.training, # Only necessary for a binary trait
  ref.geno = data.validation[,snps],
  total.snps.genome.wide = 500000, # Total SNPs across all chromosomes
  n.mil = 0.2,
  seed = 1 # For re-producibility. If not set a random seed is used
)

# Generate predictions
out.of.sample.predictions <- 
  data.validation[,jampred.res.bin$snps] %*% 
  jampred.res.bin$step2.posterior.mean.snp.weights

# Predictive r2
cor(out.of.sample.predictions, data.validation[,"d"])^2 # Should be 0.18

###########################################################
# --- Simple demonstation of syntax for binary traits --- #
###########################################################
# NB: See third example for modelling different chromosomes and sparsities

# Run JAMPred
snps <- chromosome.snps[[1]] # Only use chromosome 1 data
jampred.res.lin <- JAMPred(
  marginal.betas = marginal.ctsbetas[snps],
  n.training = n.training,
  ref.geno = data.validation.cts[,snps],
  total.snps.genome.wide = 500000, # Total SNPs across all chromosomes
  n.mil = 0.2,
  seed = 1 # For re-producibility. If not set a random seed is used
)

# Generate predictions
out.of.sample.predictions <- 
  0 + # NB: A trait mean could be set here (otherwise it is assumed the outcome is mean-centred)
  data.validation.cts[,jampred.res.lin$snps] %*% 
  jampred.res.lin$step2.posterior.mean.snp.weights

# Predictive r2
cor(out.of.sample.predictions, data.validation.cts[,"y"])^2 # Should be 0.15

##################################################################################
# --- Using JAMPred to model multiple chromosomes under different sparsities --- #
##################################################################################

# Run JAMPred, looping over sparsities and chromosomes
jampred.res.loop <- list() # List to hold all the results
for (lambda in c(0.01, 0.1, 1, 10)) { # Loop over lambda choices (see the paper)
  jampred.res.loop[[paste(lambda)]] <- list()
  for (chr in c(1:2)) { # Loop over chromosomes
    jampred.res.loop[[paste(lambda)]][[chr]] <- JAMPred(
      marginal.betas = marginal.logors[chromosome.snps[[chr]]],
      n.training = n.training,
      marginal.logor.ses = marginal.logor.ses,
      p.cases.training = n.cases.training/n.training,
      ref.geno = data.validation[,chromosome.snps[[chr]]],
      total.snps.genome.wide = 500000,
      n.mil = 0.2,
      beta.binom.b.lambda = lambda,
      initial.block.size = 100, # Default size of SNP blocks
      seed = 1 # For re-producibility. If not set a random seed is used
    )
  }
}

# Generate predictions for each lambda, summing predictive scores over chromosoms
out.of.sample.predictions <- list()
for (lambda in c(0.01, 0.1, 1, 10)) {
  out.of.sample.predictions[[paste(lambda)]] <- 0
  for (chr in c(1:2)) {
    out.of.sample.predictions[[paste(lambda)]] <-
      out.of.sample.predictions[[paste(lambda)]] +
      data.validation[,jampred.res.loop[[paste(lambda)]][[chr]]$snps] %*% 
      jampred.res.loop[[paste(lambda)]][[chr]]$step2.posterior.mean.snp.weights
  }
}

# Predictive r2
sapply(out.of.sample.predictions, function(preds) cor(preds,data.validation[,"d"])^2)
# Predictive correlation 0.37 and highest at lambda = 0.01
