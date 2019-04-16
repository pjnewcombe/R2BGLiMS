# --- Load libraries
library(R2BGLiMS)
data("JAMPRED_Example")

##################################
# --- Binary outcome example --- #
##################################

# Run JAMPred
jam.pred.res <- JAMPred(
  marginal.betas = marginal.logors,
  n.training = n.training,
  marginal.logor.ses = marginal.logor.ses, # Default NULL - only necessary if passing log-ORs for a binary trait
  n.cases.training = n.cases.training, # Default NULL - only necessary if passing log-ORs for a binary trait
  ref.geno = data.validation[,snps], # Only the SNP data is used (only need to provide SNPs)
  total.snps.genome.wide = 500000, # Total SNPs genomewide being analysed
  n.mil = 0.2
)

# Generate predictions
out.of.sample.predictions <- 
  data.validation[,jam.pred.res$snps] %*% 
  t(t(jam.pred.res$step2.posterior.mean.snp.weights))

# Predictive r2
cor(out.of.sample.predictions, data.validation[,"d"])^2

######################################
# --- Continuous outcome example --- #
######################################

# Run JAMPred
jam.pred.res <- JAMPred(
  marginal.betas = marginal.ctsbetas,
  n.training = n.training,
  ref.geno = data.validation.cts[,snps], # Only SNP data is used (only need to provide SNPs)
  total.snps.genome.wide = 500000, # Total SNPs genomewide being analysed
  n.mil = 0.2
)

# Generate predictions
out.of.sample.predictions <- 
  0 + # NB: A trait mean could be set here (otherwise it is assumed the outcome is mean-centred)
  data.validation.cts[,jam.pred.res$snps] %*% 
  t(t(jam.pred.res$step2.posterior.mean.snp.weights))

# Predictive r2
cor(out.of.sample.predictions, data.validation.cts[,"y"])^2
