library(R2BGLiMS)

#############################################
# --- Example 1) One region with RJMCMC --- #
#############################################

# By default 1 million iterations are run
library(R2BGLiMS) # Load package
data(JAM_Example) # Load example data
jam.results <- JAM(
  marginal.betas=marginal.betas[snps.region1],
  X.ref=X.ref.region1,
  model.space.prior = list("a"=1, "b"=length(snps.region1), "Variables"=snps.region1),
  n=n)
jam.results@posterior.summary.table # SNP 5 was simulated as causal
ManhattanPlot(jam.results)

###############################################################
# --- Example 2) Two regions simultaneously (with RJMCMC) --- #
###############################################################

two.regions.results <- JAM(
  marginal.betas=marginal.betas,
  X.ref=list(X.ref.region1, X.ref.region2),
  model.space.prior = list(
    "a"=1, "b"=length(snps.region1)+length(snps.region2), 
    "Variables" = c(snps.region1, snps.region2)),
  n.mil=1,
  n=n)
two.regions.results@posterior.summary.table # SNPs 5 and 16 were simulated as causal
ManhattanPlot(two.regions.results)

