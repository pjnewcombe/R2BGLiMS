# --- Load libraries
library(R2BGLiMS)
data("JAMMR_Example")

# APOSTOLOS SHORT DESCRIPTION OF LOADED DATA

########################################################################
### --- Run JAM-MR with a single value of w and independent data --- ###
########################################################################

# Implement the algorithm.
jammr.results1 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 5000, eafs = rep(0.3, 50),
                        trait.var = 5, iter = 1e6, w = 5000, jam.seed = 22)

# Diagnostics.
jammr.results1$causal   ## Approximately 0.5 which is the correct value.
jammr.results1$se
jammr.results1$se.adj

jammr.results1$snp.probs   ## Probabilities of inclusion per SNP.
jammr.results1$snp.probs[which.pl == 0]   ## Probabilities for valid SNPs.
jammr.results1$snp.probs[which.pl == 1]   ## Probabilities for pleiotropic SNPs.

jammr.results1$model.matrix   ## The top models visited by JAMMR.

# Plotting.
JAMMR_ManhattanPlot(jammr.results1)
JAMMR_SensitivityPlot(jammr.results1)

# Run without the trait variance.
jammr.results1.1 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 5000, eafs = rep(0.3, 50),
                          iter = 1e6, w = 5000, jam.seed = 22)
# Different results because trait variance estimate here was completely arbitrary.

#########################################################################
### --- Run JAM-MR with a single w and dependent genetic variants --- ###
#########################################################################

# Simulate a reference genetic matrix for correlated-SNP analysis.
G.sim <- matrix(rbinom(500000, 2, 0.3), 10000, 50)
# This is just for testing purposes. The genetic matrix is based
# on an independence assumption so the results should be very
# similar to those obtained previously.

# Run the algorithm.
jammr.results2 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 5000, G.matrix = G.sim,
                        trait.var = 5, iter = 1e6, w = 5000, jam.seed = 23)

# Diagnostics.
jammr.results2$causal
jammr.results2$se
jammr.results2$snp.probs

# Plotting.
JAMMR_ManhattanPlot(jammr.results2)
JAMMR_SensitivityPlot(jammr.results2)

#########################################################################
### --- Run JAM-MR with multiple values of w and independent data --- ###
#########################################################################

# Implement the algorithm.
jammr.results3 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 5000, eafs = rep(0.3, 50),
                        trait.var = 5, iter = 1e6, jam.seed = 24, 
                        w = c(0, 500, 1000, 2000, 5000, 10000, 20000, 50000))

# Diagnostics.
jammr.results3$causal
jammr.results3$se
jammr.results3$se.adj
jammr.results3$snp.probs
jammr.results3$w

jammr.results3$all.causal
jammr.results3$all.se
jammr.results3$all.w
jammr.results3$all.probs

# Plotting.
JAMMR_ManhattanPlot(jammr.results3)
JAMMR_ManhattanPlot(jammr.results3, only.best = TRUE)
JAMMR_SensitivityPlot(jammr.results3, adjust.se = TRUE, show.n.snps = TRUE)

########################################################################
### --- Run JAM-MR with multiple values of w and correlated data --- ###
########################################################################

# Implement the algorithm.
jammr.results4 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 5000, G.matrix = G.sim,
                        trait.var = 5, iter = 1e6, jam.seed = 24, 
                        w = c(0, 500, 1000, 2000, 5000, 10000, 20000, 50000))

# Diagnostics.
jammr.results4$causal
jammr.results4$se
jammr.results4$se.adj
jammr.results4$snp.probs
jammr.results4$w

jammr.results4$all.causal
jammr.results4$all.se
jammr.results4$all.w
jammr.results4$all.probs

# Plotting.
JAMMR_ManhattanPlot(jammr.results4)
JAMMR_ManhattanPlot(jammr.results4, only.best = TRUE)
JAMMR_SensitivityPlot(jammr.results4, adjust.se = TRUE, show.n.snps = TRUE)
