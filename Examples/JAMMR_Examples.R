
### Using an artificial dataset.
set.seed(21)
P <- 50
bx.sim <- rnorm(P, 0.15, 0.05)
sx.sim <- abs(rnorm(P, 0.02, 0.008))
which.pl <- rbinom(P, 1, 0.3)   ## Which SNPs are pleiotropic.
by.sim <- bx.sim * 0.5 + rnorm(P, 0, 0.05)   ## True causal effect is 0.5.
by.sim[which.pl == 1] <- by.sim[which.pl == 1] + rep(0.25, sum(which.pl))   ## Pleiotropic effects.
sy.sim <- abs(rnorm(P, 0.04, 0.01))

### Visualize the data.
plot(bx.sim, by.sim, type = "p", pch = 19, xlim = c(0.00, 0.30), ylim = c(-0.10, 0.45), main = "Bx-By Plot")
for (i in 1:50) lines(c(bx.sim[i], bx.sim[i]), c(by.sim[i] - 1.96 * sy.sim[i], by.sim[i] + 1.96 * sy.sim[i]))
for (i in 1:50) lines(c(bx.sim[i] - 1.96 * sx.sim[i], bx.sim[i] + 1.96 * sx.sim[i]), c(by.sim[i], by.sim[i]))
abline(a = 0, b = median(by.sim / bx.sim), col = "brown", lty = 2)

### Run JAM-MR with a single value of w and independent data.
jammr.results1 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 10000, eafs = rep(0.3, 50),
                         iter = 1e6, w = 10000, jam.seed = 22)
jammr.results1$causal
jammr.results1$se

jammr.results1$snp.probs   ## Probabilities of inclusion per SNP.
jammr.results1$snp.probs[which.pl == 0]   ## Probabilities for valid SNPs.
jammr.results1$snp.probs[which.pl == 1]   ## Probabilities for pleiotropic SNPs.

### Run with multiple values of w.
jammr.results1.1 <- JAMMR(bx.sim, sx.sim, by.sim, sy.sim, N1 = 10000, eafs = rep(0.3, 50),
                           iter = 1e6, n.grid = 6, jam.seed = 22)
jammr.results1.1$causal
jammr.results1.1$se
jammr.results1.1$all.causal
jammr.results1.1$all.se
jammr.results1.1$snp.probs



### Using a real dataset from the MendelianRandomization package.
library(MendelianRandomization)

### Run existing MR methods.
mr_ivw( mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse) )
mr_egger( mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse) )
mr_median( mr_input(bx = ldlc, bxse = ldlcse, by = chdlodds, byse = chdloddsse) )

### Run JAM-MR.
jammr.results <- JAMMR(ldlc, ldlcse, chdlodds, chdloddsse, N1 = 17723, eafs = lipid_eaf, iter = 1e6,
                        n.grid = 6, grid.limits = c(100, 20000), n.models = 0, jam.seed = 4189)
jammr.results$causal
jammr.results$se
jammr.results$snp.probs
