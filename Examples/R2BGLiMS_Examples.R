library(R2BGLiMS)

# NB: See ?JAM for JAM examples.
#
# Examples below are:
# 1) Logistic regression
# 2) Weibull regression survival analysis
# 3) Case-cohort logistic regression (with hazard ratio estimation using Prenctice weighted Cox regression)
# 4) Linear regression
# 5) Linear regression using a conjugate posterior (with coefficients analytically intergreated out)
# 6) Fixed priors on effects, rather than the default common prior with unknown variance
# 7) Using two model space partitions

##########################################
# --- Example 1) Logistic regression --- #
##########################################

utils::data(biopsy, package = "MASS") # Example logistic dataset
covariate.names <- paste0("V",c(1:9))
# Recommend standardising predictors to justify default common hierarchical prior on effects
for (v in covariate.names) {biopsy[,v] <- scale(biopsy[,v])} 
biopsy$class <- as.integer(biopsy$class) - 1
biopsyResults <- R2BGLiMS( # Takes a few minutes to run
  likelihood="Logistic",
  data=biopsy,
  outcome.var="class",
  model.space.priors=list("a"=1, "b"=length(covariate.names), "Variables"=covariate.names), # Beta-binomial(1,P) model space prior
  extra.arguments=list("AlphaPriorMu"=log(mean(biopsy$class)/(1-mean(biopsy$class)))) # Recommend centering intercept prior on logit(event rate)
)
# Simple convergence diagnostic
plot(biopsyResults@mcmc.output[,"LogLikelihood"], type="l") # Looks ok
# Summary plot of selection probabilities
ManhattanPlot(biopsyResults) # Evidence of several independent effects
# Summary table
biopsyResults@posterior.summary.table
TopModels(biopsyResults)

###########################################################
# --- Example 2) Weibull regression survival analysis --- #
###########################################################

utils::data(VA, package = "MASS")
predictors <- c("treat","age","Karn","diag.time","prior")
for (v in predictors) {VA[,v] <- scale(as.numeric(VA[,v]))} # Normalise predictors
VA$stime <- VA$stime/max(VA$stime)# Recommend scaling survival times to between 0 and 1
va.results.weibull <- R2BGLiMS(
  likelihood="Weibull",
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=length(predictors),"Variables"=predictors)) # Beta-binomial(1,P) model space prior
)
plot(va.results.weibull@mcmc.output[,"LogLikelihood"], type="l") # Looks ok
ManhattanPlot(va.results.weibull) # Clear signal at Kern
va.results.weibull@posterior.summary.table
TopModels(va.results.weibull)

######################################################
# --- Example 3) Case-cohort logistic regression --- #
######################################################

# --- Step 1: Logistic RJMCMC
data("CaseCohortExample") # Only V1,V2,V3,V4,V5 have true effects.
for (v in covariate.names) { data.cc[,v] <- data.cc[,v] - mean(data.cc[,v]) }
logistic.cc.res <- R2BGLiMS(
  likelihood="Logistic",
  data=data.cc,
  outcome.var="event",
  model.space.priors=list(list("a"=1,"b"=length(covariate.names),"Variables"=covariate.names)),
  n.mil=0.2,
  extra.arguments=list("AlphaPriorMu"=log(mean(data.cc$event)/(1-mean(data.cc$event))))
)
plot(logistic.cc.res@mcmc.output[,"LogLikelihood"], type="l") # Looks ok
ManhattanPlot(logistic.cc.res) # Identifies 4 out of 5 true effects
logistic.cc.res@posterior.summary.table

# --- Step 2: Effect estimation using the Prentice model
library(survival)
top.vars <- names(which(logistic.cc.res@posterior.summary.table[covariate.names,"BF"]>5))
prentice.model.formula <- as.formula(paste("Surv(times, event) ~ ",paste(top.vars, collapse="+")))
prentice.res <- cch(
  prentice.model.formula,
  data = data.cc,
  subcoh = ~subcohort,
  id=~ID,
  cohort.size=n.complete.cohort, 
  method="Prentice")
summary(prentice.res)

#######################################################
# --- Example 4) Linear regression with full MCMC --- #
#######################################################

data("LinearModelExample") # True effect at V1 only
lm.results <- R2BGLiMS(
  likelihood="Linear", 
  data=data.cts.outcome, 
  outcome.var="y", 
  confounders="confounder", # Example of a variable to include always
  model.space.priors = list("a"=1, "b"=length(covariate.names), "Variables"=covariate.names) # Beta-binomial(1,P) prior on model space
)
plot(lm.results@mcmc.output[,"LogLikelihood"], type="l") # Looks ok
ManhattanPlot(lm.results)
lm.results@posterior.summary.table

#############################################################################
# --- Example 5) Linear regression using a conjugate marginal posterior --- #
#############################################################################

# Coefficients are analytically integrated out
# Only models are sampled so mixing should be better
data("LinearModelExample") # True effect at V1 only
lm.conjugate.results <- R2BGLiMS(
  likelihood="LinearConj", 
  data=data.cts.outcome, 
  outcome.var="y",
  confounders="confounder", # Example of a variable to include always
  tau=max(nrow(data.cts.outcome),ncol(data.cts.outcome)^2), # Tau recommended to be maximum of P^2 and N
  model.space.priors = list("a"=1, "b"=length(covariate.names), "Variables"=covariate.names), # Beta-binomial(1,P)prior on model space
  g.prior=T)
# Summary plot of selection probabilities
ManhattanPlot(lm.conjugate.results) # More decisive than above
TopModels(lm.conjugate.results)
# Posterior sample of parameter in top model
posterior.sample <- NormInvGamPosteriorSample(
  data=data.cts.outcome, outcome.var ="y", confounders="confounder",
  model=c("V1"),tau=100)
# Median and credible interval
round(quantile(posterior.sample[,"V1"],c(0.5, 0.025, 0.975)),2) # Clearly excludes 0

##############################################
# --- Example 6) Fixed priors on effects --- #
##############################################

# By default, a common prior with unkown variance is assumed for all effects
# Might want to explicitly specify these priors instead the priors (e.g. if there is prior information)
data("LinearModelExample") # True effect at V1 only
lm.results.fixed.priors <- R2BGLiMS(
  likelihood="Linear", 
  data=data.cts.outcome, 
  outcome.var="y", 
  confounders="confounder", # Example of a variable to include always
  beta.priors = data.frame(
    cbind(
      rep(0,11), # Normal pior means
      rep(1,11)), # Normal prior SDs
    row.names=c("confounder",covariate.names)),
  model.space.priors = list("a"=1, "b"=length(covariate.names), "Variables"=covariate.names) # Beta-binomial(1,P) prior on model space
)
plot(lm.results.fixed.priors@mcmc.output[,"LogLikelihood"], type="l") # Looks ok
ManhattanPlot(lm.results.fixed.priors)
lm.results.fixed.priors@posterior.summary.table

#################################################
# --- Example 7) Two model space partitions --- #
#################################################

# Can provide different levels of prior sparsity to different groups of covariates
# Demonstrated below with two groups
data("LinearModelExample") # True effect at V1 only
lm.results.two.model.space.partitions <- R2BGLiMS(
  likelihood="Linear", 
  data=data.cts.outcome, 
  outcome.var="y", 
  confounders="confounder", # Example of a variable to include always
  model.space.priors = list(
    list("a"=1, "b"=1, "Variables"=covariate.names[1:5]), # Generous 50/50 prior proportion
    list("a"=1, "b"=1000, "Variables"=covariate.names[6:10]) # Very sparse!
    )
)
plot(lm.results.two.model.space.partitions@mcmc.output[,"LogLikelihood"], type="l") # Looks ok
ManhattanPlot(lm.results.two.model.space.partitions)
lm.results.two.model.space.partitions@posterior.summary.table
