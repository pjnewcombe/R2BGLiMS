library(R2BGLiMS)

### --- Linear regression with full MCMC of all parameters
set.seed(1)
# Simulate data
X <- matrix(rnorm(1000, 0, 1), 100,10) # n=100, 10 covariates
confounder <- rnorm(100, 0, 1) 
b <- c(1, rep(0,9), 4) # Only the first covariate has an effect
y <- rnorm(100, cbind(X,confounder)%*%b,1)
data <- data.frame(cbind(X,confounder,y))
# Run model treating prior variance of the effects as unknown
# Poisson model space prior on number of effects
poisson.model.space.prior=list( list("Rate"=0.1, "Variables"=paste("V",c(1:10),sep="")) ) 
lm.hierarchical.prior.results <- R2BGLiMS(
  likelihood="Gaussian", data=data, outcome.var="y", confounders="confounder",
  model.space.priors = poisson.model.space.prior
)
PrettyResultsTable(lm.hierarchical.prior.results)
# A numeric representation of the results is stored in the slot 'posterior.summary.table':
lm.hierarchical.prior.results@posterior.summary.table
# Re-analyse using fixed priors
# (not enough covariates to estimate effect variance)
beta.priors <- data.frame(
  cbind(rep(0,11),rep(10,11)),
  row.names=c("confounder",paste("V",c(1:10),sep="")) )
lm.fixed.prior.results <- R2BGLiMS(
  likelihood="Gaussian", data=data, outcome.var="y", confounders="confounder",
  model.space.priors = poisson.model.space.prior,
  beta.priors=beta.priors)
PrettyResultsTable(lm.fixed.prior.results) # Mixing seems a little poorer
# Re-analyse using the conjugate Gaussian model
# Only models are sampled: effects are analytically integrated out
lm.conjugate.results <- R2BGLiMS(
  likelihood="GaussianConj", data=data, outcome.var="y",confounders="confounder",
  model.space.priors=poisson.model.space.prior,
  enumerate.up.to.dim=0,
  tau=100, # Tau recommended to be maximum of P^2 and N
  g.prior=T)
PrettyResultsTable(lm.conjugate.results)
TopModels(lm.conjugate.results)
# Obtain posterior sample of parameter V1
posterior.sample <- NormInvGamPosteriorSample(
  data=data, outcome.var ="y", confounders="confounder",
  model=c("V1"),tau=100)
# Median and credible interval
round(quantile(posterior.sample[,"V1"],c(0.5, 0.025, 0.975)),2)
# Clearly excludes 0

### --- Logistic regression with two model space prior partitions
# Example dataset from MASS
utils::data(biopsy, package = "MASS")
# Normalise predictors in order to use (default) common unknown priors on their effect sizes
for (v in paste("V",c(1:9),sep="")) {
  biopsy[,v] <- (biopsy[,v]-mean(biopsy[,v],na.rm=T))/sd(biopsy[,v],na.rm=T)
}
biopsyResults <- R2BGLiMS( # Takes a few minutes to run
  likelihood="Logistic", data=biopsy, outcome.var="class",
  model.space.priors=list(
    list("Rate"=0.5, "Variables"=paste("V",c(1:3),sep="")),
    list("Rate"=0.1, "Variables"=paste("V",c(4:9),sep=""))
  )
)
PrettyResultsTable(biopsyResults)
TopModels(biopsyResults)

### --- Survival analysis with priors specified on effect sizes
# Example dataset from MASS
utils::data(VA, package = "MASS")
# Generate priors on log-effect sizes
predictors <- c("treat","age","Karn","diag.time","prior")
beta.priors <- data.frame(
  cbind(rep(0,length(predictors)),rep(10,length(predictors))),
  row.names=predictors)
# Analysis using the Weibull likelihood
# Notice how survival analysis requires additionally flagging a "times.var"
va.results.weibull <- R2BGLiMS(
  likelihood="Weibull",
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=1,"Variables"=predictors)),
  beta.prior=beta.priors
)
PrettyResultsTable(va.results.weibull)
TopModels(va.results.weibull)
# Analysis using the Cox likelihood
va.results.cox <- R2BGLiMS(
  likelihood="Cox",
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=1,"Variables"=predictors)),
  beta.prior=beta.priors
)
PrettyResultsTable(va.results.cox)
TopModels(va.results.cox)
