library(R2MHRJ)

### --- Logistic regression with two model components
utils::data(biopsy, package = "MASS")
beta.priors <- data.frame(
  cbind(rep(0,9),rep(10,9)),
  row.names=paste("V",c(1:9),sep=""))
logistic.res <- R2MHRJ(
  data=biopsy, outcome.var="class",
  model.space.priors=list(
    list("Rate"=0.5, "Variables"=paste("V",c(1:3),sep="")),
    list("Rate"=0.1, "Variables"=paste("V",c(4:9),sep=""))
  )
)
AutocorrelationPlot(logistic.res)
ManhattanPlot(logistic.res)
PrettyResultsTable(logistic.res)
PrettyModelTable(logistic.res)

### --- Survival analysis with priors provided
utils::data(VA, package = "MASS")
predictors <- c("treat","age","Karn","diag.time","cell","prior")
beta.priors <- data.frame(
  cbind(rep(0,length(predictors)),rep(10,length(predictors))),
  row.names=predictors)
surv.res <- R2MHRJ(
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=1,"Variables"=predictors)),
  beta.prior=beta.priors
)

# Two
surv.res <- R2MHRJ(
  data=VA, outcome.var="status", times.var="stime",
  model.space.priors=list(
    list("Rate"=0.1,"Variables"=predictors[1:3]),
    list("Rate"=0.5,"Variables"=predictors[4:6])
  ),
  beta.prior=beta.priors
)

# Results
ManhattanPlot(surv.res)
PrettyResultsTable(surv.res)
PrettyModelTable(surv.res)



#attach(nki70)

### --- Setup variables
confounders <- "year"
predictors <- c("sex"       "age"       "year"      "thickness" "ulcer")
preds1 <- colnames(nki70)[8:50]
preds2 <- colnames(nki70)[51:77]

### --- Setup fixed priors (these variables probably aren't exchangeable)
n.var <- length(preds1)+length(preds2)+1
beta.priors <- cbind(rep(0, n.var), rep(1, n.var))
rownames(beta.priors) <- c(confounders, preds1, preds2)

################
### --- EXAMPLE: Survival, beta-binomial model prior with one component, fixed beta priors
################

surv.res <- mhrj(
  data=nki70,
  outcome.var="event",
  times.var="time",
  confounders=confounders,
  model.space.priors=list(
    list("a"=1, "b"=1, "Variables"=c(preds1,preds2) )
  ),
  beta.priors=beta.priors,
  n.mil=1
)

### --- Convergence checks
AutocorrelationPlot(surv.res)
# ChainPlots(surv.res) # The pdf option should be used, since multiple pages are required

### --- Predictor Results. Notice how Age, specified as a confounder, has NA for post prob and BF
par(mfrow=c(1,1)) # plots do not tidy up after themselves
ManhattanPlot(surv.res)
ResultsTable(surv.res)
PrettyResultsTable(surv.res)

### --- Model Results
PrettyModelTable(surv.res) # List of top models + posterior probabilities
ModelSpaceSummary(surv.res) # Model size

names(surv.res)
surv.res$run.time # Run time
# surv.res$results contains all saved MCMC output. Each row is an iteration, and a 0 value
# indicates exclusion from the model at that iteration.
# This can, e.g., be used to plot the posterior distribution of model size
hist(
  apply(surv.res$results[,c(preds1,preds2)], MAR=1, function(x) sum(x!=0)),
  prob=T, breaks=20,
  xlab="No. included predictors", main="Posterior model size"
  )

################
### --- EXAMPLE: Logistic, Poisson model priors in two components, shared beta priors with unknown SD
################

# For logistic regression, simply do not specify times.var.
# Two have two model space prior components, model.space.priors is now a list of two lists.
# There are three options regarding the specification of priors for the beta coefficients:
# 1 - do nothing. N(0, 1e6) priors are assigned to the betas of any confounders, and predictors included
# in model selection are assigned a common prior with unknown variance (which in turn has a InverseGamma
# hyperprior)
# 2 - provide beta.priors including confounders only. This overrides the default use of N(0,1e6) confounder
# priors, but the predictors for which model selection is performed for are still assigned a common prior
# with unknown variance
# 3 - provide beta.priors including all confounders and predctors. This overrides the use of a common
# prior with unkown variance across predictors for which model selection is performed for.
# In the example below we provide an informative prior for the confounder only.

beta.prior.age <- cbind(log(0.9), 0.00001) # Ridoculously informative prior on age log-OR
rownames(beta.prior.age) <- "Age"
logistic.res <- mhrj(
  data=nki70,
  outcome.var="event",
  confounders=confounders,
  model.space.priors=list(
    list("Rate"=0.1, "Variables"=preds1),
    list("Rate"=0.7, "Variables"=preds2)
  ),
  beta.priors=beta.prior.age, # Specified for confounders only (not specifu)
  n.mil=1
)
# NOTE: Results are meaningless since this is survival data (and variables are not normalised so a common
# prior does not make sense)
# However, you can see the effect of a two part model space prior, with markedly different rates
AutocorrelationPlot(logistic.res)
