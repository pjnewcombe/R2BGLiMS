library(R2BGLiMS)

### --- Logistic regression with two model space prior partitions
# Example dataset from MASS
utils::data(biopsy, package = "MASS")
# Normalise predictors in order to use (default) common unknown priors on their effect sizes
for (v in paste("V",c(1:9),sep="")) {
  biopsy[,v] <- (biopsy[,v]-mean(biopsy[,v],na.rm=T))/sd(biopsy[,v],na.rm=T)
}
biopsyResults <- R2BGLiMS(
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
predictors <- c("treat","age","Karn","diag.time","cell","prior")
beta.priors <- data.frame(
  cbind(rep(0,length(predictors)),rep(10,length(predictors))),
  row.names=predictors)
# Notice how survival analysis requires additionally flagging a "times.var"
VAResults <- R2BGLiMS(
  likelihood="Weibull",
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=1,"Variables"=predictors)),
  beta.prior=beta.priors
)
PrettyResultsTable(VAResults)
TopModels(VAResults)
