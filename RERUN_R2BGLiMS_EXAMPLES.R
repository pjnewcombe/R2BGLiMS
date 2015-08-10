library(R2BGLiMS)

### --- Logistic regression with two model components
utils::data(biopsy, package = "MASS")
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
# Univariate
uni <- glm(formula(class~V1+V2+V3+V4+V5+V6+V7+V8+V9),biopsy,family=binomial)
summary(uni)
exp(uni$coefficients) # Seems consistent!!!

save(biopsyResults,
     file="~/Dropbox/Work Projects/R Packages/R2BGLiMS/data/biopsyResults.rda")

### --- Survival analysis with priors provided
# Be
utils::data(VA, package = "MASS")
predictors <- c("treat","age","Karn","diag.time","cell","prior")
beta.priors <- data.frame(
  cbind(rep(0,length(predictors)),rep(10,length(predictors))),
  row.names=predictors)
VAResults <- R2BGLiMS(
  likelihood="Weibull",
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=1,"Variables"=predictors)),
  beta.prior=beta.priors
)
save(VAResults,
     file="~/Dropbox/Work Projects/R Packages/R2BGLiMS/data/VAResults.rda")
