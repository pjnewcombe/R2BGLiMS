#' Calculates predicted risks for new data according to a RJMCMC logistic or Weibull regression analysis (for the latter
#' risks are provied for within a user specified time interval)
#' @export
#' @title Calculate predicted risks for new data
#' @name PredictedRisks
#' @inheritParams ResultsTable
#' @inheritParams R2BGLiMS
#' @param t If a Weibull model has been fitted, the length of time to calculate risks over
#' @param use.top.n.pred Generate predictions using only the predictors with the top N marginal posterior probabilities (that appear together). Default NULL to use all predictors
#' 
#' @return A list including a vector of event risks, ordered as the rows of `data'
#' @author Paul Newcombe
PredictedRisks <- function(
  data,
  results,
  times.var=NULL,
  t = NULL,
  use.top.n.pred=NULL
  ) {
  
  pred.risks <- list()
  vars <- colnames(results$results)[colnames(results$results)%in%colnames(data)]
  
  if (!is.null(use.top.n.pred)) {
    res.tab <- ResultsTable(results)
    top.pred <- names(sort(res.tab[results$predictors,"PostProb"],dec=T))[1:use.top.n.pred]
    top.mods.top.pred <- PrettyModelTable(results, which.vars=top.pred)
    model.used <- names(top.mods.top.pred[1,top.mods.top.pred[1,]=="X"])
    # Trim posterior
    n.save.orig <- nrow(results$results)
    top.model.selec <- apply(results$results[,model.used], MAR=1, function(x) (sum(x!=0)==length(model.used))  )
    results$results <- results$results[top.model.selec, ]
    # Log
    vars <- c(results$predictors.fix,model.used)
    pred.risks$top.pred <- top.pred
    pred.risks$model.used <- model.used
    pred.risks$its.used <- n.save.orig/nrow(results$results)
  }
    
  # Calculate predicted survival probabilties, conditional on t
  linpred <- results$results[,"alpha"] + as.matrix(data[,vars])%*%t(results$results[,vars])    
  if (results$args["Likelihood"]=="Weibull") {
    lambda <- exp(linpred)
    k <- exp(results$results[,"LogWeibullScale"])
    risks.its <- 1 - exp( -lambda*(t^k) )
    # I THINK under this parameterisation, lambda is not raised to the k
    # My parameterisation is the same as here:
    # http://www.physicsforums.com/showthread.php?t=649120
    # so the expectation is lambda^(-1/k)*gamma(1+1/k)    
    times.its <- gamma(1+1/k)*(lambda^(-1/k))
  } else if (results$args["Likelihood"]=="Logistic") {
    risks.its <- exp(linpred)/(1+exp(linpred))
    times.its <- NULL
  }
  pred.risks$risks <- apply(risks.its,MAR=1,mean)
  pred.risks$risks.cri <- cbind(
    "mean"=pred.risks$risks,
    "lower"=apply(risks.its,MAR=1,function(x) quantile(x, probs=0.025)),
    "upper"=apply(risks.its,MAR=1,function(x) quantile(x, probs=0.975))
  )
  if (!is.null(times.its)) {
    pred.risks$times <- apply(times.its,MAR=1,mean)
    pred.risks$times.cri <- cbind(
      "mean"=pred.risks$times,
      "lower"=apply(times.its,MAR=1,function(x) quantile(x, probs=0.025)),
      "upper"=apply(times.its,MAR=1,function(x) quantile(x, probs=0.975))
    )    
  }
  
  return(pred.risks)  
}

