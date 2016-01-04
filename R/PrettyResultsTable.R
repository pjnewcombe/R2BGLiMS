#' Creates a `pretty' summary results table from a Reversible Jump results object
#' @export
#' @title A pretty summary results table for reports. NOTE: This function outputs a table formatted with character strings.
#' A numeric representation of the results are stored in the slot 'posterior.summary.table'.
#' @name PrettyResultsTable
#' @inheritParams ManhattanPlot
#' @param round.digits.betas Number of decimal places to include for effect estimates. (Default is 2)
#' @param round.digits.postprob Number of decimal places to include for posterior probabilities. (Default is 2)
#' @param round.digits.bf Number of decimal places to include for Bayes factors. (Default is 1)
#' @param normalised.sds If covariates were normalised by their standard deviation, provide a named
#' vector of the standard deviations on the original scale using this argument. Effects in the resulting 
#' table will then be interpretable according to unit changes on each covariate's original scale, rather 
#' than according to standard deviation changes.
#' 
#' @return A nice summary results table.
#' @author Paul Newcombe
#' @example Examples/PrettyResultsTable_Examples.R 
PrettyResultsTable <- function(
  results,
  round.digits.betas=2,
  round.digits.postprob=2,
  round.digits.bf=1,
  normalised.sds=NULL
  ) {
  
  # Extract pre-calculated results table
  res.tab <- as.data.frame(results@posterior.summary.table, stringsAsFactors=F)

  # Normalise if provided
  cols.to.normalise <- c("Median", "CrI_Lower", "CrI_Upper", "Median_Present", "CrI_Lower_Present", "CrI_Upper_Present")
  if (!is.null(normalised.sds) ) {
    for (v in intersect(names(normalised.sds), colnames(results@bglims.rjmcmc.output)) ) {
      res.tab[v, cols.to.normalise] <-  res.tab[v, cols.to.normalise]/normalised.sds[v]   # If the original scale was huge, then the effects are much smaller for a unit on the original scale
    }
  }
  
  # Order and rename variables for table
  effect.estimate.cols <- c("Median", "CrI_Lower", "CrI_Upper", "Median_Present", "CrI_Lower_Present", "CrI_Upper_Present")
  res.tab.keep <- c("alpha")
  res.tab.new.names <- c("Intercept")
  if (results@likelihood %in% c("Cox") & c("alpha") %in% row.names(res.tab) ) { # re-log the intercept
    res.tab <- res.tab[-which(row.names(res.tab)=="alpha"),] # Delete interecept row
  }
  if ("LogWeibullScale" %in% rownames(res.tab)) {
    res.tab.keep <- c("LogWeibullScale", res.tab.keep)
    res.tab.new.names <- c("Scale", res.tab.new.names)
    res.tab["LogWeibullScale",effect.estimate.cols] <- exp(res.tab["LogWeibullScale",effect.estimate.cols])
  }
  if ("LogGaussianResidual" %in% rownames(res.tab)) {
    res.tab.keep <- c("LogGaussianResidual", res.tab.keep)
    res.tab["LogGaussianResidual",] <- exp(res.tab["LogGaussianResidual",])
    res.tab.new.names <- c("Residual", res.tab.new.names)
  }
  if (results@bglims.arguments$nBetaHyperPriorComp>0) {
    res.tab.keep <- c(res.tab.keep, paste("LogBetaPriorSd",c(1:results@bglims.arguments$nBetaHyperPriorComp), sep="") )
    if (results@bglims.arguments$nBetaHyperPriorComp==1) {
      res.tab.new.names <- c(res.tab.new.names, "log(beta) Hyperprior SD")
    } else {
      res.tab.new.names <- c(res.tab.new.names, paste("log(beta) Hyperprior SD - component",c(1:results@bglims.arguments$nBetaHyperPriorComp)) )
    }
  }
  predictors <- rownames(res.tab[!rownames(res.tab)%in%c(
    "LogWeibullScale", "LogGaussianResidual", "alpha",
    paste("LogBetaPriorSd",c(1:results@bglims.arguments$nBetaHyperPriorComp), sep=""),
    "LogLikelihood"),])
  res.tab.keep <- c(res.tab.keep,predictors)
  res.tab.new.names <- c(res.tab.new.names,predictors)
  
  # Re-order/rename variables
  res.tab <- res.tab[res.tab.keep,]
  rownames(res.tab) <- res.tab.new.names
  
  # Round and replace with >0.01
  beta.cols <- c("Median", "CrI_Lower", "CrI_Upper", "Median_Present", "CrI_Lower_Present", "CrI_Upper_Present")
  res.tab[,beta.cols] <- round(res.tab[,beta.cols], round.digits.betas)
  pretty.tab <- as.data.frame(cbind(
    "median"=paste(format(res.tab[,"Median"], trim=TRUE, sci=F)),
    "95% CrI"= paste("(", format(res.tab[,"CrI_Lower"], trim=TRUE, sci=F), ", ",
      format(res.tab[,"CrI_Upper"], trim=TRUE, sci=F), ")", sep=""),
    "median*"=paste(format(res.tab[,"Median_Present"], trim=TRUE, sci=F)),
    "95% CrI*"= paste("(", format(res.tab[,"CrI_Lower_Present"], trim=TRUE, sci=F), ", ",
                     format(res.tab[,"CrI_Upper_Present"], trim=TRUE, sci=F), ")", sep=""),
    "Posterior Probability"=paste(format(round(res.tab[,"PostProb"], round.digits.postprob), sci=F) ),
    "Bayes Factor"=paste(format(round(res.tab[,"BF"], round.digits.bf), sci=F) )
    ))
  
  # Replace those that come out at 0
  beta.col.names <- c("median", "95% CrI", "median*", "95% CrI*")
  to.replace <- paste("0.",paste(rep(0,round.digits.betas), collapse=""), sep="")
  replace.with <- paste("<0.",paste(rep(0,round.digits.betas-1), collapse=""), "1", sep="")
  for (v in beta.col.names) {
    pretty.tab[,v] <- gsub(to.replace, replace.with, pretty.tab[,v])
  }
  to.replace <- paste("0.",paste(rep(0,round.digits.postprob), collapse=""), sep="")
  replace.with <- paste("<0.",paste(rep(0,round.digits.postprob-1), collapse=""), "1", sep="")
  pretty.tab[,"Posterior Probability"] <- gsub(to.replace, replace.with, pretty.tab[,"Posterior Probability"])
  to.replace <- paste("0.",paste(rep(0,round.digits.bf), collapse=""), sep="")
  replace.with <- paste("<0.",paste(rep(0,round.digits.bf-1), collapse=""), "1", sep="")
  pretty.tab[,"Bayes Factor"] <- gsub(to.replace, replace.with, pretty.tab[,"Bayes Factor"])
  
  # Measure to include in column titles
  measure <- "Effect"
  if (results@likelihood %in% c("Logistic", "JAM", "JAM_MCMC")) {
    measure <- "OR"
  } else if (results@likelihood %in% c("Cox", "Weibull")) {
    measure <- "HR"
  }
  colnames(pretty.tab)[colnames(pretty.tab) %in% beta.col.names] <- paste(measure, beta.col.names)    

  # Insert NAs for variables included at all times
  na.string <- paste(paste(rep(" ",round.digits.postprob),collapse=""), "NA", sep="")
  pretty.tab[which(pretty.tab[,"Posterior Probability"]==na.string), c(1,2) ] <- na.string
  
  # Add names
  rownames(pretty.tab) <- rownames(res.tab)
  
  # Only keep posterior probablities and Bayes Factors for conjugate models
  if (results@likelihood %in% c("JAM", "GaussianConj")) {
    pretty.tab <- pretty.tab[predictors, c("Posterior Probability", "Bayes Factor")]
  }

  return(pretty.tab)  
}	
