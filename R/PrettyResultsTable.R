#' @include ResultsTable.R
NULL

#' Creates a `pretty' summary results table from a Reversible Jump results object
#' @export
#' @title A pretty summary results table for reports - note that entries are formatted as character strings so,
#' for processing results use \code{\link{ResultsTable}}.
#' @name PrettyResultsTable
#' @inheritParams ResultsTable
#' @param round.digits.betas Number of decimal places to include for effect estimates. (Default is 2)
#' @param round.digits.postprob Number of decimal places to include for posterior probabilities. (Default is 2)
#' @param round.digits.bf Number of decimal places to include for Bayes factors. (Default is 1)
#' @param measure Abbreviated measure to include in column titles, e.g. "HR". (Default is "OR")
#' @param post.prob.cut Optionally can provide a  posterior probability threshold for including
#' predictors in the table for which model selection has been performed for.
#' @param abbreviated.names Use abbreviated names for Posterior Probability and Bayes Factors columns? (Default is FALSE)
#' @param remove.col.white.spaces Should the column titles have no white spaces? (Default is FALSE)
#' 
#' @return Prints a nice summary results table 
#' @author Paul Newcombe
#' @example Examples/PrettyResultsTable_Examples.R 
PrettyResultsTable <- function(
  results,
  round.digits.betas=2,
  round.digits.postprob=2,
  round.digits.bf=1,
  measure="OR",
  vars.to.include=NULL,
  var.dictionary=NULL,
  normalised.sds=NULL,
  post.prob.cut=NULL,
  abbreviated.names=FALSE,
  remove.col.white.spaces=FALSE
  ) {
  # Normalise if provided
  if ( !is.null(normalised.sds) ) {
    for (v in intersect(names(normalised.sds), colnames(results$results)) ) {
      results$results[,v] <- results$results[,v]/normalised.sds[v]   # If the original scale was huge, then the effects are much smaller for a unit on the original scale
    }
  }
  
  # Establish whether a results table, or a results object has been provided
  if (is.list(results)) {
    res.tab <- as.data.frame(ResultsTable(results, var.dictionary=var.dictionary), stringsAsFactors=F)    
  } else {
    res.tab <- results
  }
  
  # Order and rename variables for table
  res.tab.keep <- c("alpha")
  res.tab.new.names <- c("Intercept")
  if (results$args$Likelihood %in% c("Weibull", "Logistic") ) { # re-log the intercept
    res.tab["alpha",] <- log(res.tab["alpha",])
  }
  if ("LogWeibullScale" %in% rownames(res.tab)) {
    res.tab.keep <- c("LogWeibullScale", res.tab.keep)
    res.tab.new.names <- c("Scale", res.tab.new.names)
  }
  if (results$args$nBetaHyperPriorComp>0) {
    res.tab.keep <- c(res.tab.keep, paste("LogBetaPriorSd",c(1:results$args$nBetaHyperPriorComp), sep="") )
    if (results$args$nBetaHyperPriorComp==1) {
      res.tab.new.names <- c(res.tab.new.names, "log(beta) Hyperprior SD")
    } else {
      res.tab.new.names <- c(res.tab.new.names, paste("log(beta) Hyperprior SD - component",c(1:results$args$nBetaHyperPriorComp)) )
    }
  }
  if(!is.null(vars.to.include)) {
    predictors <- vars.to.include
  } else {
    predictors <- rownames(res.tab[!rownames(res.tab)%in%c(
      "LogWeibullScale", "alpha",
      paste("LogBetaPriorSd",c(1:results$args$nBetaHyperPriorComp), sep=""),
      "LogLikelihood"),])
  }
  res.tab.keep <- c(res.tab.keep,predictors)
  res.tab.new.names <- c(res.tab.new.names,predictors)
  
  # Re-order/rename variables
  res.tab <- res.tab[res.tab.keep,]
  rownames(res.tab) <- res.tab.new.names
  
  # Filter on posterior probability
  if (!is.null(post.prob.cut)) {
    res.tab <- res.tab[
      (is.na(res.tab[,"PostProb"]) | res.tab[,"PostProb"]>=post.prob.cut),]
  }
  
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
  if (!is.null(measure)) {
    colnames(pretty.tab)[colnames(pretty.tab) %in% beta.col.names] <- paste(measure, beta.col.names)    
  }
  
  # Abbreviate names
  if (abbreviated.names == TRUE) {
    colnames(pretty.tab)[c(5:6)] <- c("Post Prob", "BF")
  }
  
  # Remove white spaces from column titles
  if (remove.col.white.spaces==TRUE) {
    colnames(pretty.tab) <- gsub(" ", "_", colnames(pretty.tab)) 
  }
  
  # Insert NAs for variables included at all times
  na.string <- paste(paste(rep(" ",round.digits.postprob),collapse=""), "NA", sep="")
  pretty.tab[which(pretty.tab[,"Posterior Probability"]==na.string), c(1,2) ] <- na.string
  
  # Add names
  rownames(pretty.tab) <- rownames(res.tab)

  return(pretty.tab)  
}	
