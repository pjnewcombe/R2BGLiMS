#' @include ResultsTable.R
NULL

#' Generates a table of the top models from a Reversible Jump Results object.
#' @export
#' @title Table of top models
#' @name TopModels
#' @inheritParams ResultsTable
#' @param vars.to.include A list of variables to restrict to when calculating the model space. Otherwise
#' models are inferred with respect to the complete model space analysed
#' @param n.top.models number of top models to include in the table (default 20)
#' @param post.prob.digits number of digits to round to for posterior percetage column (default 1)
#' @param for.latex Set to true if passing to xtable, and a doublebackslash will be included
#' in front of the percentage sign, and $\\bullet$ will be used as the present symbol
#' @param present.character character to indicate when a variable is present
#' @param post.prob.label Can specify how the post prob should be labelled (if width is a problem)
#' @return Ordered table of best models, with corresponding posterior probabilities.
#' @author Paul Newcombe
#' @example Examples/TopModels_Examples.R 
TopModels <- function(
  results,
  vars.to.include=NULL,
  n.top.models=20,
  post.prob.digits=1,
  for.latex=FALSE,
  present.character="X",
  post.prob.label="Posterior Probability",
  remove.empty.cols=TRUE) {
  
  if (for.latex) {
    per.sign <- "\\%"
    present.character <- "$\\bullet$"
  } else {
    per.sign <- "%"
  }
  s.print.option <- paste("%.",post.prob.digits,"f", sep="")
  
  likelihood.type <- results$args$Likelihood
	results.table <- ResultsTable(results)

  ### Identify model at each iteration
  if(!is.null(vars.to.include)) {
    var.names <- colnames(results$results)
    var.names <- var.names[var.names %in% vars.to.include]
  } else {
    var.indices <- c((results$args$startRJ+2+as.integer(likelihood.type=="Weibull")):(results$args$V+1+as.integer(likelihood.type=="Weibull")) )    
    var.names <- colnames(results$results)[var.indices]
  }
  all.models <- apply(
    results$results[,var.names],
    MAR=1,
    function(r) paste( as.integer(r!=0), collapse="_")
  )
  
  ### --- Make table
  models.table <- sort(table(all.models),d=T)/nrow(results$results)
  n.top.models <- min(n.top.models, length(models.table))
  models.table <- models.table[1:n.top.models]
  models.tab.str <- strsplit(names(models.table), split="_")
  models.tab <- NULL
  for (m in 1:length(models.tab.str)) {
    models.tab <- rbind(models.tab, as.integer(models.tab.str[[m]]) )
  }
  colnames(models.tab) <- var.names
  models.tab <- cbind(
    models.tab,
    "Post Prob"=paste( sprintf(s.print.option, models.table*100), per.sign, sep="")
  )
  rownames(models.tab) <- NULL
  
  ### --- Make prettier
  models.tab[,var.names] <- gsub("0","",x=models.tab[,var.names])
  models.tab[,var.names] <- gsub("1",present.character,x=models.tab[,var.names], fixed=TRUE)
  colnames(models.tab)[ncol(models.tab)] <- post.prob.label
  if (remove.empty.cols) {
    empty.cols <- NULL
    for (v in colnames(models.tab)[colnames(models.tab)!=post.prob.label]) {
      if ( length(grep(present.character, models.tab[,v], fixed=TRUE))==0) {
        empty.cols <- c(empty.cols,which(colnames(models.tab)==v))
      }
    }
    if (!is.null(empty.cols)) {
      models.tab <-models.tab[,-empty.cols]        
    }
  }
    
	return(models.tab)		
}	
