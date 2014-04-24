#' Generates a PDF of chain plots for all variables from a Reversible Jump results object
#' @export
#' @title Coefficient history trace plots from a Reversible Jump results object
#' @name ChainPlots
#' @inheritParams ResultsTable
#' @param predictors.only If true only chain plots for predictors are plotted (and intercept,
#'  likelihood etc. are excluded). Default is FALSE.
#' @param plot.file PDF path
#' @param add.to.title Optional character string to add to the start of each plot title
#' @param add.panel.labels Add lower case letter labels to panel titles. Default FALSE.
#' Note: Does not work when providing a var.dictionary
#' @param par.mfrow Allow passing the mfrow layout vector to par. Default par(mfrow=c(5,2))
#' @return NA
#' @author Paul Newcombe
#' @example Examples/ChainPlots_Examples.R 
ChainPlots <- function(
  results,
  vars.to.include=NULL,
  var.dictionary=NULL,
  predictors.only=FALSE,
  plot.file=NULL,
  add.to.title=NULL,
  par.mfrow=NULL,
  add.panel.labels=FALSE) {

  if (!is.null(var.dictionary)) {
    add.panel.labels <- FALSE
  }
  
  if(!is.null(plot.file)) {
    pdf(file=plot.file,
        paper="special",
        width=8,
        height=11)
  }

  # Rename non-variables
  cols.keep <- c("LogLikelihood", "alpha")
  cols.new.names <- c("Log-likelihood", "Intercept")
  if ("LogWeibullScale" %in% colnames(results$results)) {
    cols.keep <- c(cols.keep, "LogWeibullScale")
    cols.new.names <- c(cols.new.names, "Scale")
  }
  if (results$args$nBetaHyperPriorComp>0) {
    cols.keep <- c(cols.keep, paste("LogBetaPriorSd",c(1:results$args$nBetaHyperPriorComp), sep="") )
    if (results$args$nBetaHyperPriorComp==1) {
      cols.new.names <- c(cols.new.names, "log(beta) Hyperprior SD")
    } else {
      cols.new.names <- c(cols.new.names, paste("log(beta) Hyperprior SD - component",c(1:results$args$nBetaHyperPriorComp)) )
    }
  }
  if(!is.null(vars.to.include)) {
    predictors <- vars.to.include
  } else {
    predictors <- colnames(results$results)[!colnames(results$results)%in%c(
      "LogWeibullScale", "alpha",
      paste("LogBetaPriorSd",c(1:results$args$nBetaHyperPriorComp), sep=""),
      "LogLikelihood")]
  }
  cols.keep <- c(cols.keep,predictors)
  cols.new.names <- c(cols.new.names,predictors)
  if (predictors.only) {
    cols.keep <- predictors
    cols.new.names <- predictors
  }
  results$results <- results$results[,cols.keep]
  colnames(results$results) <- cols.new.names
  
  # Plot
  if (!is.null(par.mfrow)) {
    par(mfrow=par.mfrow)
  } else {
    par(mfrow=c(5,2))
  }
  letter.ind <- 0
  for (v in colnames(results$results)) {
    letter.ind <- letter.ind+1
    main.title <- add.to.title
    if (add.panel.labels) {
      main.title <- paste(letters[letter.ind],") ",main.title,sep="")
    }
    main.title <- paste(main.title, v, sep="")    
    if (sum(!is.na(results$results[,v]))==0) {
      main.title <- paste(main.title,"*",sep="")
    }
    results$results[,v][results$results[,v]==0] <- NA
    if (!is.null(var.dictionary)) {
      # Overrides everything previously
      if (v %in% names(var.dictionary)) {
        main.title <- var.dictionary[v]
      }
    }
    if (sum(!is.na(results$results[,v]))>0) {
      plot(results$results[,v],
           type="l",
           ylab="Value",
           xlab="Saved Iteration",
           main=main.title )				
    } else {
      plot(rep(0,nrow(results$results)),
           type="n",
           ylab="Value",
           xlab="Saved Iteration",
           main=main.title )
    }    
  }
  	
	if(!is.null(plot.file)) {
	  dev.off()
	}
	
}	
