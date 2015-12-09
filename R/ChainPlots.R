#' Generates a PDF of chain plots for all variables from a Reversible Jump results object
#' @export
#' @title Coefficient history trace plots from a Reversible Jump results object
#' @name ChainPlots
#' @inheritParams ManhattanPlot
#' @param predictors.only If true only chain plots for predictors are plotted (and intercept,
#'  likelihood etc. are excluded). Default is FALSE.
#' @param file PDF path
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
  file=NULL,
  add.to.title=NULL,
  par.mfrow=NULL,
  add.panel.labels=FALSE) {

  if (!is.null(var.dictionary)) {
    add.panel.labels <- FALSE
  }
  
  if(!is.null(file)) {
    pdf(file=file,
        paper="special",
        width=8,
        height=11)
  }

  ### --- Rename key model parameters (likelihood etc) and add dimension
  results@bglims.rjmcmc.output <- cbind("ModelDimension"=apply(results@bglims.rjmcmc.output, MAR=1, function(x) sum(x!=0) ), results@bglims.rjmcmc.output)
  cols.keep <- c("ModelDimension", "LogLikelihood", "alpha")
  cols.new.names <- c("Model Dimension", "Log-likelihood", "Intercept")
  if ("LogWeibullScale" %in% colnames(results@bglims.rjmcmc.output)) {
    cols.keep <- c(cols.keep, "LogWeibullScale")
    cols.new.names <- c(cols.new.names, "Scale")
  } else if ("LogGaussianResidual" %in% colnames(results@bglims.rjmcmc.output)) {
    cols.keep <- c(cols.keep, "LogGaussianResidual")
    cols.new.names <- c(cols.new.names, "Residual")
  }
  if (results@bglims.arguments$nBetaHyperPriorComp>0) {
    cols.keep <- c(cols.keep, paste("LogBetaPriorSd",c(1:results@bglims.arguments$nBetaHyperPriorComp), sep="") )
    if (results@bglims.arguments$nBetaHyperPriorComp==1) {
      cols.new.names <- c(cols.new.names, "log(beta) Hyperprior SD")
    } else {
      cols.new.names <- c(cols.new.names, paste("log(beta) Hyperprior SD - component",c(1:results@bglims.arguments$nBetaHyperPriorComp)) )
    }
  }
  if(!is.null(vars.to.include)) {
    predictors <- vars.to.include
  } else {
    predictors <- colnames(results@bglims.rjmcmc.output)[!colnames(results@bglims.rjmcmc.output)%in%c(
      "LogWeibullScale", "LogGaussianResidual", "alpha",
      paste("LogBetaPriorSd",c(1:results@bglims.arguments$nBetaHyperPriorComp), sep=""),
      "LogLikelihood")]
  }
  cols.keep <- c(cols.keep,predictors)
  cols.new.names <- c(cols.new.names,predictors)
  if (predictors.only) {
    cols.keep <- predictors
    cols.new.names <- predictors
  }
  results@bglims.rjmcmc.output <- results@bglims.rjmcmc.output[,cols.keep]
  colnames(results@bglims.rjmcmc.output) <- cols.new.names
  
  ### --- Plot
  if (!is.null(par.mfrow)) {
    par(mfrow=par.mfrow)
  } else {
    par(mfrow=c(5,2))
  }
  letter.ind <- 0
  for (v in colnames(results@bglims.rjmcmc.output)) {
    letter.ind <- letter.ind+1
    main.title <- add.to.title
    if (add.panel.labels) {
      main.title <- paste(letters[letter.ind],") ",main.title,sep="")
    }
    main.title <- paste(main.title, v, sep="")    
    if (sum(!is.na(results@bglims.rjmcmc.output[,v]))==0) {
      main.title <- paste(main.title,"*",sep="")
    }
    results@bglims.rjmcmc.output[,v][results@bglims.rjmcmc.output[,v]==0] <- NA
    if (!is.null(var.dictionary)) {
      # Overrides everything previously
      if (v %in% names(var.dictionary)) {
        main.title <- var.dictionary[v]
      }
    }
    if (sum(!is.na(results@bglims.rjmcmc.output[,v]))>0) {
      plot(results@bglims.rjmcmc.output[,v],
           type="l",
           ylab="Value",
           xlab="Saved Iteration",
           main=main.title )				
    } else {
      plot(rep(0,nrow(results@bglims.rjmcmc.output)),
           type="n",
           ylab="Value",
           xlab="Saved Iteration",
           main=main.title )
    }    
  }
  	
	if(!is.null(file)) {
	  dev.off()
	}
	
}	
