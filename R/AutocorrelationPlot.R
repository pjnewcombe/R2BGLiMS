#' Generates an autocorrelation plot from a Java MCMC results class object
#' 
#' @export
#' @title Autocorrelation plot for model selection. Model inclusion of each predictor is
#' represented across all iterations as a black and white heatmap.
#' @name AutocorrelationPlot
#' @inheritParams ManhattanPlot
#' @param plot.title Optionally specify a title for the plot.
#' @param file Path to a PDF to print plot to.
#' @param include.var.names.on.x.axis Include variable names perpendicularly below the x-axis. Turn off if it looks too crowded.
#' @param include.post.probs Include a barplot of posterior probabilities at the top.
#' @param cex.x.axis.ticks Control the size of tick marks for the covariates along the X-axis. Decrease if there are many.
#' @param cex.y.axis.ticks Control the size of tick marks for the Y-axis.
#' @param cex.y.axis.labels Control the size of the y-labels. Decrease if these are not fitting on.
#' @param mar.for.varnames Margin size for the bottom axis. Default is 8 so there is room for perpendicular variable names. 
#' Decrease if not using variable names.
#' @return NA
#' @author Paul Newcombe
#' @example Examples/AutocorrelationPlot_Examples.R
AutocorrelationPlot <- function(
  results,
  plot.title=NULL,
  file=NULL,
  include.var.names.on.x.axis=TRUE,
  covariates.to.include=NULL,
  var.dictionary=NULL,
  include.post.probs=TRUE,
  cex.x.axis.ticks=1,
  cex.y.axis.ticks=1,
  cex.y.axis.labels=1,
  mar.for.varnames=8) {
    
  # Save initial graphical parameters
  original.par <- par(no.readonly = TRUE)
  
  # Setup for plot
  if (is.null(covariates.to.include)) {
    # If variable lsit not provided, take from the model space prior
    covariates.to.include <- unlist(lapply(results@model.space.priors, function(x) x$Variables))
  }
  results@mcmc.output <- results@mcmc.output[,covariates.to.include]
  n.var <- ncol(results@mcmc.output)
  x <- c(1:n.var)
  y <- seq(
    from=(ceiling(results@burnin.fraction*results@bglims.arguments$iterations)+results@bglims.arguments$thin),
    to=results@bglims.arguments$iterations,
    by=results@bglims.arguments$thin)
  z <- results@mcmc.output!=0
  z <- t(z)+0 # Makes numeric
  
  # Initiate Plot
  if(!is.null(file)) {
    pdf(file=file,
        paper="special",
        width=8,
        height=11)
  }
  if (include.post.probs) {
    par(mfrow=c(6,1))
    layout(matrix(c(1,2,2,2,2,2), 6, 1, byrow = TRUE))
    par(mar=c(0, 4*cex.y.axis.labels, 4, 2) + 0.1)
    
    # 1 - Barplot  
    barplot(
      height=apply(results@mcmc.output,MAR=2,function(x) sum(x!=0)/length(x) ),
      ylim=c(0,1),
      main=plot.title,
      axes=FALSE,
      ylab="Post Prob",
      xlab="",
      names.arg=rep("",ncol(results@mcmc.output)),
      space=0,
      xlim=c(0,ncol(results@mcmc.output)),
      xaxs="i",
      col="black", cex.lab=cex.y.axis.labels)
    axis(2,
         at=seq(from=0,to=1,by=0.25),
         labels=c("","0.25","0.5","0.75","1"), cex.axis=cex.y.axis.ticks)
    mar.top <- 0.1
  } else {
    mar.top <- 4
  }
  
  # 2 - Autocorrelation plot
  if (include.var.names.on.x.axis) {
    par(mar=c(mar.for.varnames, 4*cex.y.axis.labels, mar.top, 2) + 0.1)
    image(x,y,z, axes=F, xlab="", ylab="Iteration", col=c("white", "black"), cex.lab=cex.y.axis.labels)
    tick.labels <- colnames(results@mcmc.output)
    if (!is.null(var.dictionary)) {
      tick.labels[tick.labels%in%names(var.dictionary)] <- var.dictionary[tick.labels[tick.labels%in%names(var.dictionary)]]
    }
    axis(side=1, at=c(1:ncol(results@mcmc.output)), labels=tick.labels, las=2, cex.axis=cex.x.axis.ticks)
  } else {
    par(mar=c(5, 4*cex.y.axis.labels, mar.top, 2) + 0.1)
    image(x,y,z, axes=F, xlab="Predictors", ylab="Iteration", col=c("white", "black"), cex.lab=cex.y.axis.labels)
    axis(side=1,at=c(0.5, (length(x)+0.5)),labels=rep("",2), cex.axis=cex.x.axis.ticks)
  }
  y.ticks <- seq(
    from=ceiling(results@burnin.fraction*results@bglims.arguments$iterations),
    to=results@bglims.arguments$iterations,
    by=0.5e6)
  y.tick.labs <- paste(y.ticks/1e6,"m",sep="")
  y.tick.labs[grep(".5m",y.tick.labs)] <- ""
  axis(side=2, at=y.ticks, labels=y.tick.labs, col.axis="black", cex.axis=cex.y.axis.ticks)
  
  # Close pdf
  if (!is.null(file)) {
    dev.off()    
  }
  
}	
