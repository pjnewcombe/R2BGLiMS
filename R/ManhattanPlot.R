#' Produces a Manhattan plot of marginal predictor evidence from a Reversible Jump Results object.
#' @export
#' @title Manhattan Plot
#' @name ManhattanPlot
#' @param results A \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} object frum running the \code{\link[R2BGLiMS]{R2BGLiMS}} command.
#' @param plot.quantity Can be "PosteriorProbability" (default), "BayesFactor", "ABF", or "pvalue". The
#' latter is implemented to allow easy comparison between analysis frameworks using the same plot style,
#' and only works if a vector of p-values is supplied. ABF is for plotting Approximate Bayes Factors
#' in favour of the null (i.e. small values correspond to evidence of association), as proposed in 
#' Wakefield, J. (2009). Genetic Epidemiology.
#' @param results.vector Alternatively can simply pass a named vector of probabilties
#' @param plot.title Optional character string to name the plot with
#' @param y.max If plotting Bayes Factors, this is an optional upper limit for the y-axis (e.g. 
#' to make comparison with a different set of results easier on the eye). Note that if p-values are
#' being plotted, then this can be used to supply the maximum -log10(p-value) on the y-axis.
#' @param cex.axis Optional size setting for the variable names relative to the x-axis. May want to
#' adjust if variable names are not legible. (default is 0.5)
#' @param add.bonferroni If plotting p-values this option draws a dashed red line at the Bonferroni
#' threshold. (default is TRUE).
#' @param top.hits If a list of "top hits" are provided here, along with the original X matrix below,
#' then the points are coloured according to their pairwise correlation with the "top hits" (default NULL)
#' @param X.mat Must be provided if top.hits are provided above - used to calculate pairwise correlations.
#' @param point.cols Optional vector of colours to use for the points (as many and in the same
#' order as the covariates)
#' @param point.labs Optional vector of labels to use for the points (as many and in the same
#' order as the covariates)
#' @param supress.x.ticks Set to TRUE to avoid adding x tick mark labels (e.g. if there are many many
#' variables)
#' @author Paul Newcombe
#' @example Examples/ManhattanPlot_Examples.R 
ManhattanPlot <- function(
  results=NULL,
  plot.quantity="PosteriorProbability",
  results.vector=NULL,
  vars.to.include=NULL,
  var.dictionary=NULL,
  plot.title="Manhattan Plot",
  y.max=NULL,
  add.bonferroni=FALSE,
  top.hits=NULL,
  X.mat=NULL,
  point.cols=NULL,
  point.labs=NULL,
  suppress.x.ticks=FALSE
  ) {
  ### --- Errors
  if (!is.null(top.hits)) {
    if(is.null(X.mat)) stop("Must provide original X matrix to indictae correlation with top hits")
  }
  
  ### --- Setup PROVIDED results.vector
  if (!is.null(results.vector)){
    if(!is.null(vars.to.include)) {
      results.vector <- results.vector[vars.to.include]
    } else {
      vars.to.include <- names(results.vector)
    }
  }
  
  ### --- Pre-process results object -> results.vector
  if (!is.null(results)) {
    # By default exclude these boring modelling parameters
    if (is.null(vars.to.include)) {
      vars.to.include <- unlist(lapply(results@model.space.priors, function(x) x$Variables))
    }
    results.table <- results@posterior.summary.table
    results.table <- results.table[vars.to.include,]
    if (plot.quantity=="PosteriorProbability") {
      results.vector <- results.table[,"PostProb"]      
    } else if (plot.quantity=="BayesFactor") {
      results.vector <- log10(results.table[,"BF"])
    }
    names(results.vector) <- rownames(results.table)
  }
  
  ### --- Correlation representative point colours (if x-matrix provided)
  if (is.null(point.cols)) {
    point.cols <- "black"    
  }
  if (!is.null(top.hits)) {
    library(gplots)
    rsqmat <- cor(X.mat[,vars.to.include], use="pairwise.complete.obs")
    cor.with.top.hit <- NULL
    for (v in vars.to.include) {
      cor.with.top.hit <- c(cor.with.top.hit, max(abs(rsqmat[v,top.hits])))
    }
    names(cor.with.top.hit) <- vars.to.include
    palette( rev(rich.colors(32)) ) # colors: 1 to 32
    point.cols <- 1+31*(1-cor.with.top.hit)
  }
  
  ### --- Setup x tick labels
  if (suppress.x.ticks) {
    x.ticks <- FALSE
    x.tick.labels <- FALSE    
  } else {
    x.ticks <- TRUE
    x.tick.labels <- names(results.vector)
    if (!is.null(var.dictionary)) {
      x.tick.labels[x.tick.labels%in%names(var.dictionary)] <- var.dictionary[x.tick.labels[x.tick.labels%in%names(var.dictionary)]]
    }
  }
  
  if (plot.quantity == "BayesFactor") {
    ## --- Bayes factors
    if (is.null(y.max)) {
      y.max <- ceiling(max(results.vector,na.rm=T))
    }
    y.ticks <- c(0:y.max)
    .ManhattanSetup(
      plot.title = plot.title,
      y.ticks = y.ticks,
      y.lab = "Bayes Factor",
      y.tick.labels=c("<=1",paste(10^y.ticks[-1],sep="")),
      x.ticks = x.ticks,
      x.tick.labels = x.tick.labels,
      results.vec = results.vector,
      point.cols = point.cols,
      point.labs = point.labs
      )
  } else if (plot.quantity == "PosteriorProbability") {
    ## --- Posterior probability
    y.ticks <- seq(0,1,by=0.2)    
    .ManhattanSetup(
      plot.title = plot.title,
      y.ticks = y.ticks,
      y.lab = "Posterior Probability",
      y.tick.labels=y.ticks,
      x.ticks = x.ticks,
      x.tick.labels = x.tick.labels,
      results.vec = results.vector,
      point.cols = point.cols,
      point.labs = point.labs
    )    
  } else if (plot.quantity == "pvalue") {
    ## --- p-values
    if (is.null(y.max)) {
      y.max <- ceiling( max(-log10(results.vector)) )
    }
    y.ticks <- c(0:y.max)    
    .ManhattanSetup(
      plot.title = plot.title,
      y.ticks = y.ticks,
      y.lab = "p-value",
      y.tick.labels=paste(10^-y.ticks,sep=""),
      x.ticks = x.ticks,
      x.tick.labels = x.tick.labels,
      results.vec = -log10(results.vector),
      point.cols = point.cols,
      point.labs = point.labs
    )
    # Bonferonni
    if (add.bonferroni) {
      abline(h=-log10(0.05/length(results.vector)), lty=2, col = "red")
    }    
  } else if (plot.quantity == "ABF") {
    ## --- Wakefield Bayes Factors
    if (is.null(y.max)) {
      y.max <- ceiling( max(-log10(results.vector)) )
    }
    y.ticks <- c(0:y.max)
    abf.tick.labs <- paste(10^-y.ticks,sep="")
    abf.tick.labs[1] <- "<=1"
    results.vector[results.vector>1] <- 1
    .ManhattanSetup(
      plot.title = plot.title,
      y.ticks = y.ticks,
      y.lab = "ABF for the Null",
      y.tick.labels=abf.tick.labs,
      x.ticks = x.ticks,
      x.tick.labels = x.tick.labels,
      results.vec = -log10(results.vector),
      point.cols = point.cols,
      point.labs = point.labs
    )
  }
}

############################
.ManhattanSetup <- function(
  plot.title,
  y.ticks,
  y.lab,
  y.tick.labels,
  results.vec,
  x.ticks = FALSE,
  x.tick.labels=FALSE,
  point.cols,
  point.labs) {
  
  ## --- Initiate plot
  plot(0,type="n",xlim=c(1, length(results.vec) ),
       ylim=c(min(y.ticks), max(y.ticks)),axes=F,
       xlab="",ylab=y.lab, main=plot.title)
  
  ## --- Set up X and Y axis  
  axis(side=1,at=c(1:length(results.vec)),labels=x.tick.labels, tick = x.ticks, las=2, cex.axis=0.5)
  axis(side=2,at=y.ticks, labels=y.tick.labels,las=2,col.axis="black")
  
  ## --- Add grey lines for all y.ticks
  for (l in y.ticks) {
    abline( h=l, lty=2, col = "grey" )
  }
  
  ## --- Plot statistics
  points(
    x=c(1:length(results.vec)),
    y=results.vec,
    pch=16,
    col=point.cols)
  if (!is.null(point.labs)) {
    text(
      x=which(point.labs!=""),
      y=results.vec[which(point.labs!="")],
      labels = point.labs[which(point.labs!="")]
      )
  }
}
