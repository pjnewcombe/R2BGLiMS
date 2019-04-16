#' Manhattan plotting function for results from \code{\link[R2BGLiMS]{JAMMR}} analyses.
#' @export
#' @title JAMMR_ManhattanPlot
#' @name JAMMR_ManhattanPlot
#' @param jammr.results A \code{\link[R2BGLiMS]{JAMMR}} results object frum running the \code{\link[R2BGLiMS]{JAMMR}} command.
#' @param only.best APOSTOLOS TO PROVIDE
#' 
#' @author Apostolos Gkatzionis

JAMMR_ManhattanPlot <- function (jammr.results, only.best = FALSE) {
  
  ## w = 1 is the same as plotting only the best run, so acknowledge that.
  if (length(jammr.results$all.w) == 1) only.best <- TRUE
  
  if (only.best) {
    ## If "only.best", we plot a single Manhattan Plot with the best-w run.
    
    ## Set up.
    P <- length(jammr.results$snp.probs)
    par(mfrow = c(1, 1))
    
    ## Plot SNP inclusion probabilities.
    plot(1:P, jammr.results$snp.probs, xlab = "SNP", ylab = "Inclusion Probability", main = "Manhattan Plot", ylim = c(0, 1), pch = 19, axes = FALSE)
    
    ## Plot axes and horizontal lines.
    abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 2)
    axis(side = 1, at = 1:P, labels = paste("SNP", c(1:P), sep = ""), las = 2, cex.axis = 0.5)
    axis(side = 2, at = seq(0, 1, length = 6))
    
  } else {
    ##  If "only.best" = FALSE, we plot everything.
    
    ## Get the number of SNPs and w's.
    P <- length(jammr.results$snp.probs)
    nw <- length(jammr.results$all.w)
    
    ## Set up par-mfrow.
    if (nw %in% c(1, 4, 9, 16, 25, 36, 49, 64, 81, 100)) {
      get.mfrow <- c(floor(sqrt(nw)), floor(sqrt(nw)))
    } else if (floor(sqrt(nw)) * (floor(sqrt(nw)) + 1) >= nw){
      get.mfrow <- c(floor(sqrt(nw)), floor(sqrt(nw)) + 1)
    } else {
      get.mfrow <- c(floor(sqrt(nw)) + 1, floor(sqrt(nw)) + 1)
    }
    par(mfrow = get.mfrow)
    
    ## Start doing the Manhattan Plots.
    for (I in 1:nw) {
      
      ## Plot SNP inclusion probabilities.
      plot(1:P, jammr.results$all.probs[I, ], xlab = "SNP", ylab = "Inclusion Probability", main = "Manhattan Plot", ylim = c(0, 1), pch = 19, axes = FALSE)
      
      ## Plot axes and horizontal lines.
      abline(h = c(0, 0.2, 0.4, 0.6, 0.8, 1), col = "grey", lty = 2)
      axis(side = 1, at = 1:P, labels = paste("SNP", c(1:P), sep = ""), las = 2, cex.axis = 0.5)
      axis(side = 2, at = seq(0, 1, length = 6))
      
    }
  }
}
