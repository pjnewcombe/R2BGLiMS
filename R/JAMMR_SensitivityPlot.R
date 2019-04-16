#' Plots causal effect estimates and confidence intervals for multiple w runs
#' from \code{\link[R2BGLiMS]{JAMMR}} analyses.
#' @export
#' @title JAMMR_SensitivityPlot
#' @name JAMMR_SensitivityPlot
#' @param jammr.results A \code{\link[R2BGLiMS]{JAMMR}} results object frum running the \code{\link[R2BGLiMS]{JAMMR}} command.
#' @param adjust.se APOSTOLOS TO PROVIDE
#' @param line.at.zero Whether to include a line at zero
#' @param show.n.snps APOSTOLOS TO PROVIDE
#' 
#' @author Apostolos Gkatzionis

JAMMR_SensitivityPlot <- function(
  jammr.results,
  adjust.se = TRUE,
  line.at.zero = TRUE,
  show.n.snps = FALSE) {
  
  ## This should only take as input a JAM-MR run.
  
  ## Store the useful parts of the JAM-MR run.
  means <- jammr.results$all.causal
  stderr <- jammr.results$all.se
  nw <- length(jammr.results$all.w)
  if (nw == 1) warning("Sensitivity analysis plot will be created for only one w value.")
  
  ## If requested, adjust standard errors.
  if (adjust.se) stderr <- stderr * 1.3178
  
  ## Compute plot limits for the y-axis.
  y.axis.limits <- c(min(means - 1.96 * stderr) - 0.1, max(means + 1.96 * stderr) + 0.1)
  
  ## Plot the points and the lines.
  plot(x = 1:nw, y = means, type = "p", axes = FALSE, main = "Sensitivity Plot", xlab = "w / N1", ylab = "Causal Effect", ylim = y.axis.limits, pch = 19)
  for (i in 1:nw) lines(c(i, i), c(means[i] - 1.96 * stderr[i], means[i] + 1.96 * stderr[i]))
  
  ## Plot the axes.
  axis(side = 1, at = 1:nw, labels = as.character(jammr.results$all.w), las = 1, cex.axis = 0.7)
  axis(side = 2)
  
  ## Additions.
  if (line.at.zero) abline(h = 0, lty = 2, col = "brown")
  if (show.n.snps) text( x = 1:nw, y = rep(y.axis.limits[1], nw), labels = round(rowSums(jammr.results$all.probs)), cex = 0.8, col = "black" )
  
}
