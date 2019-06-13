#' Returns list of covariates between which there are the strongest pairwise correlations
#' @export
#' @title Strongest Pairwise Correlation
#' @name StrongestPairwiseCorrelation
#' @param data Data matrix (or data frame), within which pairwise correlations will be calculated.
#' @param covariates Character vector indicating which columns of data should be considered. By 
#' default all columns of data will be considered.
#' 
#' @return A character vector containing the names of the covariates between which the strongest
#' pairwise correlation was observed.
#' 
#' @author Paul Newcombe
StrongestPairwiseCorrelation <- function(
  data=NULL,
  covariates=colnames(data)
  ) {
  
  cat("Calculating correlation matrix...\n")
  cor.mat <- cor(data[,covariates], use="pairwise.complete.obs")^2
  cat("...done.\n")
  
  diag(cor.mat) <- 0
  max.pairwise.cors <- apply(cor.mat, MAR=1, function(x) max(x, na.rm=T))
  max.max.pairwise.cors <- max(max.pairwise.cors)
  cat("The maximum pairwise correlation was",max.max.pairwise.cors,"\n")
  correlated.predictors <- names(which(max.pairwise.cors==max.max.pairwise.cors))
  
  return(correlated.predictors)
}
