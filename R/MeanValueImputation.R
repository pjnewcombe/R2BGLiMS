#' Replaces any missing data with mean values of the covariate. NOTE: This is a VERY crude
#' imputation method, for exploratory purposes only.
#' @export
#' @title Mean value imputation
#' @name MeanValueImputation
#' @param data Data set to impute missingness in.
#' @param which.vars Which covaraites to perform the mean value imputation for. If left NULL this
#' is performed for all covariates in the dataset.
#' @return data Filled in dataset.
#' @author Paul Newcombe
MeanValueImputation <- function(
  data,
  which.vars=NULL
  ) {
  
  # Perform for all variables if a list is not provided
  if (is.null(which.vars)) {
    which.vars <- colnames(data)
  }
  
  # Mean value imputation
  for (v in which.vars) {
    if (sum(is.na(data[,v]))>0) {
      data[which(is.na(data[,v])),v] <- mean(data[,v], na.rm=T)      
    }
  }
  
  return(data)
}
