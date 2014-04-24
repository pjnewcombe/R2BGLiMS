#' Generates splines to append to original dataset, for spline modelling. Knots are created at quantiles,
#' and binary indicators of being greater than or equal to the knot are returned. NOTE: if a covariate
#' does not have enough unique quantiles it will be excluded from the spline modelling, and
#' returned in it's original form.
#' @export
#' @title Generates splines.
#' @name GenerateSplines
#' @param data Data set to append splines to
#' @param which.vars Which variables to generate splines for. If left at default NULL
#' @param n.knots How many splines to generate for each variable (knots will be placed at quantiles). NOTE:
#' if a covariate does not have enough unique quantiles it will be excluded from the spline modelling, and
#' returned in it's original form.
#' @return spline.data A dataset containing the splines. Can be merged to the original data using cbind, and
#' the names of the spline covariates can be obtained using colnames.
#' @author Paul Newcombe
GenerateSplines <- function(
  data,
  which.vars=NULL,
  n.knots=6
  ) {
  
  if (is.null(which.vars)) {
    which.vars <- colnames(data)
  }
  all.vars <- which.vars

  # Figure out quantiles to place knots at
  quantiles <- matrix(NA,length(which.vars),n.knots,dimnames=list(which.vars,NULL))
  for (v in which.vars) {
    quantiles[v,] <- quantile(data[,v], p=seq(from=0,to=1,length.out=n.knots+2)[-c(1,n.knots+2) ],na.rm=T)
  }
  
  # Check distinct quantiles
  n.distinct <- apply(quantiles[which.vars,], MAR=1, function(x) length(unique(x)))
  # Exclude any variables without enough distinct quantiles from spline modelling - avoids colinearity
  which.vars <- which.vars[!which.vars%in%which.vars[which(n.distinct<n.knots)]]
  
  # Check bottom is not minimum
  mins <- apply(data[,which.vars], MAR=2, min)
  # Exclude any variables whose bottom spline is the minimum - avoids colinearity
  which.vars <- which.vars[!which.vars%in%which.vars[which(mins==quantiles[which.vars,1])]]

  # Append splines
  spline.vars <- NULL
  for (q in c(1:n.knots)) {
    data.add <- data[,which.vars]
    for (v in which.vars) {
      data.add[,v] <- as.integer(data.add[,v]>=quantiles[v,q])
    }
    colnames(data.add) <- paste(colnames(data.add),"_geQ",q,sep="")
    data <- cbind(data, data.add)
    spline.vars <- c(spline.vars, colnames(data.add))
  }
  
  # Pull out splines (and original covariates for any variables excluded from spline modelling)
  keep.vars <- c(all.vars[!all.vars%in%which.vars], spline.vars)
  spline.data <- data[,keep.vars]
  
  # Return
  return(spline.data)
}
