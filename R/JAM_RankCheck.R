#' Checks the rank of a reference genotype matrix, which is required to run JAM.
#' @export
#' @title JAM rank check.
#' @name JAM_RankCheck
#' @inheritParams R2BGLiMS
#' 
#' @return TRUE/FALSE of whether the matrix is full rank (hopefully it is and TRUE is returned). If X.ref is a list
#' of reference matrices for different regions, a vector of boolen TRUE/FALSE values is returned, one for each LD
#' block.
#' 
#' @seealso See also \code{\link{JAM}}.
#' 
#' @author Paul Newcombe

JAM_RankCheck <- function(
  X.ref=NULL
) {
  
  # If data.frame make matrix to assess length
  if (is.data.frame(X.ref)) {
    X.ref <- as.matrix(X.ref) # convert to matrix
  }
  
  # If is not a list after above, make it a list of length 1
  if (!is.list(X.ref)) {
    X.ref <- list(X.ref) # convert to list if there is a single block
  }
  
  # Iniitiate vector of full rank true/false
  full.rank.vector <- rep(FALSE, length(X.ref))
  
  for (ld.block in 1:length(X.ref)) {
    cat (paste("Taking the QR decomposition of block ",ld.block,"...\n",sep="") )
    qr.decomp <- qr(X.ref[[ld.block]])
    if (qr.decomp$rank == ncol(X.ref[[ld.block]])) {
      full.rank.vector[ld.block] <- TRUE
      cat ("...Done.\nBlock",ld.block,"is full rank.\n")
    } else {
      cat ("...Done.\nBlock",ld.block,"is NOT full rank. Try further pruning?\n")
    }
  }
  
  return(full.rank.vector)
}