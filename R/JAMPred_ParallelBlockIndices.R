#' For a given number of CPU cores and LD blocks, this function returns a list 
#' indicating which blocks should be analysed on each CPU core.
#' @export
#' @title List of which LD blocks are to be analysed on which CPU core
#' @name JAMPred_ParallelBlockIndices
#' @param n.cores Number of CPU cores to use
#' @param n.blocks.to.analyse Total number of LD blocks to analyse
#' 
#' @return A list of which blocks are to be analysed on each CPU core.
#' 
#' @author Paul Newcombe

JAMPred_ParallelBlockIndices <- function(
  n.cores,
  n.blocks.to.analyse
) {
  
  # --- Split blocks across parallel jobs
  block.to.core.ratio <- n.blocks.to.analyse/n.cores
  
  # --- Excess number of blocks in last group
  excess <- n.blocks.to.analyse - floor(block.to.core.ratio)*n.cores
  
  if (excess == 0 ) {
    # Number of blocks divides perfectly into number of cores
    parallel.block.indices <- lapply(
      c(1:n.cores), function(i) seq(from=(i-1)*block.to.core.ratio+1, to=i*block.to.core.ratio)
    )
    
  } else if (excess > 0) {
    parallel.block.indices <- lapply(
      c(1:(n.cores-excess)), function(i) seq(from=(i-1)*floor(block.to.core.ratio)+1, to=i*floor(block.to.core.ratio))
    )
    parallel.block.indices.excess <- lapply(
      c(1:excess), function(i) seq(from=(i-1)*ceiling(block.to.core.ratio)+1, to=i*ceiling(block.to.core.ratio)) + max(unlist(parallel.block.indices))
    )
    parallel.block.indices[((n.cores-excess)+1):n.cores] <- parallel.block.indices.excess
  }
  
  return(parallel.block.indices)
}
