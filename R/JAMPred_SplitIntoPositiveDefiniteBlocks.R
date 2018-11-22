#' For a list of snps and corresponding individual level reference genotypes, this function can be used to split
#' the SNPs into blocks of positive definite correlation structure, i.e. for use with JAMPred. After starting with
#' an initial partitioning of the SNPs, either provided according to LD information, or naively formed based on 
#' a user specified initial block size, this function identifies any blocks which do not have infertible genotype
#' matrices and splits them until all blocks have postivie definite genotype matrices in the reference data.
#' @export
#' @title Splits a list of SNPs into positive definite blocks.
#' @name JAMPred_SplitIntoPositiveDefiniteBlocks
#' @param snps List of SNPs to be split into blocks
#' @param initial.block.size Block size for the initial partitioning. Default is 100. This is ignored when initial.snps.blocks is given.
#' @param initial.snps.blocks A user specficied initial partitioning of the SNPs. This must be a list, each element of which is a vector of SNP
#' names that correspond to the columns of ref.geno. This could be defined, for example, according to LD information.
#' @param ref.geno Reference genotype matrix. Rows are people, columns are SNPs (named according to the vector snps) and elements are between 
#' 0 and 2 (can be dosages).
#' 
#' @return A list each element of which is a vector of SNP names.
#' 
#' @author Paul Newcombe

JAMPred_SplitIntoPositiveDefiniteBlocks <- function(
  snps,
  initial.block.size = 100,
  initial.snps.blocks = NULL,
  ref.geno) {
  
  if (is.null(initial.snps.blocks)) {
    # --- Initial split into blocks by the target block size
    block.splits <- seq(from=1, to=length(snps), by=initial.block.size)
    
    # Code below requires final element of splits to be nSNPs+1 (does not always happen)
    if (max(block.splits)<length(snps)) block.splits <- c(block.splits, (length(snps)+1) )
    if (max(block.splits)==length(snps)) block.splits[length(block.splits)] <- length(snps) + 1
    snps.blocks <- sapply(c(1: (length(block.splits)-1) ),  function(b) snps[block.splits[b]:(block.splits[b+1]-1)],simplify=FALSE) # Disable simplify so always a list
  } else {
    snps.blocks <- initial.snps.blocks
  }
  
  # --- Check RANK
  block.rank.checks <- unlist(lapply(snps.blocks, function(b) JAM_RankCheck(X.ref = ref.geno[,b]) ))
  if (sum(block.rank.checks) == length(block.rank.checks)) {
    cat("All blocks are full rank")  
  } else {
    # --- Deal with blocks which are not full rank
    problem.blocks <- snps.blocks[!block.rank.checks]
    # Try further splitting problem blocks in a while loop until they are all full rank
    while (sum(unlist(lapply(problem.blocks, function(b) JAM_RankCheck(X.ref = ref.geno[,b]) ))) != length(problem.blocks) ) {
      still.not.full.rank <- which(unlist(lapply(problem.blocks, function(b) JAM_RankCheck(X.ref = ref.geno[,b]) )) == FALSE)
      for (b in still.not.full.rank) {
        snps.in.this.block <- problem.blocks[[b]]
        if (length(snps.in.this.block)>=4) {
          problem.blocks[[b]] <- snps.in.this.block[1:round(length(snps.in.this.block)/2)] # Split the problem block
          problem.blocks[[length(problem.blocks)+1]] <- snps.in.this.block[!snps.in.this.block%in%problem.blocks[[b]]]
        } else if (length(snps.in.this.block)==3) {
          problem.blocks[[b]] <- snps.in.this.block[1:2] # Just take first 2 SNPs and remove the 3rd (can not have a 1 SNP block)
        } else if (length(snps.in.this.block)==2) {
          problem.blocks[[b]] <- NULL # Delete the block if it is still not full rank at 2 SNPs
        }
      }
    }
    # --- Put the problem blocks back
    snps.blocks <- snps.blocks[block.rank.checks]
    snps.blocks[(length(snps.blocks)+1):(length(snps.blocks)+length(problem.blocks))] <- problem.blocks
    # --- Re-order
    snps.blocks <- snps.blocks[order(unlist(lapply(snps.blocks, function(b) which(snps==b[1]))))]
    snps <- unlist(snps.blocks) # might have dropped some SNPs
  }
  
  return(snps.blocks)
}
