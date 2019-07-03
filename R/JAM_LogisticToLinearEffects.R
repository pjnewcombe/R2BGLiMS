#' Given marginal log-ORs and corresponding MAFs, this function infers the effect estimates that would have
#' been obtained running one-at-a-time linear regressions of the binary outcome. This conversion is required for
#' JAM since it is based on a linear regression model.
#' @export
#' @title Logistic to linear effect conversion for JAM
#' @name JAM_LogisticToLinearEffects
#' @param log.ors Vector of SNP summary log-odds ratios
#' @param log.or.ses Vector of the standard errors of the SNP log-odds ratios
#' @param snp.genotype.sds Vector of SNP genotype standard deviations. Please supply this or mafs. This is 
#' the preferred option since an assumption of Hardy-Weinberg Equilibrium is not required when
#' reverse-stamdardising the effects during the transform.
#' @param mafs Vector of SNP minor allele frequencies. snp.sds if available would be the preferred option.
#' @param n Size of dataset in which the log-odds ratios were calculated
#' @param p.cases Proportion of cases in the dataset in which the log-odds ratios were calculated
#' 
#' @return A vector of effects on the linear scale
#' 
#' @author Paul Newcombe

JAM_LogisticToLinearEffects <- function(
  log.ors = NULL,
  log.or.ses = NULL,
  snp.genotype.sds = NULL,
  mafs = NULL,
  n = NULL,
  p.cases = NULL
) {
  
  # Standardised least squares estimate is signed z-score/sqrt(n)
  standardised.least.squares.effect <- (log.ors/log.or.ses)/sqrt(n)
  
  if (!is.null(mafs) & !is.null(snp.genotype.sds)) {
    cat("snp.sds and mafs were provided. snp.genotype.sds will be used as the preferred option.\n")
    mafs <- NULL
  }
  
  if (!is.null(mafs)) {
    # Multiply by trait SD for effect on trait scale and divide by SNP SD for per allele effect
    linear.beta.hats <- standardised.least.squares.effect*sqrt(p.cases*(1-p.cases))/sqrt(2*mafs*(1-mafs))
  }
  
  if (!is.null(snp.genotype.sds)) {
    # Multiply by trait SD for effect on trait scale and divide by SNP SD for per allele effect
    linear.beta.hats <- standardised.least.squares.effect*sqrt(p.cases*(1-p.cases))/snp.genotype.sds
  }
  
  return(linear.beta.hats)
}
