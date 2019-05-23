#' Given marginal log-ORs and corresponding MAFs, this function infers the effect estimates that would have
#' been obtained running one-at-a-time linear regressions of the binary outcome. This conversion is required for
#' JAM since it is based on a linear regression model.
#' @export
#' @title Logistic to linear effect conversion for JAM
#' @name JAM_LogisticToLinearEffects
#' @param log.ors Vector of SNP summary log-odds ratios
#' @param log.or.ses Vector of the standard errors of the SNP log-odds ratios
#' @param mafs Vector of SNP minor allele frequencies
#' @param n Size of dataset in which the log-odds ratios were calculated
#' @param p.cases Proportion of cases in the dataset in which the log-odds ratios were calculated
#' 
#' @return A vector of effects on the linear scale
#' 
#' @author Paul Newcombe

JAM_LogisticToLinearEffects <- function(
  log.ors = NULL,
  log.or.ses = NULL,
  mafs = NULL,
  n = NULL,
  p.cases = NULL
) {
  
  # Standardised least squares estimate is signed z-score/sqrt(n)
  standardised.least.squares.effect <- (log.ors/log.or.ses)/sqrt(n)
  
  # Multiply by trait SD for effect on trait scale and divide by SNP SD for per allele effect
  linear.beta.hats <- standardised.least.squares.effect*sqrt(p.cases*(1-p.cases))/sqrt(mafs*(1-mafs)) 
  
  return(linear.beta.hats)
}
