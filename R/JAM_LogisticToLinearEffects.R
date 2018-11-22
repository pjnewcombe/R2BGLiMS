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
#' @param n.cases Number of cases in the dataset in which the log-odds ratios were calculated
#' 
#' @return A vector of effects on the linear scale
#' 
#' @author Paul Newcombe

JAM_LogisticToLinearEffects <- function(
  log.ors = NULL,
  log.or.ses = NULL,
  mafs = NULL,
  n = NULL,
  n.cases = NULL
) {
  
  snp.sds <- sqrt(mafs*(1-mafs))
  p.cases <- n.cases/n
  
  # Standardised effects (z-scores)
  standardised.beta.hats <- log.ors/(log.or.ses*sqrt(n.training))
  
  # Divide by SNP SDs to get allelic effects
  beta.hats.from.logistic.ps <- standardised.beta.hats/snp.sds # Divide standardised linear effects by SNP standard deviations
  
  # Adjust for fraction of cases
  beta.hats.from.logistic.ps <- beta.hats.from.logistic.ps*sqrt(p.cases*(1-p.cases)) # Multiply by trait SD for effect on trait scale
  
  return(beta.hats.from.logistic.ps)
}
