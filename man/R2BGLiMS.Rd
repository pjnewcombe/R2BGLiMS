\name{R2BGLiMS}
\alias{R2BGLiMS}
\title{Call BGLiMS from R}
\usage{
  R2BGLiMS(likelihood = NULL, data = NULL,
    outcome.var = NULL, times.var = NULL, xTx = NULL,
    t = NULL, cluster.var = NULL, confounders = NULL,
    model.selection = TRUE, model.space.priors = NULL,
    beta.priors = NULL, g.prior = FALSE, n.mil = 1,
    seed = 1, alt.initial.values = FALSE,
    results.path = NULL, results.label = "R2BGLiMS",
    args.file = NULL, debug.path = NULL)
}
\arguments{
  \item{likelihood}{Type of model to fit. Current options
  are "Logistic" (for binary data), "Weibull" (for survival
  data), "Gaussian" (for continuous data), and
  "GaussianMarg" (for analysis of univariate associations
  from Gaussian linear regressions).}

  \item{data}{Matrix or dataframe containing the data to
  analyse. Rows are indiviuals, and columns are variables,
  named in concordance with the following options.}

  \item{outcome.var}{Which column in data contains the
  binary outcome/survival status variable (default
  "Disease")}

  \item{times.var}{If survival data, the column in data
  which contains the event times (default NULL)}

  \item{xTx}{GaussianMarg ONLY: List containing each
  block's plug-in estimate for X'X.}

  \item{t}{GaussianMarg ONLY: Vector of quantities
  calculated from the summary statistics.}

  \item{cluster.var}{If hierarchical data and random
  intercepts are required, the column in data contains the
  clustering variable (default NULL)}

  \item{confounders}{vector of confounders to fix in the
  model at all times, i.e. exclude from model selection
  (default NULL)}

  \item{model.selection}{Whether to use model selection
  (default is TRUE)}

  \item{model.space.priors}{Must be specified if
  model.selection is set to TRUE. A list with as many
  elements as desired model space components (one or more).
  Each element of the list is also a list, corresponding to
  a particular model space component listing the prior
  parameters and covariates in that component. Either a
  Poisson prior on model space can be used, in which case
  an element "Rate" must be included, or Beta-Binomial
  model space priors can be used, in which case elements
  "a" and "b" must be included corresponding to the
  beta-binomial hyperparameters. Finally, an element
  "Variables" must be included - a vector of variable names
  in that component. Note that when choosing a
  beta-binomial model space prior, higher values of "b"
  relative to "a" increase sparsity, whereas higher values
  of "a" encourage the inclusion of more covariates. a ~
  number of prior successes, b~number of prior failures.}

  \item{beta.priors}{Optional matrix or data.frame
  containing two columns (1st:mean, 2nd:SD) describing
  normal priors for the beta coefficients of the
  confounders or confounders and predictors (i.e. log-ORs
  in the case of logistic regression). Rows of beta.priors
  must be named with the corresponding variable names in
  the data. NOTE: For predictors included in the model
  selection search, the interpretation of these priors is
  on the effect conditional on inclusion in the model.
  There are three options: 1) Do not specify in which case
  N(0, 1e6) priors are assigned to the betas of any
  confounders, and the effects of any predictors included
  in model selection are assigned a common normal prior
  with unknown variance, which in turn assigned a
  InversGamma(0.001, 0.001) hyperprior. 2) Provide a
  beta.priors matrix including confounders only. This
  overrides the default use of N(0,1e6) confounder priors,
  but the predictors for which model selection is performed
  for are still assigned a common prior with unknown
  variance. 3) provide a beta.priors matrix including all
  confounders and predictors, thereby avoiding the use of a
  common prior with unkown variance across predictors for
  which model selection is performed for.}

  \item{g.prior}{Whether to use a g-prior for the beta's -
  i.e. a multivariate normal with correlation structure
  proportional to X'X^-1 (default = FALSE).}

  \item{n.mil}{Number of million iterations to run (default
  is 1)}

  \item{seed}{Which random number seed to use in the RJMCMC
  sampler (default is 1)}

  \item{alt.initial.values}{Whether to use the alternative
  set of initial values in the arguments file (default is
  FALSE)}

  \item{results.path}{Optional path if wish to to save
  algorithm results output (when set to the default NULL
  all output is written to a temporary directory and
  deleted).}

  \item{results.label}{Optional label to for algorithm
  output files (if you have specified results.path)}

  \item{args.file}{Optional path to an alternative
  arguments file to use instead of the default (default is
  used if left as NULL).}

  \item{debug.path}{Optional path to save the data and
  results files (rather than as temporary files) to help in
  debugging.}
}
\value{
  A Reversible Jump results object is returned. This is a
  list of two elements: "args" which records various
  modelling arguments used in the analysis, and "results" -
  a matrix containing all saved posterior samples from the
  analysis. Columns correspond to parameters and each row
  contains values from a particular itertation, with 0
  indicating exclusion from the model.

  The function \code{\link{ResultsTable}} can be used to
  print summary posterior results for all parameters. Other
  functions for summarising results are listed under "see
  also".
}
\description{
  Calls BGLiMS - a Java package for fitting GLMs under
  Bayesian model selection. NOTE: The predictors to explore
  with model selection are specified via the
  model.space.priors argument - see examples. By Default a
  common, and unknown, prior SD is used for all predictors
  under model selection, which in turn is assigned an
  Inverse-Gamma hyper-prior. Fixed normal priors may be
  user specified for all predictor (and confounder)
  coefficients via the beta.priors argument.
}
\examples{
library(R2BGLiMS)

### --- Logistic regression with two model space prior partitions
# Example dataset from MASS
utils::data(biopsy, package = "MASS")
# Normalise predictors in order to use (default) common unknown priors on their effect sizes
for (v in paste("V",c(1:9),sep="")) {
  biopsy[,v] <- (biopsy[,v]-mean(biopsy[,v],na.rm=T))/sd(biopsy[,v],na.rm=T)
}
biopsyResults <- R2BGLiMS(
  likelihood="Logistic", data=biopsy, outcome.var="class",
  model.space.priors=list(
    list("Rate"=0.5, "Variables"=paste("V",c(1:3),sep="")),
    list("Rate"=0.1, "Variables"=paste("V",c(4:9),sep=""))
  )
)
PrettyResultsTable(biopsyResults)
TopModels(biopsyResults)

### --- Survival analysis with priors specified on effect sizes
# Example dataset from MASS
utils::data(VA, package = "MASS")
# Generate priors on log-effect sizes
predictors <- c("treat","age","Karn","diag.time","cell","prior")
beta.priors <- data.frame(
  cbind(rep(0,length(predictors)),rep(10,length(predictors))),
  row.names=predictors)
# Notice how survival analysis requires additionally flagging a "times.var"
VAResults <- R2BGLiMS(
  likelihood="Weibull",
  data=VA,
  outcome.var="status",
  times.var="stime",
  model.space.priors=list(list("a"=1,"b"=1,"Variables"=predictors)),
  beta.prior=beta.priors
)
PrettyResultsTable(VAResults)
TopModels(VAResults)
}
\author{
  Paul Newcombe
}
\seealso{
  For summary results of covariates see
  \code{\link{ResultsTable}},
  \code{\link{PrettyResultsTable}} and
  \code{\link{ManhattanPlot}}. For posterior model space
  summaries see \code{\link{ModelSizes}} and
  \code{\link{TopModels}}. For convergence check plots see
  \code{\link{ChainPlots}} and
  \code{\link{AutocorrelationPlot}}.
}

