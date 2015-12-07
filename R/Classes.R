#' @name R2BGLiMS_Results-class
#' @aliases R2BGLiMS_Results,R2BGLiMS_Results-class
#' 
#' @title The R2BGLiMS_Results class
#' 
#' @description Container for results from running R2BGLiMS as generated
#' by the function \code{\link[R2BGLiMS]{R2BGLiMS}}. 
#' 
#' @slot likelihood Likelihood type of model. See \code{\link[R2BGLiMS]{R2BGLiMS}} for the options.
#' @slot posterior.summary.table Posterior summaries of all parameters
#' @slot enumerate.up.to.dim If posterior inference was made by exhaustively enumerating and assessing models one by one, modles were considered
#' up to this dimension.
#' @slot n.iterations Number of iterations which the RJMCMC was run for.
#' @slot thin Ith iterations which were saved.
#' @slot model.space.prior List defining the model space prior. See \code{\link[R2BGLiMS]{R2BGLiMS}}.
#' @slot confounders Vectors of variables fixed in the model, and excluded from model selection.
#' @slot run.times A list containing run times broken down into different processes.
#' @slot n.covariate.blocks.for.jam The number of partitioned LD blocks used for JAM.
#' @slot bglims.arguments The arguments passed to the Java BGLiMS function. This is a named list - the different arguments are a range of datatypes.
#' @slot bglims.rjmcmc.output The Reversible Jump MCMC output from BGLiMS. Columns are parameters, rows are iterations.
#'   
#' @author Paul J. Newcombe \email{paul.newcombe@@mrc-bsu.cam.ac.uk}
setClass("R2BGLiMS_Results",
         representation = representation(
           likelihood = "character",
           posterior.summary.table = "matrix",
           enumerate.up.to.dim = "numeric",
           enumerated.posterior.inference = "list",
           n.iterations = "numeric",
           thin = "numeric",
           model.space.priors = "list",
           confounders = "vector",
           run.times = "list",
           n.covariate.blocks.for.jam = "numeric",
           bglims.arguments = "list",
           bglims.rjmcmc.output = "data.frame"),
         validity = function(object){
           errors <- character()           
           if (length(errors) == 0) TRUE else errors
         }
)
