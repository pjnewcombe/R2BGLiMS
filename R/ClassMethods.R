#' @name R2BGLiMS_Results-methods
#' @aliases show, R2BGLiMS_Results-method
#' 
#' @title 'show' method for \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} objects
#' 
#' @description Prints a summary of the \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} object.
#' 
#' @return Prints a summary of the \code{\link[R2BGLiMS]{R2BGLiMS_Results-class}} object. 
#' 
#' @author Paul J. Newcombe \email{paul.newcombe@@mrc-bsu.cam.ac.uk}
setMethod("show",
          signature = "R2BGLiMS_Results",
          definition = function(object){
            cat("An object of class ",class(object),".\n",sep="")
            cat("Data were analysed under the",object@likelihood,"likelihood.\n")
            cat("Model selection was performed for",length(unlist(lapply(object@model.space.priors, function(x) x$Variables))),"covariates.\n")
            if (object@enumerate.up.to.dim>0) {
              cat("Inference was obtained by exhaustively enumerating models up to dimension",object@enumerate.up.to.dim,"\n")              
            } else {
              cat ("Inference was obtained by running RJMCMC for",object@n.iterations,"iterations.\n")
            }
            cat("A posterior summary table of all parameters is contained in the slot 'posterior.summary.table'")
          })
