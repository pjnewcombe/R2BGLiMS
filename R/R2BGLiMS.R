#' @include ResultsTable.R
NULL

#' Calls BGLiMS - a Java package for fitting GLMs under Bayesian model selection. NOTE: The predictors to explore with model
#' selection are specified via the model.space.priors argument - see examples. By Default a common, and unknown, prior SD
#' is used for all predictors under model selection, which in turn is assigned an Inverse-Gamma hyper-prior.
#' Fixed normal priors may be user specified for all predictor (and confounder) coefficients via the beta.priors argument.
#' 
#' @export
#' @title Call BGLiMS from R
#' @name R2BGLiMS
#' @param likelihood Type of model to fit. Current options are 
#' "Logistic" (for binary data), 
#' "Weibull" (for survival data), 
#' "Cox" (for survival data), 
#' "RocAUC" (to optimise ROC AUC),
#' "Gaussian" (for linear regression), 
#' "GaussianConj" (linear regression exploiting conjugate results), 
#' "GaussianMarg" (for linear regression summary statistics) 
#' and "GaussianMargConj" (for analysis under a marginal conjugate linear regression model).
#' @param data Matrix or dataframe containing the data to analyse. 
#' Rows are indiviuals, and columns contain the variables and outcome.
#' If modelling summary statistics specify xTx instead (see below).
#' @param outcome.var Name of outcome variable in data. For survival data see times.var below.
#' If modelling summary statistics specify z instead (see below).
#' @param times.var SURVIVAL DATA ONLY Name of column in data which contains the event times.
#' @param confounders Optional vector of confounders to fix in the model at all times, i.e. exclude from model selection.
#' @param model.selection Whether to use model selection (default is TRUE). NB: Even if set to FALSE, please provide a 
#' dummy model.space.priors argument (see below). This will not be used mathematically, but the package requires taking the list of variable
#' names from the "Variables" element therein.
#' @param model.space.priors Must be specified if model.selection is set to TRUE.
#' A list with as many elements as desired model space `components'.
#' Each element of the list is also a list, listing the prior parameters and 
#' covariates in a particular model space component. 
#' Either a Poisson prior on model space can be used, in which case an
#' element "Rate" must be included, or Beta-Binomial model space priors can
#' be used, in which case elements "a" and "b" must be included corresponding
#' to the beta-binomial hyperparameters. Higher values of "b" relative to "a" 
#' encourage sparsity.
#' A vector "Variables" must be included listing the variable names in that component. 
#' @param beta.priors There are three options: 
#' 1) Leave as null, in which case N(0, 1e6) priors are placed on the effects
#' of confounders and any predictors included in model selection are assigned a
#' common normal prior with unknown variance.
#' 2) Provide fixed priors for the confounders only. 
#' 3) Provide fixed priors for all covariates. 
#' beta.priors (a matrix/data.frame) can be used to provide fixed priors; rows must be named with the corresponding 
#' variable names in the data, and include Guassian prior means and sds in the first and 
#' second columns respectively.
#' @param g.prior CONJUGATE LINEAR REGRESSION ONLY: Gaussian conjugate models ONLY: Whether to use a g-prior for the beta's - i.e. a multivariate normal 
#' with correlation structure proportional to sigma^2*X'X^-1 or to use independence priors (default = FALSE).
#' @param model.tau CONJUGATE LINEAR REGRESSION ONLY: Gaussian conjugate models ONLY: Whether to model tau or not (default FALSE). If set to true,
#' then a Zellner-Siow prior is used, centred on the value provide by tau. The value provided in tau is also used as the
#' initial value.
#' @param tau CONJUGATE LINEAR REGRESSION ONLY: Value to use for sparsity parameter tau (tau*sigma^2 parameterisation).
#' Default residual.var.n. If modelling this parameter, this value is used to center the Zellner-Siow prior
#' and as an intial value.
#' @param enumerate.up.to.dim CONJUGATE LINEAR REGRESSION ONLY: Whether to calculate posterior scores
#' for every possible model of dimension up to the integer specified (as an alternative to RJMCMC). 
#' Currenly maximum allowed dimension is 5. 
#' Leaving at 0 means this is not carried out. 
#' @param xTx SUMMARY STATISTICS ONLY: List containing each block's external plug-in estimate for X'X.
#' @param z SUMMARY STATISTICS ONLY: t(X)*y vector, calculated from the summary statistics.
#' @param n.mil Number of million iterations to run (default is 1)
#' @param seed Which random number seed to use in the RJMCMC sampler.
#' @param do.chain.plot Whether to produce a PDF containing chain plots. Default FALSE.
#' @param results.path Optional path if wish to to save algorithm output.
#' @param results.label Optional label for algorithm output files (if you have specified results.path).
#' @param extra.arguments A named list of any additional arguments for BGLiMS. Type "data(DefaultArguments)" and look in the 
#' "default.arguments" list to see the names (which must match) and default values of available extra arguments. Currently these 
#' are not documented - please contact the package author, Paul Newcombe for details.
#' @param initial.model Optionally, an initial model can be provided as a vector of 0's and 1's. Default is NULL
#' and the null model is used. If set to 1, the saturated model is used.
#' @param max.model.dim Optional specification of maximum model dimension (default -1 means no maximum is set).
#' @param debug.path Optional path to save the data and results files to (rather than as temporary files), for aid in
#' debugging.
#' 
#' @return A Reversible Jump results object is returned. This is a list of two elements: "args" which records various modelling
#' arguments used in the analysis, and "results" - a matrix containing all saved posterior samples from the analysis. Columns
#' correspond to parameters and each row contains values from a particular itertation, with 0 indicating exclusion from the model.
#' 
#' The function \code{\link{PrettyResultsTable}} can be used to print summary posterior results for all parameters. Other functions
#' for summarising results are listed under "see also".
#' 
#' @seealso For summary results of covariates see \code{\link{ResultsTable}},
#' \code{\link{PrettyResultsTable}} and \code{\link{ManhattanPlot}}. For posterior model space
#' summaries see \code{\link{ModelSizes}} and \code{\link{TopModels}}. For convergence check
#' plots see \code{\link{ChainPlots}} and \code{\link{AutocorrelationPlot}}.
#' 
#' @author Paul Newcombe
#' 
#' @example Examples/R2BGLiMS_Examples.R

R2BGLiMS <- function(
  likelihood=NULL,
  data=NULL,
  outcome.var=NULL,
  times.var=NULL,
  confounders=NULL,
  model.selection=TRUE,
  model.space.priors=NULL,
  beta.priors=NULL,
  g.prior=FALSE,
  model.tau=FALSE,
  tau=max(nrow(data),ncol(data)^2,length(z)^2),
  enumerate.up.to.dim=0,
  xTx=NULL,
  z=NULL,
  n.mil=1,
  seed=1,
  do.chain.plot=FALSE,
  results.path=NULL,
  results.label="R2BGLiMS",
  extra.arguments=NULL,
  initial.model=NULL,
  max.model.dim=-1,
  debug.path=NULL
) {
  
  ###########################
  ### --- Old options --- ###
  ###########################
  cluster.var <- NULL # Very out of date; probably broken
  alt.initial.values <- FALSE # Now done using the extra.arguments option
  
  ##############################
  ##############################
  ### --- Error messages --- ###
  ##############################
  ##############################
  
  ### --- Java installation error messages
  try.java <- try(system("java -version"), silent=TRUE)
  if (try.java!=0) stop("Java is not installed and is required to run BGLiMS.\nPlease install from java.com/download.")  
  
  ### --- Basic input checks
  if (is.null(likelihood)) stop("No likelihood, i.e. the model type, has been specified; please specify as Logistic,
                                Weibull, Cox, Gaussian, GaussianConj, RocAUC, GaussianMarg or GaussianMargConj")
  if (!is.null(likelihood)) {
    if (!likelihood %in% c("Logistic", "Weibull", "Cox", "Gaussian", "GaussianConj", "GaussianMarg", "GaussianMargConj", "RocAUC", "RocAUC_Testing")) {
      stop("Likelihood must be specified as Logistic, Weibull, Cox, Gaussian, GaussianConj, RocAUC, GaussianMarg, or GaussianMargConj")
    }
  }
  if (is.null(data)&is.null(xTx)) stop("The data to analyse has not been specified")
  if (is.null(outcome.var)&is.null(z)) stop("An outcome variable has not been specified")
  
  ### --- Likelihood specific checks
  if (likelihood %in% c("Logistic", "Weibull", "RocAUC", "RocAUC_Testing")) {
    # This check is not done for Cox - since with uncensored data the outcome is not binary
    if (is.factor(data[,outcome.var])) {
      data[,outcome.var] <- as.integer(data[,outcome.var])-1
    } else if (is.character(data[,outcome.var])) {
      data[,outcome.var] <- as.integer(as.factor(data[,outcome.var]))-1    
    }    
    if (length(table(data[,outcome.var]))!=2) stop("Outcome variable must be binary")    
  }
  if (length(grep("Conj", likelihood)) > 0) {
    if (is.null(tau)) stop("Must specify a value for tau, the variable selection coefficient. Recommend max(N, P^2).")    
  }
  
  ### --- Enumeration error messages
  if (enumerate.up.to.dim>0) {
    if (model.tau) stop ("Tau must be fixed to enumerate model specific posterior
                         scores.")
    if ((enumerate.up.to.dim>5)) stop ("Currenly only possible to enumerate models
                                     up to dimension 5") # If change edit help above    
  }  
  
  ### --- Model space prior error messages
  if (is.null(model.space.priors)) stop("Must specify a prior over the model space (even if model selection
                                        will not be used this is used to provide the list of covariate
                                        names")
  if (!is.null(model.space.priors)) {
    if (!is.list(model.space.priors)) stop("model.space.priors must be a list of list(s). If only a single model space component is required, this must still be a list of a list")
    # Check structure
    for (c in 1:length(model.space.priors)) {
      no.prior.family <- TRUE
      if ("a"%in%names(model.space.priors[[c]])&"b"%in%names(model.space.priors[[c]])) {
        no.prior.family <- FALSE        
      }
      if ("Rate"%in%names(model.space.priors[[c]])) {
        no.prior.family <- FALSE        
      }
      if (no.prior.family) {
        stop("Each model.space.prior list(s) must contain either named a and b elements to sepecify a beta-binomial prior, or a named Rate element to specify a Poisson prior")
      }
      if (!"Variables"%in%names(model.space.priors[[c]])) {
        stop("Each model.space.prior list(s) must contain an element named Variables, defining which covariates to search over")
      }
      if (!likelihood %in% c("GaussianMarg", "GaussianMargConj")) {
        if (sum(model.space.priors[[c]]$Variables%in%colnames(data))!=length(model.space.priors[[c]]$Variables)) {
          stop(paste("Not all variables in model space component",c,"are present in the data"))
        }
        # Sort out any factors
        for (v in model.space.priors[[c]]$Variables) {
          if (is.factor(data[,v])) {
            data[,v] <- as.integer(data[,v])
          } else if (is.character(data[,v])) {
            data[,v] <- as.integer(as.factor(data[,v]))-1
          }
        }        
      }
    }
  }
  
  ### --- Confounders error messages
  if (!is.null(confounders)) {
    if (sum(confounders%in%colnames(data))!=length(confounders)) stop("One or more confounders are not present in the data")
    if (sum(confounders%in%unlist(model.space.priors))>0) stop("One or more confounders are also declared in the model space")
    for (v in confounders) {
      if (is.factor(data[,v])) {
        data[,v] <- as.integer(data[,v])
      } else if (is.character(data[,v])) {
        data[,v] <- as.integer(as.factor(data[,v]))        
      }
    }
  }
  
  ### -- Beta prior error messages
  if (!is.null(beta.priors)) {
    if (length(grep("Conj",likelihood))>0) {stop("Fixed priors for the coefficients can not be specified for the conjugate model; the prior
                                                 is take as a function of X'X")}
    beta.priors.not.mat <- TRUE
    if (is.data.frame(beta.priors)) {beta.priors.not.mat <- FALSE}
    if (is.matrix(beta.priors)) {beta.priors.not.mat <- FALSE}
    if (beta.priors.not.mat) stop("beta.priors must be a matrix or data frame")
    if (ncol(beta.priors)!=2) stop("beta.priors must have two columns - 1st for means, 2nd for SDs")
    if ( is.null(rownames(beta.priors)) ) stop("Rows of beta.priors must be named with corresponding variable names")
    if (!likelihood %in% c("GaussianMarg", "GaussianMargConj")) {
      if (sum(rownames(beta.priors)%in%colnames(data))!=nrow(beta.priors)) stop("One or more variables in beta.priors are not present in the data")
    }
  }  
  
  # Setup file paths
  pack.root <- path.package("R2BGLiMS")
  bayesglm.jar <- paste(pack.root,"/BGLiMS/BGLiMS.jar",sep="")
  
  ###########################
  ###########################
  ### --- Prior setup --- ###
  ###########################
  ###########################
  
  # Deal with confounders for conjugate model
  if (!is.null(confounders)&(likelihood=="GaussianConj")) {
    cat("Obtaining confounder adjusted residuals...\n")
    data <- .ConfounderAdjustedResiduals(data, outcome.var, confounders)
    confounders <- NULL
  }
    
  # Establish model space prior component sizes
  n.mod.space.comps <- length(model.space.priors)
  mod.space.comp.sizes <- NULL
  if (model.selection) {
    for (c in 1:n.mod.space.comps) {
      mod.space.comp.sizes <- c(mod.space.comp.sizes, length(model.space.priors[[c]]$Variables))      
    }
  }
  # Establish model space prior distribution family
  if ("Rate" %in% names(model.space.priors[[1]]) ){
    modSpaceDistribution <- 0    # Poisson
  } else {
    modSpaceDistribution <- 1    # Beta-binomial
  }
  # Predictors to keep in the model at all times go at the front
  predictors <- confounders
  for (i in 1:n.mod.space.comps) {
    predictors <- c(predictors, model.space.priors[[i]][["Variables"]])
  }
  if (!is.null(beta.priors)) {
    # Establish whether just for confounders or for all predictors
    # Check order of rows
    beta.priors <- as.data.frame(beta.priors)
    if (nrow(beta.priors)==length(confounders)) {
      # Fixed priors have been provided for confounders only.
      beta.priors <- beta.priors[confounders,] # Ensure order is correct
    } else if (nrow(beta.priors)==length(predictors)) {
      # Fixed priors have been provided for everything turn beta.sd.prior off
      beta.priors <- beta.priors[predictors,] # Ensure order is correct
    } else {
      stop("beta.priors matrix must contain either all confounders or all predictors (including any confounders)")
    }
  } else if (!is.null(confounders)) {
    # If confoudners have been specified but no beta.priors, create - 
    # must have fixed priors for confounders
    beta.priors <- data.frame(
      cbind(rep(0,length(confounders)),rep(1000,length(confounders))),
      row.names=confounders)    
  }
  
  ### --- Write BGLiMS Arguments
  #arguments.file <- paste(pack.root,"/BGLiMS/Arguments_Package.txt",sep="")      
  t1 <- proc.time()["elapsed"]
  clean.up.data <- FALSE
  clean.up.results <- FALSE
  clean.up.arguments <- FALSE
  now <-format(Sys.time(), "%b%d%H%M%S") # Used to ensure unique names in the temp directory
  # Set arguments
  load(paste(pack.root,"/data/DefaultArguments.rda",sep=""))
  if (!is.null(extra.arguments)) {
    for (arg in names(extra.arguments)) {
      cat("Setting user specified argument ",arg,"\n")    
      default.arguments[[arg]] <- extra.arguments[[arg]]
    }    
  }
  # Setup arguments file
  if (!is.null(debug.path)) {
    arguments.path <- paste(debug.path,"/",results.label,"_Arguments_",now, sep="")
  } else {
    arguments.path <- paste(tempdir(),"/",results.label,"_Arguments_",now, sep="")      
    clean.up.arguments <- TRUE
  }
  try(system(paste("mkdir '",arguments.path,"'", sep="")))
  arguments.root <- paste(arguments.path, results.label, sep="/")
  arguments.file <- paste(arguments.root,"_Arguments.txt",sep="")
  # Write arguments
  write(paste(names(default.arguments)[1],default.arguments[[1]]), file = arguments.file)    
  for (arg in names(default.arguments)[-1]) {
    write(paste(arg,format(default.arguments[[arg]],sci=F)), file = arguments.file, append = T)    
  }
  rm(default.arguments)
  
  ### --- Write data
  data.file <- NULL
  if (is.null(data.file)) {
    cat("\nWriting temporary data files...\n")
    if (!is.null(debug.path)) {
      data.path <- paste(debug.path,"/",results.label,"_Data_",now, sep="")      
    } else {
      clean.up.data <- TRUE
      data.path <- paste(tempdir(),"/",results.label,"_Data_",now, sep="")      
    }
    system(paste("mkdir '",data.path,"'", sep=""))
    data.file <- paste(data.path, "/",results.label,".txt", sep="")
    .WriteData(
      data.file=data.file,
      likelihood=likelihood,
      data=data,
      predictors=predictors,
      confounders=confounders, 
      outcome.var=outcome.var,
      times.var=times.var,
      xTx=xTx,
      z=z,
      g.prior=g.prior,
      tau=tau,
      model.tau=model.tau,
      enumerate.up.to.dim=enumerate.up.to.dim,
      cluster.var=cluster.var,
      beta.priors=beta.priors,
      model.space.priors=model.space.priors,
      initial.model=initial.model
      )
    cat("\n...finished writing temporary data files.\n")
  }
  t2 <- proc.time()["elapsed"]
  write.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("\nData written in",hrs,"hrs",mins,"mins and",secs,"seconds"))
  
  ### --- Generate results root filenames
  if (is.null(results.path)) { # Will write to a temporary directory
    if (!is.null(debug.path)) {
      results.path <- paste(debug.path,"/",results.label,"_Results_",now, sep="")
    } else {
      results.path <- paste(tempdir(),"/",results.label,"_Results_",now, sep="")      
      clean.up.results <- TRUE
    }
  }
  try(system(paste("mkdir '",results.path,"'", sep="")))
  results.root <- paste(results.path, results.label, sep="/")
  results.file <- paste(results.root,".txt",sep="")
  plot.file <- paste(results.root,".pdf",sep="")    
  
  ### --- Generate commands
  n.thin <- max(100*n.mil,1)
  n.its <- n.mil*1e6
  n.output <- round(n.its/10)
  comm <- paste(
    "java -jar \"", bayesglm.jar, "\" \"",
    arguments.file, "\" \"", data.file, "\" \"",
    results.file, "\" ",
    format(n.its,sci=F)," ",0," ",
    format(n.thin,sci=F)," ",
    format(n.output,sci=F)," ",
    seed," ",
    as.integer(model.selection)," ",
    as.integer(g.prior)," ",
    as.integer(alt.initial.values)," ",
    as.integer(max.model.dim)," ",
    n.mod.space.comps," ",
    modSpaceDistribution, sep=""
  )
  for (i in 1:n.mod.space.comps) {
    if (modSpaceDistribution==0) {
      comm <- paste(comm, format(model.space.priors[[i]][["Rate"]],sci=F) )      
    } else if (modSpaceDistribution==1) {
      comm <- paste(comm, format(model.space.priors[[i]][["a"]],sci=F) )      
      comm <- paste(comm, format(model.space.priors[[i]][["b"]],sci=F) )      
    }
  }
  if (n.mod.space.comps>1) {
    for (i in 1:(n.mod.space.comps-1) ) {
      comm <- paste(comm, length(model.space.priors[[i]][["Variables"]]) )
    }
  }
  
  ### --- Run commands
  cat("\nCalling BGLiMS...\n\n")
  cat("\nCommand: \n")
  cat(comm, "\n")
  jobs <- list()  
  t1 <- proc.time()["elapsed"]
  system(comm)
  t2 <- proc.time()["elapsed"]
  bglims.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("\nBGLiMS analysis finished in",hrs,"hrs",mins,"mins and",secs,"seconds"))
  
  ### --- Results processing
  cat("\nProcessing BGLiMS output...")
  t1 <- proc.time()["elapsed"]
  results <- .ReadResults(
    results.file=results.file,
    first.n.mil.its=NULL)
  if (do.chain.plot) {ChainPlots(results, plot.file=plot.file)}
  t2 <- proc.time()["elapsed"]
  results.processing.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("\nResults processed in",hrs,"hrs",mins,"mins and",secs,"seconds"))
  
  # Clean up
  cat("\nCleaning up...")
  if(clean.up.arguments) { system(paste("rm -rf ",arguments.path, sep="")) }
  if(clean.up.data) { system(paste("rm -rf ",data.path, sep="")) }
  if(clean.up.results) { system(paste("rm -rf ",results.path, sep="")) } else {
    results[["raw.results.file"]] <- results.file
  }
  
  # Append extra information to results object
  results$args$model.selection <- model.selection
  results$args$mcmc.seed <- seed
  results$args$model.space.priors <- model.space.priors
  results$args$confounders <- confounders
  results$args$predictors <- predictors
  results$run.times <- list()
  results$run.times$write.time <- write.time
  results$run.times$bglims.time <- bglims.time
  results$run.times$results.processing.time <- results.processing.time
  if (enumerate.up.to.dim>0) {
    results$approx.probs <- EnumeratedApproxPostProbs(results)
  }
  
  # Append rsults table to results object
  results$results.table <- ResultsTable(results)
  
  return(results)
}
