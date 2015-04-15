#' Calls BGLiMS - a Java package for fitting GLMs under Bayesian model selection. NOTE: The predictors to explore with model
#' selection are specified via the model.space.priors argument - see examples. By Default a common, and unknown, prior SD
#' is used for all predictors under model selection, which in turn is assigned an Inverse-Gamma hyper-prior.
#' Fixed normal priors may be user specified for all predictor (and confounder) coefficients via the beta.priors argument.
#' 
#' @export
#' @title Call BGLiMS from R
#' @name R2BGLiMS
#' @param likelihood Type of model to fit. Current options are "Logistic" (for binary data), "Weibull" (for survival data), 
#' "Gaussian" (for continuous data), "GaussianConj" (linear regression exploiting conjugate results), "GaussianMarg" (for analysis of univariate associations from Gaussian linear 
#' regressions) and "GaussianMargConj" (for analysis under a marginal conjugate linear regression model).
#' @param data Matrix or dataframe containing the data to analyse. Rows are indiviuals, and columns are variables, named in concordance
#' with the following options.
#' @param outcome.var Which column in data contains the binary outcome/survival status variable (default "Disease")
#' @param times.var If survival data, the column in data which contains the event times (default NULL)
#' @param xTx GaussianMarg and GaussianMargConj ONLY: List containing each block's plug-in estimate for X'X.
#' @param t GaussianMarg and GaussianMargConj ONLY: Vector of quantities calculated from the summary statistics.
#' @param sigma2_invGamma_a Guassian models ONLY: Inverse-Gamma parameter one for the residual
#' precision. For the conjugate model this parameter is intergrated out, and this may be provided as in the
#' Bottolo and Richardson 2010 notation. For an informative prior for the non-conjugate model
#' (taking into account Java parameterisation) choose N/2.
#' @param sigma2_invGamma_b Guassian models ONLY: Inverse-Gamma parameter one for the residual
#' precision. For the conjugate model this parameter is intergrated out, and this may be provided as in the
#' Bottolo and Richardson 2010 notation. For an informative prior for the non-conjugate model
#' (taking into account Java parameterisation) choose N/(2*variance estimate).
#' @param g.prior Gaussian conjugate models ONLY: Whether to use a g-prior for the beta's - i.e. a multivariate normal 
#' with correlation structure proportional to sigma^2*X'X^-1 or to use independence priors (default = FALSE).
#' @param tau Gaussian conjugate models ONLY: Value to use for sparsity parameter tau (tau*sigma^2 parameterisation).
#' Default residual.var.n. If modelling this parameter, this value is used to center the Zellner-Siow prior
#' and as an intial value.
#' @param model.tau Gaussian conjugate models ONLY: Whether to model tau or not (default FALSE). If set to true,
#' then a Zellner-Siow prior is used, centred on the value provide by tau. The value provided in tau is also used as the
#' initial value.
#' @param tau.proposal.sd Gaussian conjugate models ONLY: When modelling tau an initial SD to use in the adaption
#' of the proposal distribution. Under the conjugate model, tau makes use of the BGLiMS parameter `betaPriorSd'. Thus
#' the parameters is used on a very different scale to normal and as such it is best to specify this seperately.
#' Defaults to 0.05.
#' @param all.model.scores.up.to.dim Gaussian conjugate models ONLY: When NOT modelling tau, whether to output the posterior scores
#' for every possible model of dimension up to the integer specified. Currenly maximum allowed dimension is 2. Setting
#' to default 0 means this is not carried out. Can be used as an alternative to running the RJMCMC.
#' @param cluster.var If hierarchical data and random intercepts are required, the column in data contains the clustering variable (default NULL)
#' @param confounders vector of confounders to fix in the model at all times, i.e. exclude from model selection (default NULL)
#' @param model.selection Whether to use model selection (default is TRUE). NB: Even if set to FALSE, please provide a 
#' dummy model.space.priors argument (see below). This will not be used mathematically, but the package requires taking the list of variable
#' names from the "Variables" element therein.
#' @param model.space.priors Must be specified if model.selection is set to TRUE.
#' A list with as many elements as desired model space components (one or more).
#' Each element of the list is also a list, corresponding to a particular
#' model space component listing the prior parameters and covariates in that
#' component. Either a Poisson prior on model space can be used, in which case an
#' element "Rate" must be included, or Beta-Binomial model space priors can
#' be used, in which case elements "a" and "b" must be included corresponding
#' to the beta-binomial hyperparameters. Finally, an element "Variables" must
#' be included - a vector of variable names in that component. Note that when choosing a
#' beta-binomial model space prior, higher values of "b" relative to "a" increase sparsity,
#' whereas higher values of "a" encourage the inclusion of more covariates. a ~ number of prior successes, b~number of prior
#' failures.
#' @param max.model.dim Optional specification of maximum model dimension (default -1 means no maximum is set)
#' @param initial.model Optionally, an initial model can be provided as a vector of 0's and 1's. Default is NULL
#' and the null model is used. If set to 1, the saturated model is used.
#' @param beta.priors Optional matrix or data.frame containing two columns (1st:mean, 2nd:SD) describing normal priors for
#' the beta coefficients of the confounders or confounders and predictors (i.e. log-ORs in the case of logistic regression).
#' Rows of
#' beta.priors must be named with the corresponding variable names in the data. NOTE: For predictors included in the model
#' selection search, the interpretation of these priors is on the effect conditional on inclusion in the model.
#' There are three options: 1) Do not specify in which case N(0, 1e6) priors are assigned to the betas of any confounders, 
#' and the effects of any predictors included in model selection are assigned a common normal prior with unknown variance,
#' which in turn assigned a InversGamma(0.001, 0.001) hyperprior. 2) Provide a beta.priors matrix including confounders only. 
#' This overrides the default use of N(0,1e6) confounder priors, but the predictors for which model selection is performed for
#' are still assigned a common prior with unknown variance. 3) provide a beta.priors matrix including all confounders and predictors,
#' thereby avoiding the use of a common prior with unkown variance across predictors for which model selection is performed for.
#' @param n.mil Number of million iterations to run (default is 1)
#' @param seed Which random number seed to use in the RJMCMC sampler (default is 1)
#' @param alt.initial.values Whether to use the alternative set of initial values in the arguments file (default is FALSE)
#' @param results.path Optional path if wish to to save algorithm results output (when set to the default NULL 
#' all output is written to a temporary directory and deleted). NOTE this must exist already. A folder will be created in this
#' path taking the name of `results.label' below
#' @param results.label Optional label to for algorithm output files (if you have specified results.path)
#' @param args.file Optional path to an alternative arguments file - e.g. "/Users/pauln/Arguments/Arguments_Alt.txt". 
#' If left as NULL, a default file located at BGLiMS/Arguments_Package.txt within the package directory is used. Note that
#' this file can be used as a template. Currently this is not documented - please contact the package author, Paul Newcombe.
#' @param do.chain.plot Whether to produce a PDF containing chain plots. Default FALSE.
#' @param debug.path Optional path to save the data and results files to (rather than as temporary files), for aid in
#' debugging.
#' 
#' @return A Reversible Jump results object is returned. This is a list of two elements: "args" which records various modelling
#' arguments used in the analysis, and "results" - a matrix containing all saved posterior samples from the analysis. Columns
#' correspond to parameters and each row contains values from a particular itertation, with 0 indicating exclusion from the model.
#' 
#' The function \code{\link{ResultsTable}} can be used to print summary posterior results for all parameters. Other functions
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
  xTx=NULL,
  t=NULL,
  sigma2_invGamma_a=NULL,
  sigma2_invGamma_b=NULL,
  g.prior=FALSE,
  tau=NULL,
  model.tau=FALSE,
  tau.proposal.sd=0.05,
  all.model.scores.up.to.dim=0,
  cluster.var=NULL,
  confounders=NULL,
  model.selection=TRUE,
  model.space.priors=NULL,
  max.model.dim=-1,
  initial.model=NULL,
  beta.priors=NULL,
  n.mil=1,
  seed=1,
  alt.initial.values=FALSE,
  results.path=NULL,
  results.label="R2BGLiMS",
  args.file=NULL,
  do.chain.plot=FALSE,
  debug.path=NULL
) {
  ### --- Error messages
  try.java <- try(system("java -version"), silent=TRUE)
  if (try.java!=0) stop("Java is not installed and is required to run BGLiMS.\nPlease install from java.com/download.")  
  if (is.null(likelihood)) stop("No likelihood, i.e. the model type, has been specified; please specify as Logistic,
                                Weibull, Gaussian, GaussianConj, GaussianMarg or GaussianMargConj")
  if (!is.null(likelihood)) {
    if (!likelihood %in% c("Logistic", "Weibull", "Gaussian", "GaussianConj", "GaussianMarg", "GaussianMargConj")) {
      stop("Likelihood must be specified as Logistic, Weibull, Gaussian, GaussianConj, GaussianMarg, or GaussianMargConj")
    }
  }
  if (is.null(data)&is.null(xTx)) stop("The data to analyse has not been specified")
  if (is.null(outcome.var)&is.null(t)) stop("A binary outcome/survival status variable has not been specified")
  if (likelihood %in% c("Logistic", "Weibull")) {
    if (is.factor(data[,outcome.var])) {
      data[,outcome.var] <- as.integer(data[,outcome.var])-1
    } else if (is.character(data[,outcome.var])) {
      data[,outcome.var] <- as.integer(as.factor(data[,outcome.var]))-1    
    }    
    if (length(table(data[,outcome.var]))!=2) stop("Outcome variable must be binary")    
  }
  if (is.null(model.space.priors)) stop("Must specify a prior over the model space (even if model selection
                                        will not be used this is used to provide the list of covariate
                                        names")
  if (!is.null(confounders)) {
    if (sum(confounders%in%colnames(data))!=length(confounders)) stop("One or more confounders are not present in the data")
    for (v in confounders) {
      if (is.factor(data[,v])) {
        data[,v] <- as.integer(data[,v])
      } else if (is.character(data[,v])) {
        data[,v] <- as.integer(as.factor(data[,v]))        
      }
    }
  }
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
  ### -- Beta priors
  if (!is.null(beta.priors)) {
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
  
  ### --- Possible future options
  # NEVER USE
  # fixed.prior.sd Prior SD that will be used in zero-centred normal priors on the betas when none is given (default = 1000)
  # data.file Location of correctly formatted .txt data file. If left as default NA then the arguments below must be provided.
  # args.file Path to an alternative arguments file, if would prefer not to use the defaults (default NULL and the default arguments file will be used)
  # first.n.mil.its how many million iterations to read in from the beginning. Useful if you would like to explore whether
  # convergence could have been achieved with less iterations (default is NULL)
  fixed.prior.sd <- 1000
  data.file <- NULL
  first.n.mil.its <- NULL
  
  # Setup file paths
  pack.root <- path.package("R2BGLiMS")
  bayesglm.jar <- paste(pack.root,"/BGLiMS/BGLiMS.jar",sep="")
  if (is.null(args.file)) {
    args.file <- paste(pack.root,"/BGLiMS/Arguments_Package.txt",sep="")    
  } else {
    cat("Using alternative arguments file",args.file,"\n")
  }
  
  ### --- Prior setup
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
      cbind(rep(0,length(confounders)),rep(fixed.prior.sd,length(confounders))),
      row.names=confounders)    
  }
  
  ### --- Write data files if not provided
  t1 <- proc.time()["elapsed"]
  clean.up.data <- FALSE
  clean.up.results <- FALSE
  now <-format(Sys.time(), "%b%d%H%M%S") # Used to ensure unique names in the temp directory
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
      t=t,
      sigma2_invGamma_a=sigma2_invGamma_a,
      sigma2_invGamma_b=sigma2_invGamma_b,
      g.prior=g.prior,
      tau=tau,
      model.tau=model.tau,
      tau.proposal.sd=tau.proposal.sd,
      all.model.scores.up.to.dim=all.model.scores.up.to.dim,
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
    system(paste("mkdir '",results.path,"'", sep=""))
  }
  results.root <- paste(results.path, results.label, sep="/")
  results.file <- paste(results.root,".txt",sep="")
  plot.file <- paste(results.root,".pdf",sep="")    
  
  ### --- Generate commands
  n.thin <- max(100*n.mil,1)
  n.its <- n.mil*1e6
  n.output <- round(n.its/10)
  comm <- paste(
    "java -jar \"", bayesglm.jar, "\" \"",
    args.file, "\" \"", data.file, "\" \"",
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
    first.n.mil.its=first.n.mil.its)
  if (do.chain.plot) {ChainPlots(results, plot.file=plot.file)}
  t2 <- proc.time()["elapsed"]
  results.processing.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("\nResults processed in",hrs,"hrs",mins,"mins and",secs,"seconds"))
  
  # Clean up
  cat("\nCleaning up...")
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
  if (all.model.scores.up.to.dim>0) {
    results$approx.probs <- EnumeratedApproxPostProbs(results)
  }
  
  return(results)
}
