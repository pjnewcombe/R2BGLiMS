#' @include HiddenFunctions.R
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
#' "JAM" (for conjugate linear regression using summary statistics, integrating out parameters) 
#' and "JAM_MCMC" (for linear regression using summary statistics, with full MCMC of all parameters).
#' @param data Matrix or dataframe containing the data to analyse. 
#' Rows are indiviuals, and columns contain the variables and outcome.
#' If modelling summary statistics specify X.ref, marginal.betas, and n instead (see below).
#' @param outcome.var Name of outcome variable in data. For survival data see times.var below.
#' If modelling summary statistics with JAM this can be left null but you must specify X.ref, marginal.beats and n instead (see below).
#' @param times.var SURVIVAL DATA ONLY Name of column in data which contains the event times.
#' @param confounders Optional vector of confounders to fix in the model at all times, i.e. exclude from model selection.
#' @param model.selection Whether to use model selection (default is TRUE). NB: Even if set to FALSE, please provide a 
#' dummy model.space.priors argument (see below). This will not be used mathematically, but the package requires taking the list of variable
#' names from the "Variables" element therein.
#' @param model.space.priors Must be specified if model.selection is set to TRUE.
#' Two options are available. 1) A fixed prior is placed on the proportion of causal 
#' covariates, and all models with the same number of covariates is equally likely. This
#' is effectively a Poisson prior over the different possible model sizes. A list must
#' be supplied for `model.space.priors` with an element "Rate", specifying the prior 
#' proportion of causal covariates, and an element "Variables" containing the list of covariates
#' included in the model search. 2) The prior proportion of causal covariates
#' is treated as unknown and given a Beta(a, b) hyper-prior, in which case 
#' elements "a" and "b" must be included in the `model.space.priors` list rather
#' than "Rate". Higher values of "b" relative to "a" will encourage sparsity.
#' NOTE: It is easy to specify different model space priors for different collections of
#' covariates by providing a list of lists, each element of which is itself a model.space.prior
#' list asm described above for a particular subset of the covariates.
#' @param beta.priors There are three options: 
#' 1) Leave as null, in which case N(0, 1e6) priors are placed on the effects
#' of confounders and any predictors included in model selection are assigned a
#' common normal prior with unknown variance.
#' 2) Provide fixed priors for the confounders only. 
#' 3) Provide fixed priors for all covariates. 
#' beta.priors (a matrix/data.frame) can be used to provide fixed priors; rows must be named with the corresponding 
#' variable names in the data, and include Guassian prior means and sds in the first and 
#' second columns respectively.
#' @param g.prior Whether to use a g-prior for the beta's - i.e. a multivariate normal 
#' with correlation structure proportional to sigma^2*X'X^-1 (set to TRUE) or to use independence priors (leave as FALSE).
#' @param model.tau Whether to model tau or not (default FALSE). If set to true,
#' then a Zellner-Siow prior is used, centred on the value provide by tau. The value provided in tau is also used as the
#' initial value.
#' @param tau Value to use for sparsity parameter tau (under the tau*sigma^2 parameterisation).
#' A recommended default is max(n, P^2) where n is the number of individuals, and P the number of predictors. 
#' @param enumerate.up.to.dim Whether to make posterior inference by exhaustively calculating
#' the posterior support for every possible model up to this dimension. Leaving at 0 to disable
#' and use RJMCMC instead. The current maximum allowed value is 5.
#' @param X.ref: Reference genotype matrix used by JAM to impute the SNP-SNP correlations. If multiple regions are to be 
#' analysed this should be a list containing reference genotype matrices for each region. Individual's genotype
#' must be coded as a numeric risk allele count 0/1/2. Non-integer values reflecting imputation may be given.
#' NB: The risk allele coding MUST correspond to that used in marginal.betas. These matrices must each be positive definite and
#' the column names must correspond to the names of the marginal.betas vector.
#' @param marginal.betas Vector of marginal effect estimates to re-analyse with JAM under multivariate models.
#' @param n: Sample size which marginal.betas were calculated in.
#' @param max.fpr ROC AUC ONLY: Maximum acceptable false positive rate (or x-axis value) to optimise a truncated ROC AUC.
#' @param min.tpr ROC AUC ONLY: Minimum acceptable true positive rate, i.e. sensitivity (or y-axis value) to optimise a truncated ROC AUC.
#' @param n.mil Number of million iterations to run (default is 1)
#' @param seed Which random number seed to use in the RJMCMC sampler.
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
#' @return A R2BGLiMS_Results class object is returned. See the slot 'posterior.summary.table' for a posterior
#' summary of all parameters.
#' 
#' The function \code{\link{PrettyResultsTable}} can be used to print summary posterior results for all parameters. Other functions
#' for summarising results are listed under "see also".
#' 
#' @seealso For a nicer summary of all covariates see 
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
  tau=NULL,
  enumerate.up.to.dim=0,
  X.ref=NULL,
  marginal.betas=NULL,
  n=NULL,
  max.fpr=1,
  min.tpr=0,
  n.mil=1,
  seed=1,
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
  alt.initial.values <- FALSE # Now done using the extra.arguments option
  
  ##############################
  ##############################
  ### --- Error messages --- ###
  ##############################
  ##############################
  
  ### --- Java installation error messages
  try.java <- try(system("java -version"), silent=TRUE)
  if (try.java!=0) stop("Java is not installed and is required to run BGLiMS.\nPlease install a Java JDK from java.com/download.")  
  
  ### --- Basic input checks
  if (is.null(likelihood)) stop("No likelihood, i.e. the model type, has been specified; please specify as Logistic,
                                Weibull, Cox, Gaussian, GaussianConj, RocAUC, JAM_MCMC or JAM")
  if (!is.null(likelihood)) {
    if (!likelihood %in% c("Logistic", "Weibull", "Cox", "Gaussian", "GaussianConj", "JAM_MCMC", "JAM", "RocAUC", "RocAUC_Testing")) {
      stop("Likelihood must be specified as Logistic, Weibull, Cox, Gaussian, GaussianConj, RocAUC, JAM_MCMC, or JAM")
    }
  }
  if (is.null(data)&is.null(X.ref)) stop("The data to analyse has not been specified")
  if (is.null(outcome.var)&is.null(marginal.betas)) stop("An outcome variable has not been specified")
  
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
  if (likelihood %in% c("GaussianConj") & is.null(tau)) {
    cat("tau was not provided, setting to the maximum of n and P^2\n")
    tau <- max(nrow(data), ncol(data)^2)
  }
  if (likelihood %in% c("JAM") & is.null(tau)) {
    cat("tau was not provided, setting to P^2\n")
    tau <- length(marginal.betas)^2
  }
  
  ### --- Enumeration error messages
  if (enumerate.up.to.dim>0) {
    if (model.tau) stop ("Tau must be fixed to enumerate model specific posterior
                         scores.")
    if ((enumerate.up.to.dim>5)) stop ("Currenly only possible to enumerate models
                                     up to dimension 5") # If change edit help above    
  }  
  
  ### --- Model space prior error messages
  if (is.null(model.space.priors)) stop("Must specify a prior over the model space.")
  if (!is.null(model.space.priors)) {
    if (!is.list(model.space.priors)) stop("model.space.priors must be a list, or list of list(s).")
    if (!is.list(model.space.priors[[1]])) {
      # Convert to list of lists if only a single model space component is provided
      model.space.priors <- list(model.space.priors)
    }
    # Check structure
    for (c in 1:length(model.space.priors)) {
      no.prior.family <- TRUE
      if ("a"%in%names(model.space.priors[[c]])&"b"%in%names(model.space.priors[[c]])) {
        no.prior.family <- FALSE
        model.space.priors[[c]]$Family <- "Poisson"
      }
      if ("Rate"%in%names(model.space.priors[[c]])) {
        no.prior.family <- FALSE        
        model.space.priors[[c]]$Family <- "BetaBinomial"
      }
      if (no.prior.family) {
        stop("Each model.space.prior list(s) must contain either named a and b elements to sepecify a beta-binomial prior, or a named Rate element to specify a Poisson prior")
      }
      if (!"Variables"%in%names(model.space.priors[[c]])) {
        stop("Each model.space.prior list(s) must contain an element named Variables, defining which covariates to search over")
      }
      if (!likelihood %in% c("JAM_MCMC", "JAM")) {
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
    if (likelihood %in% c("GaussianConj", "JAM")) {stop("Fixed priors for the coefficients can not be specified for the conjugate model; the prior
                                                 is take as a function of X'X")}
    beta.priors.not.mat <- TRUE
    if (is.data.frame(beta.priors)) {beta.priors.not.mat <- FALSE}
    if (is.matrix(beta.priors)) {beta.priors.not.mat <- FALSE}
    if (beta.priors.not.mat) stop("beta.priors must be a matrix or data frame")
    if (ncol(beta.priors)!=2) stop("beta.priors must have two columns - 1st for means, 2nd for SDs")
    if ( is.null(rownames(beta.priors)) ) stop("Rows of beta.priors must be named with corresponding variable names")
    if (!likelihood %in% c("JAM_MCMC", "JAM")) {
      if (sum(rownames(beta.priors)%in%colnames(data))!=nrow(beta.priors)) stop("One or more variables in beta.priors are not present in the data")
    }
  }  
  
  ### --- X.ref error message for JAM
  if (!is.null(X.ref)) {
    if (is.data.frame(X.ref)) {
      X.ref <- matrix(X.ref) # convert to matrix
    }
    if (!is.list(X.ref)) {
      X.ref <- list(X.ref) # convert to list if there is a single block
    }
    if (is.null(marginal.betas)) { stop("For analysis with JAM you must provide a vector of marginal summary statistics") }
    if (is.null(n)) { stop("You must specificy the number of individuals the marginal effect estimates were calculated in.") }
    if (sum(unlist(lapply(X.ref, function(x) !is.numeric(x) )))>0) {stop("Reference genotype matrices must be numeric, coded as risk allele countsin the 0 to 2 range")}
    if (max(unlist(X.ref))>2 | min(unlist(X.ref)) < 0) {stop("Reference genotype matrices must be coded coded as risk allele counts in the 0 to 2 range")}
    if (sum(names(marginal.betas) %in% unlist(lapply(X.ref, colnames))) < length(marginal.betas)) {stop("Reference genotype matrices do not include all SNPs in the marginal.betas vector")}
  }
  
  # Setup file paths/sytem command information
  if (.Platform$OS.type == "windows") {
    fsep <- "\\"
    del.command <- "rmdir /s /q"
  } else {
    fsep <- "/"
    del.command <- "rm -rf"
  }  
  pack.root <- path.package("R2BGLiMS")
  bayesglm.jar <- file.path(pack.root, "BGLiMS", "BGLiMS.jar", fsep=fsep)
  if (!is.null(debug.path)) {
    debug.path <- file.path(debug.path, fsep=fsep)
    main.path <- debug.path
    clean.up.data <- FALSE
    clean.up.results <- FALSE
    clean.up.arguments <- FALSE
  } else {
    main.path <- tempdir()
    clean.up.data <- TRUE
    clean.up.results <- TRUE
    clean.up.arguments <- TRUE
  }
  if (!is.null(results.path)) { # Must be run after the above
    results.path <- file.path(results.path, fsep=fsep)    
    clean.up.results <- FALSE
  }
  
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
  t1 <- proc.time()["elapsed"]
  now <-format(Sys.time(), "%b%d%H%M%S") # Used to ensure unique names in the temp directory
  # Set arguments
  load(file.path(pack.root, "data", "DefaultArguments.rda", fsep=fsep))
  if (!is.null(extra.arguments)) {
    for (arg in names(extra.arguments)) {
      cat("Setting user specified argument ",arg,"\n")    
      default.arguments[[arg]] <- extra.arguments[[arg]]
    }    
  }
  # Setup arguments file
  arguments.path <- file.path(main.path, paste(results.label,"_Arguments_",now, sep=""), fsep=fsep)
  try(system(paste("mkdir '",arguments.path,"'", sep="")))
  arguments.file <- file.path(arguments.path, paste(results.label,"_Arguments.txt",sep=""), fsep=fsep)
  # Write arguments
  write(paste(names(default.arguments)[1],default.arguments[[1]]), file = arguments.file)    
  for (arg in names(default.arguments)[-1]) {
    write(paste(arg,format(default.arguments[[arg]],sci=F)), file = arguments.file, append = T)    
  }
  rm(default.arguments)
  
  #########################
  #########################
  #########################
  ### --- JAM Setup --- ###
  #########################
  #########################
  #########################
  
  if (likelihood %in% c("JAM", "JAM_MCMC")) {
    ### --- Generate X'X, after normalising X
    xTx <- list()
    for (ld.block in 1:length(X.ref)) {
      # Normalise X
      X.normalised <- apply(X.ref[[ld.block]], MAR=2, function(x) x-mean(x))
      # Calculate X'X
      xTx[[ld.block]] <- t(X.normalised) %*% X.normalised
    }
    
    ### --- Generate Xy for JAM
    z <- rep(NA,length(marginal.betas))
    names(z) <- names(marginal.betas)
    for (ld.block in 1:length(X.ref)) {
      mafs <- apply(X.ref[[ld.block]], MAR=2, mean)/2 # Take MAFs from reference X
      for (snp in colnames(X.ref[[ld.block]])) {
        # Group counts under HWE
        n1 <- n*(1-mafs[snp])*mafs[snp]*2
        n2 <- n*mafs[snp]*mafs[snp]
        # Group means from beta-hat (mean-centred)
        y0 <- -(n1 * marginal.betas[snp] + n2 * 2 * marginal.betas[snp])/n
        y1 <- y0 + marginal.betas[snp]
        y2 <- y0 + 2*marginal.betas[snp]
        z[snp] <- y1 * n1 + 2 * y2 * n2
      }
    }    
  }

  ### --- Write data
  data.file <- NULL
  if (is.null(data.file)) {
    cat("\nWriting temporary data files...\n")
    data.path <- file.path(main.path, paste(results.label,"_Data_",now, sep=""), fsep=fsep)      
    system(paste("mkdir '",data.path,"'", sep=""))
    data.file <- file.path(data.path, paste(results.label,".txt",sep=""), fsep=fsep)
    .WriteData(
      data.file=data.file,
      likelihood=likelihood,
      data=data,
      outcome.var=outcome.var,
      times.var=times.var,
      confounders=confounders, 
      predictors=predictors,
      model.space.priors=model.space.priors,
      beta.priors=beta.priors,
      g.prior=g.prior,
      model.tau=model.tau,
      tau=tau,
      enumerate.up.to.dim=enumerate.up.to.dim,
      xTx=xTx,
      z=z,
      max.fpr=max.fpr,
      min.tpr=min.tpr,
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
    results.path <- file.path(main.path, paste(results.label,"_Results_",now,sep=""), fsep=fsep)
  }
  try(system(paste("mkdir '",results.path,"'", sep="")))
  results.file <- file.path(results.path, paste(results.label,".txt",sep=""), fsep=fsep)
  plot.file <- file.path(results.path, paste(results.label,".pdf",sep=""), fsep=fsep)
  
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
  
  ##################################
  ##################################
  ##################################
  ### --- Results processing --- ###
  ##################################
  ##################################
  ##################################
  
  cat("\n\nReading in the BGLiMS results file...\n")  
  t1 <- proc.time()["elapsed"]
  
  ### --- Read BGLiMS arguments
  bglims.arguments <- as.list(read.table(results.file, header=T, sep=" ", nrows=1))
  bglims.arguments$Likelihood <- as.character(bglims.arguments$Likelihood)
  bglims.arguments$ModelSpacePriorFamily <- as.character(bglims.arguments$ModelSpacePriorFamily)  

  n.lines.until.rjmcmc.output <- 3 # There are always three lines of meta-data
  enumerated.posterior.inference <- list("No enumeration was done.")
  if (enumerate.up.to.dim>0) {
    ### -- Read enumerated posterior scores
    enumerated.model.posterior.scores <- NULL
    for (max.model.dimension in c(0:enumerate.up.to.dim)) { # From 0 to account for the null model
      n.models.this.dimension <- choose(bglims.arguments$V - bglims.arguments$startRJ, max.model.dimension)
      model.scores.this.dimension <- read.table(
        results.file,
        skip = n.lines.until.rjmcmc.output,
        header=FALSE,
        nrows=n.models.this.dimension)    
      n.lines.until.rjmcmc.output <- n.lines.until.rjmcmc.output+n.models.this.dimension
      enumerated.model.posterior.scores <- rbind(enumerated.model.posterior.scores, model.scores.this.dimension)    
    }
    colnames(enumerated.model.posterior.scores) <- c("Model", "PosteriorScore")
    enumerated.model.posterior.scores$Model <- as.character(enumerated.model.posterior.scores$Model)
    enumerated.posterior.inference <- .ApproxPostProbsByModelEnumeration(enumerated.model.posterior.scores, model.space.priors, enumerate.up.to.dim)
  }
  
  ### --- Read MCMC output
  n.rows.written <- bglims.arguments$iterations/bglims.arguments$thin
  bglims.rjmcmc.output <- read.table(
    results.file,
    skip = n.lines.until.rjmcmc.output,
    header=TRUE,
    nrows=n.rows.written)
  Lhalf <- round(nrow(bglims.rjmcmc.output)/2)     	 # Burnin is a half	
  bglims.rjmcmc.output <- bglims.rjmcmc.output[(Lhalf+1):nrow(bglims.rjmcmc.output),]   # Remove burnin

  cat("... finished reading BGLiMS results.\n")  
  
  ### --- Summary table
  posterior.summary.table <- matrix(NA,ncol(bglims.rjmcmc.output),8)
  colnames(posterior.summary.table) = c("PostProb","Median","CrI_Lower","CrI_Upper",
    "Median_Present","CrI_Lower_Present","CrI_Upper_Present","BF")
  rownames(posterior.summary.table) <- colnames(bglims.rjmcmc.output)
  # Prior probabilties - used for Bayes Factors below
  prior.probs <- rep(NA, nrow(posterior.summary.table))
  names(prior.probs) <- rownames(posterior.summary.table)
  for (c in 1:length(model.space.priors)) {
    if ("Rate" %in% names(model.space.priors[[c]]) ) {
      prior.probs[model.space.priors[[c]]$Variables] <- .ModelSpaceSpecProb(length(model.space.priors[[c]]$Variables), model.space.priors[[c]]$Rate)
    } else {
      prior.probs[model.space.priors[[c]]$Variables] <- model.space.priors[[c]]$a/(model.space.priors[[c]]$a + model.space.priors[[c]]$b)                
    }
  }
  # Fill in the summary table
  for (v in rownames(posterior.summary.table)) {
    posterior.summary.table[v,c("CrI_Lower", "Median", "CrI_Upper")] <- quantile(bglims.rjmcmc.output[,v],c(0.025, 0.5, 0.975))
    posterior.summary.table[v,c("CrI_Lower_Present", "Median_Present", "CrI_Upper_Present")] <- quantile(bglims.rjmcmc.output[,v][bglims.rjmcmc.output[,v]!=0],c(0.025, 0.5, 0.975) )
    if (v %in% unlist(lapply(model.space.priors, function(x) x$Variables))) {
      posterior.summary.table[v,"PostProb"] <- length( bglims.rjmcmc.output[,v][bglims.rjmcmc.output[,v]!=0] ) / nrow(bglims.rjmcmc.output)      
      if (enumerate.up.to.dim>0) {
        # Replace with enumeration probs
        posterior.summary.table[v,"PostProb"] <- enumerated.posterior.inference$marg.probs[v]
      }
      posterior.summary.table[v,"BF"] <- .BayesFactor(prior.probs[v], posterior.summary.table[v,"PostProb"])
      if (likelihood %in% c("Weibull", "Cox", "Logistic", "RocAUC", "RocAUC_Testing") ) {
        # Exponentiate quantiles
        posterior.summary.table[v,c("CrI_Lower", "Median", "CrI_Upper","CrI_Lower_Present", "Median_Present","CrI_Upper_Present")] <- exp(posterior.summary.table[v,c("CrI_Lower", "Median", "CrI_Upper","CrI_Lower_Present", "Median_Present","CrI_Upper_Present")])
      }
    }
  }
  
  t2 <- proc.time()["elapsed"]  
  results.processing.time <- t2-t1
  
  ########################
  ########################
  ### --- Clean up --- ###
  ########################
  ########################
  
  cat("\nCleaning up...")
  if(clean.up.arguments) { system(paste(del.command,arguments.path)) }
  if(clean.up.data) { system(paste(del.command,data.path)) }
  if(clean.up.results) { system(paste(del.command,results.path)) }    
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("\nResults processed in",hrs,"hrs",mins,"mins and",secs,"seconds"))
  
  #########################
  #########################
  ### --- Run times --- ###
  #########################
  #########################
  
  run.times <- list()
  run.times$write.time <- write.time
  run.times$bglims.time <- bglims.time
  run.times$results.processing.time <- results.processing.time
  
  ################################
  ################################
  ### --- Create S4 Object --- ###
  ################################
  ################################
  if (is.null(confounders)) {
    confounders <- c("None")
  }
  results <- new(
    "R2BGLiMS_Results",
    likelihood = likelihood,
    posterior.summary.table = posterior.summary.table,
    enumerate.up.to.dim = enumerate.up.to.dim,
    enumerated.posterior.inference = enumerated.posterior.inference,
    n.iterations = n.its,
    thin = bglims.arguments$thin,
    model.space.priors = model.space.priors,
    confounders = confounders,
    run.times=run.times,
    n.covariate.blocks.for.jam = 1,
    bglims.arguments=bglims.arguments,
    bglims.rjmcmc.output=bglims.rjmcmc.output
    )
  
  return(results)
}
