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
#' "CLogLog" complementary log-log link for binary data, 
#' "Weibull" (for survival data), 
#' "Cox" (for survival data), 
#' "CaseCohort_Prentice" (for case-cohort survival data with Prentice weighting), 
#' "CaseCohort_Barlow" (for case-cohort survival data with Barlow weighting), 
#' "RocAUC" (to optimise ROC AUC),
#' "RocAUC_Anchoring" (to optimise ROC AUC using the anchoring formulation),
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
#' @param subcohort.var CASE-COHORT DATA ONLY Name of column in data which contains a binary indicator of membership in the
#' propsectviely sampled (and population representative) sub-cohort.
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
#' @param beta.prior.partitions By default, covariates without informative priors are 
#' ascribed a common Gaussian prior with a Unif(0.01,2) hyper-prior on the effect 
#' standard deviation. beta.prior.partitions can be used to partition the covariates 
#' into different prior groups, within each of which exchangeability can be assumed 
#' and a common prior with an independently estimated variance is used. 
#' beta.prior.partitions must be a list with as many elements as desired 
#' partitions. Each element must in turn be a list containing the following named 
#' elements: "Variables" - a list of covariates to be included in this partition, 
#' and "UniformA" and "UniformB" the Uniform hyper parameters. I.e., for a particular
#' partition, for v in "Variables", 
#' beta[v] | sigma_beta ~ N(0, sigma_beta^2) and 
#' sigma_beta ~ Unif("UniformA", "UniformB").
#' @param dirichlet.alphas.for.roc.model This can either be a single number, if you wish the same concentration parameter 
#' to be used for every covariate. Alternatively a vector can be supplied, with an element corresponding to every variable.
#' @param g.prior Whether to use a g-prior for the beta's, i.e. a multivariate normal 
#' with correlation structure proportional to sigma^2*X'X^-1, which is thought to aid
#' variable selection in the presence of strong correlation. By default this is enabled.
#' @param tau Value to use for sparsity parameter tau (under the tau*sigma^2 parameterisation).
#' When using the g-prior, a recommended default is max(n, P^2) where n is the number of individuals, and P is the number of predictors. 
#' @param enumerate.up.to.dim Whether to make posterior inference by exhaustively calculating
#' the posterior support for every possible model up to this dimension. Leaving at 0 to disable
#' and use RJMCMC instead. The current maximum allowed value is 5.
#' @param X.ref Reference genotype matrix used by JAM to impute the SNP-SNP correlations. If multiple regions are to be 
#' analysed this should be a list containing reference genotype matrices for each region. Individual's genotype
#' must be coded as a numeric risk allele count 0/1/2. Non-integer values reflecting imputation may be given.
#' NB: The risk allele coding MUST correspond to that used in marginal.betas. These matrices must each be positive definite and
#' the column names must correspond to the names of the marginal.betas vector.
#' @param yTy.ref Estimate of the trait variance * N; for use with JAMv2 (Daniel Ahfock's extension).
#' @param n.for.jam For JAMv2 (Daniel Ahfock's extension), the N of the data must be provided.
#' @param marginal.betas Vector of marginal effect estimates to re-analyse with JAM under multivariate models.
#' @param n: Sample size which marginal.betas were calculated in.
#' @param subcohort.sampling.fraction: CaseCohort_Barlow ONLY: The sampling fraction of the sub-cohort from the full cohort, in order
#' to calculate weights for use with the Barlow Case-Cohort pseudo-likelihood (Barlow 1994) REF
#' @param casecohort.pseudo.weight: CaseCohort ONLY: Multiplier for the pseudo log-likelihood relative to
#' the prior.
#' @param max.fpr ROC AUC ONLY: Maximum acceptable false positive rate (or x-axis value) to optimise a truncated ROC AUC.
#' @param min.tpr ROC AUC ONLY: Minimum acceptable true positive rate, i.e. sensitivity (or y-axis value) to optimise a truncated ROC AUC.
#' @param n.iter Number of iterations to run (default is 1e6)
#' @param n.mil.iter Number of million iterations to run. Can optionally be used instead of n.iter for convenience, which it will overide if specified.
#' @param thinning.interval Every nth iteration to store (i.e. for the Java algorithm to write to a file and then read into R). By default this is the
#' number of iterations divided by 1e4 (so for 1 million iterations every 100th is stored.) Higher values (so smaller posterior sample) can lead to 
#' faster runtimes for large numbers of covariates.
#' @param seed Which random number seed to use in the RJMCMC sampler. If none is provided, a random seed is picked between 1 and 2^16.
#' @param results.label Optional label for algorithm output files (if you have specified save.path).
#' @param extra.arguments A named list of any additional arguments for BGLiMS. Type "data(DefaultArguments)" and look in the 
#' "default.arguments" list to see the names (which must match) and default values of available extra arguments. Currently these 
#' are not documented - please contact the package author, Paul Newcombe for details.
#' @param initial.model Optionally, an initial model can be provided as a vector of 0's and 1's. Default is NULL
#' and the null model is used. If set to 1, the saturated model is used.
#' @param max.model.dim Optional specification of maximum model dimension (default -1 means no maximum is set).
#' @param save.path Optional path to save BGLiMS's data and results files. These are usually written as temporary files, and deleted
#' after running R2BGLiMS. However, this option might help for debugging.
#' @param extra.java.arguments A character string to be passed through to the java command line. E.g. to specify a
#' different temporary directory by passing "-Djava.io.tmpdir=/Temp".
#' 
#' @return An R2BGLiMS_Results class object is returned. See the slot 'posterior.summary.table' for a posterior summary of all parameters. 
#' See slot 'mcmc.output' for a matrix containing the raw MCMC output from the saved posterior samples (0 indicates a covariate is excluded 
#' from the model in a particular sample. Some functions for summarising results are listed under "see also".
#' 
#' @seealso Summary results are stored in the slot posterior.summary.table. See \code{\link{ManhattanPlot}} for a visual 
#' summary of covariate selection probabilities. For posterior model space summaries see \code{\link{TopModels}}. For 
#' convergence checks see \code{\link{ChainPlots}} and \code{\link{AutocorrelationPlot}}.
#' 
#' @author Paul Newcombe
#' 
#' @example Examples/R2BGLiMS_Examples.R

R2BGLiMS <- function(
  likelihood=NULL,
  data=NULL,
  outcome.var=NULL,
  times.var=NULL,
  subcohort.var=NULL,
  confounders=NULL,
  model.selection=TRUE,
  model.space.priors=NULL,
  beta.priors=NULL,
  beta.prior.partitions=NULL,
  dirichlet.alphas.for.roc.model=0.01,
  g.prior=TRUE,
  tau=NULL,
  enumerate.up.to.dim=0,
  X.ref=NULL,
  yTy.ref=NULL,
  n.for.jam=NULL,
  marginal.betas=NULL,
  subcohort.sampling.fraction=NULL,
  casecohort.pseudo.weight=1,
  n=NULL,
  max.fpr=1,
  min.tpr=0,
  n.iter=1e6,
  n.mil.iter=NULL,
  thinning.interval=NULL,
  seed=NULL,
  extra.arguments=NULL,
  initial.model=NULL,
  max.model.dim=-1,
  results.label=NULL,
  save.path=NULL,
  extra.java.arguments=NULL
) {
  
  ###########################
  ### --- Old options --- ###
  ###########################
  alt.initial.values <- FALSE # Now done using the extra.arguments option
  model.tau <- FALSE # Stochastic modelling of Tau. Too experimental to offer as an option for now.
  
  ##############################
  ##############################
  ### --- Error messages --- ###
  ##############################
  ##############################
  
  ### --- Number of iterations
  if (!is.null(n.mil.iter)) {
    cat("Number of million iterations specified. Overriding n.iter.\n")
    n.iter <- n.mil.iter*1e6
  }
  cat(n.iter,"iterations will be run...\n")
  
  ### --- Java installation error messages
  try.java <- try(system("java -version"), silent=TRUE)
  if (try.java!=0) stop("Java is not installed and is required to run BGLiMS.\nPlease install a Java JDK from java.com/download.")  
  
  ### --- Basic input checks
  if (is.null(likelihood)) stop("No likelihood, i.e. the model type, has been specified; please specify as Logistic, CLogLog,
                                Weibull, Cox, CaseCohort_Prentice, CaseCohort_Barlow, Gaussian, GaussianConj, RocAUC, RocAUC_Anchoring, JAM_MCMC or JAM")
  if (!is.null(likelihood)) {
    if (!likelihood %in% c("Logistic", "CLogLog", "Weibull", "Cox", "CaseCohort_Prentice", "CaseCohort_Barlow", "Gaussian", "GaussianConj", "JAM_MCMC", "JAM", "JAMv2", "RocAUC", "RocAUC_Anchoring")) {
      stop("Likelihood must be specified as Logistic, CLogLog, Weibull, Cox, CaseCohort_Prentice, CaseCohort_Barlow, Gaussian, GaussianConj, RocAUC, JAM_MCMC, or JAM")
    }
  }
  if (is.null(data)&is.null(X.ref)) stop("The data to analyse has not been specified")
  if (is.null(outcome.var)&is.null(marginal.betas)) stop("An outcome variable has not been specified")
  
  ### --- Likelihood specific checks
  if (likelihood %in% c("CaseCohort_Prentice","CaseCohort_Barlow")) {
    # Check sub-cohort covariate given
    if (is.null(subcohort.var)) stop("For the Case-Cohort model must specify which column of data contains the sub-cohort indicator")
    if (length(table(data[,subcohort.var]))!=2) stop("Subcohort membership indicator must be binary")    
  }
  if (likelihood %in% c("CaseCohort_Barlow")) {
    if (is.null(subcohort.sampling.fraction)) stop("For the Barlow Case-Cohort model must specify the subcohort sampling fraction")
  }
  if (likelihood %in% c("Logistic", "CLogLog", "Weibull", "RocAUC", "RocAUC_Anchoring")) {
    # This check is not done for Weibull, Cox or CaseCohort - since with uncensored data the outcome is not binary
    if (is.factor(data[,outcome.var])) {
      data[,outcome.var] <- as.integer(data[,outcome.var])-1
    } else if (is.character(data[,outcome.var])) {
      data[,outcome.var] <- as.integer(as.factor(data[,outcome.var]))-1    
    }    
    if (length(table(data[,outcome.var]))>2) {
      stop("Outcome variable must be binary")
    } else if ((likelihood != "Weibull") & (length(table(data[,outcome.var]))==1)) {
      # Single class allowed for survival models (may all be survivors)
      stop("Outcome variable only has one class")
    }
  }
  if (likelihood %in% c("GaussianConj", "JAM", "JAMv2") & is.null(tau)) {
    if (g.prior) {
      if (likelihood %in% c("GaussianConj")) {
        cat("tau was not provided. Since the g-prior is in use, setting to the recommended maximum of n and P^2\n")
        tau <- max(nrow(data), ncol(data)^2)        
      } else if (likelihood %in% c("JAM", "JAMv2")) {
        cat("tau was not provided. Since the g-prior is in use, setting to the recommended maximum of n and P^2\n")
        tau <- max(n,length(marginal.betas)^2)
      }
    } else {
      stop("Please choose a value for tau. Note that you have selected to use independent priors.
           Did you mean to use the g-prior?")
    }
  }
  if (likelihood %in% c("RocAUC")) {
    if (model.selection) {if(is.null(initial.model)){stop("Must specify an inital model for ROC AUC model 
                                                     selection")}}
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
      if (!likelihood %in% c("JAM_MCMC", "JAM", "JAMv2")) {
        if (sum(model.space.priors[[c]]$Variables%in%colnames(data))!=length(model.space.priors[[c]]$Variables)) {
          stop(paste("Not all variables in model space component",c,"are present in the data"))
        }
        # Sort out any factors - MOVE TO DATA FORMATTING?
        for (v in model.space.priors[[c]]$Variables) {
          if (is.factor(data[,v])) {
            data[,v] <- as.integer(data[,v])
          } else if (is.character(data[,v])) {
            data[,v] <- as.integer(as.factor(data[,v]))-1
          }
        }        
      }
    }
    # Check partitions are dichotomous
    if (length(unique(unlist(lapply(model.space.priors, function(x) x$Variables))))
        != sum(unlist(lapply(model.space.priors, function(x) length(x$Variables))))) {
      stop("There is overlap in the covariates between some of the model space prior partitions")
    }
    # Check confounders do not appear
    if ( sum(unique(unlist(lapply(model.space.priors, function(x) x$Variables))) %in% confounders) > 0) {
      stop("Some of the confounders also appear in the model space prior.")
    }
  }
  
  ### --- Confounders error messages
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
  
  ### -- Beta prior error messages
  if (!is.null(beta.priors)) {
    if (likelihood %in% c("GaussianConj", "JAM", "JAMv2")) {stop("Fixed priors for the coefficients can not be specified for the conjugate model; the prior
                                                 is take as a function of X'X")}
    beta.priors.not.mat <- TRUE
    if (is.data.frame(beta.priors)) {beta.priors.not.mat <- FALSE}
    if (is.matrix(beta.priors)) {beta.priors.not.mat <- FALSE}
    if (beta.priors.not.mat) stop("beta.priors must be a matrix or data frame")
    if (ncol(beta.priors)!=2) stop("beta.priors must have two columns - 1st for means, 2nd for SDs")
    if ( is.null(rownames(beta.priors)) ) stop("Rows of beta.priors must be named with corresponding variable names")
    if (!likelihood %in% c("JAM_MCMC", "JAM", "JAMv2")) {
      if (sum(rownames(beta.priors)%in%colnames(data))!=nrow(beta.priors)) stop("One or more variables in beta.priors are not present in the data")
    }
  } else if (length(confounders)>0) {
    if (!likelihood %in% c("GaussianConj")) { # Confounders are not modelled for GaussianConj
      warning("beta.priors were not provided for the confounders (which should not be
            treated as exchangeable with covariates under model selection). Setting to N(0,1e6).")
      beta.priors <- data.frame(
        cbind(rep(0,length(confounders)),rep(1000,length(confounders))),
        row.names=confounders)
    }
  }
  
  ### --- Beta prior partition error messages
  if (!is.null(beta.prior.partitions)) {
    if (likelihood %in% c("GaussianConj", "JAM", "JAMv2", "RocAUC")) {
      stop("the beta.prior.partitions option is not available for the conjugate models.")
    }
    if (!is.list(beta.prior.partitions)) {beta.prior.partitions <- list(beta.prior.partitions)}
    if (!is.list(beta.prior.partitions[[1]])) stop("beta.prior.partitions must be a list or list of list(s).")
    # Check structure
    for (c in 1:length(beta.prior.partitions)) {
      # Check Hyper-prior structure
      no.hyperprior.parameters <- TRUE
      if ("UniformA"%in%names(beta.prior.partitions[[c]]) & "UniformB"%in%names(beta.prior.partitions[[c]])) {
        no.hyperprior.parameters <- FALSE
        beta.prior.partitions[[c]]$Family <- "Uniform"
        if (!"Init" %in% names(beta.prior.partitions[[c]])) { # Take as mean of Uniform limits
          beta.prior.partitions[[c]]$Init <- mean(c(beta.prior.partitions[[c]]$UniformA, beta.prior.partitions[[c]]$UniformB))
        } else {
          if (beta.prior.partitions[[c]]$Init<beta.prior.partitions[[c]]$UniformA | beta.prior.partitions[[c]]$Init>beta.prior.partitions[[c]]$UniformB) {
            stop(paste("Initial SD for covariate prior partition",c,"is outside the Uniform range."))
          }
        }
      }
      if ("GammaA"%in%names(beta.prior.partitions[[c]]) &"GammaB"%in%names(beta.prior.partitions[[c]])) {
        no.hyperprior.parameters <- FALSE
        beta.prior.partitions[[c]]$Family <- "Gamma"
        if (!"Init" %in% names(beta.prior.partitions[[c]])) {
          beta.prior.partitions[[c]]$Init <- 0.1
        }
      }
      if (no.hyperprior.parameters) stop(paste("Hyper-parameters UniformA/UniformB or GammaA/GammaB not found for covariate prior partition",c))
      # Check Variables
      if (!"Variables"%in%names(beta.prior.partitions[[c]])) {
        stop(paste("Each covariate prior partition must contain an element named Variables defining the covariates in the partition.
             Not supplied for partition",c) )
      }
      if (!likelihood %in% c("JAM_MCMC", "JAM", "JAMv2")) {
        if (sum(beta.prior.partitions[[c]]$Variables%in%colnames(data))!=length(beta.prior.partitions[[c]]$Variables)) {
          stop(paste("Not all variables in covariate prior partition",c,"are present in the data"))
        }
      }
      if (sum(beta.prior.partitions[[c]]$Variables%in%rownames(beta.priors))>0) {
        stop(paste("Informative priors have also been specified for some of the variables in
                   covariate prior partition",c))
      }
    }
    # Check partitions are dichotomous
    if (length(unique(unlist(lapply(beta.prior.partitions, function(x) x$Variables))))
        != sum(unlist(lapply(beta.prior.partitions, function(x) length(x$Variables))))) {
      stop("There is overlap in the covariates between some of the prior partitions")
    }
    # Check confounders do not appear
    if ( sum(unique(unlist(lapply(beta.prior.partitions, function(x) x$Variables))) %in% confounders) > 0) {
      stop("Some of the confounders also appear in the covariate prior partitions. These should not be treated
           as exchangeable with covariates under model selection.")
    }
  } else {
    if (!likelihood %in% c("GaussianConj", "JAM", "JAMv2", "RocAUC")) {
      all.covariates <- unique(c(unlist(lapply(model.space.priors,function(x) x$Variables)),rownames(beta.priors),confounders))
      # Create default single beta.prior.partitions with Uniform(0.01, 2) - SMMR
      if (is.null(beta.priors)) {
        beta.prior.partitions=list(list("UniformA"=0.01, "UniformB"=2, "Variables"=all.covariates, "Family"="Uniform", "Init"=1))
      } else if (nrow(beta.priors)!=length(all.covariates)) {
        beta.prior.partitions=list(list("UniformA"=0.01, "UniformB"=2,"Variables"=all.covariates[!all.covariates%in%rownames(beta.priors)], "Family"="Uniform", "Init"=1))
      }
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
    if (length(X.ref)==1) {
      #qr.decomp <- qr(X.ref[[1]])
      #if (qr.decomp$rank < ncol(X.ref[[1]])) stop ("The reference matrix is not full rank. 
      #                                      Ideally a larger reference sample should be used, 
      #                                      or you could try pruning correlated SNPs.")
    } else {
      #for (ld.block in 1:length(X.ref)) {
      #  qr.decomp <- qr(X.ref[[ld.block]])
      #  if (qr.decomp$rank < ncol(X.ref[[ld.block]])) stop (
      #    paste("The reference matrix for block",ld.block,"is not full rank.
      #        Ideally a larger reference sample should be used, 
      #        or you could try pruning correlated SNPs.")
      #  )
      #}
    }
    if (is.null(marginal.betas)) { stop("For analysis with JAM you must provide a vector of marginal summary statistics") }
    if (is.null(n)) { stop("You must specificy the number of individuals the marginal effect estimates were calculated in.") }
    if (sum(unlist(lapply(X.ref, function(x) !is.numeric(x) )))>0) {stop("Reference genotype matrices must be numeric, coded as risk allele countsin the 0 to 2 range")}
    if (max(unlist(X.ref))>2 | min(unlist(X.ref)) < 0) {stop("Reference genotype matrices must be coded coded as risk allele counts in the 0 to 2 range")}
    if (sum(names(marginal.betas) %in% unlist(lapply(X.ref, colnames))) < length(marginal.betas)) {stop("Reference genotype matrices do not include all SNPs in the marginal.betas vector")}
  }
  
  ##########################
  ##########################
  ### --- File setup --- ###
  ##########################
  ##########################
  
  if (is.null(seed)) {seed <- sample.int(2^16, size = 1)}
  now <-format(Sys.time(), "%b%d%H%M%S") # Used to ensure unique names in the temp directory
  if (is.null(save.path)) {
    clean.up.bglims.files <- TRUE
    arguments.file <- paste(tempfile(),now,"Arguments.txt",sep="_")
    data.file <- paste(tempfile(),now,"Data.txt",sep="_")
    results.file <- paste(tempfile(),now,"Results.txt",sep="_")
  } else {
    clean.up.bglims.files <- FALSE
    if (is.null(results.label)) {results.label <- now}
    arguments.file <- file.path(save.path,paste(results.label,"Arguments.txt",sep="_"))
    data.file <- file.path(save.path,paste(results.label,"Data.txt",sep="_"))
    results.file <- file.path(save.path,paste(results.label,"Results.txt",sep="_"))
  }
  
  # Setup file paths/sytem command information
  pack.root <- path.package("R2BGLiMS")
  bayesglm.jar <- file.path(pack.root, "BGLiMS", "BGLiMS.jar")
  
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
      stop("Currently informative priors can only be supplied (using beta.priors) for either all confounders or everything.")
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
  # Set arguments
  load(file.path(pack.root, "data", "DefaultArguments.rda"))
  if (!is.null(extra.arguments)) {
    for (arg in names(extra.arguments)) {
      cat("Setting user specified argument",arg,"\n")    
      default.arguments[[arg]] <- extra.arguments[[arg]]
    }    
  }
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
  
  if (likelihood %in% c("JAM", "JAMv2", "JAM_MCMC")) {
    ### --- Generate X'X, after normalising X
    xTx <- list()
    for (ld.block in 1:length(X.ref)) {
      # Normalise X
      X.normalised <- apply(X.ref[[ld.block]], MAR=2, function(x) x-mean(x))
      # Calculate X'X
      xTx[[ld.block]] <- t(X.normalised) %*% X.normalised
    }
    
    ### --- Generate z = X'y for JAM
    z <- NULL
    for (ld.block in 1:length(X.ref)) {
      snps.in.block <- colnames(X.ref[[ld.block]])
      z <- c(z, JAM_PointEstimates(
        marginal.betas = marginal.betas[snps.in.block],
        X.ref=X.ref[[ld.block]],
        n=n, just.get.z=TRUE) )
    }
  }

  ### --- Write data
  .WriteData(
    data.file=data.file,
    likelihood=likelihood,
    data=data,
    outcome.var=outcome.var,
    times.var=times.var,
    subcohort.var=subcohort.var,
    confounders=confounders, 
    predictors=predictors,
    model.space.priors=model.space.priors,
    beta.priors=beta.priors,
    beta.prior.partitions=beta.prior.partitions,
    dirichlet.alphas.for.roc.model=dirichlet.alphas.for.roc.model,
    g.prior=g.prior,
    model.tau=model.tau,
    tau=tau,
    enumerate.up.to.dim=enumerate.up.to.dim,
    xTx=xTx,
    z=z,
    yTy.ref=yTy.ref,
    n.for.jam=n.for.jam,
    subcohort.sampling.fraction=subcohort.sampling.fraction,
    casecohort.pseudo.weight=casecohort.pseudo.weight,
    max.fpr=max.fpr,
    min.tpr=min.tpr,
    initial.model=initial.model
  )  
  t2 <- proc.time()["elapsed"]
  write.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("Data written in",hrs,"hrs",mins,"mins and",secs,"seconds.\n"))
    
  ### --- Generate commands
  if (is.null(thinning.interval)) {
    thinning.interval <- max(n.iter/1e4, 1) # If 1 million iterations, save every 100th
  }
  n.iter.report.output <- round(n.iter/10) # How often to report progress to console
  if (!is.null(extra.java.arguments)) {extra.java.arguments <- paste(extra.java.arguments," ",sep="")}
  comm <- paste(
    "java ",extra.java.arguments,"-jar \"", bayesglm.jar, "\" \"",
    arguments.file, "\" \"", data.file, "\" \"",
    results.file, "\" ",
    format(n.iter,sci=F)," ",0," ",
    format(thinning.interval,sci=F)," ",
    format(n.iter.report.output, sci=F)," ",
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
  cat("Calling BGLiMS Java algorithm with command: \n")
  cat(comm, "\n")
  jobs <- list()  
  t1 <- proc.time()["elapsed"]
  system(comm)
  t2 <- proc.time()["elapsed"]
  bglims.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("BGLiMS Java computation finished in",hrs,"hrs",mins,"mins and",secs,"seconds.\n"))
  
  ##################################
  ##################################
  ##################################
  ### --- Results processing --- ###
  ##################################
  ##################################
  ##################################
  
  cat("Reading BGLiMS posterior output...\n")  
  t1 <- proc.time()["elapsed"]
  
  ### --- Read BGLiMS arguments
  bglims.arguments <- as.list(read.table(results.file, header=T, sep=" ", nrows=1))
  bglims.arguments$Likelihood <- as.character(bglims.arguments$Likelihood)
  bglims.arguments$ModelSpacePriorFamily <- as.character(bglims.arguments$ModelSpacePriorFamily)
  n.lines.until.rjmcmc.output <- 3 # There are always three lines of meta-data
  enumerated.posterior.inference <- list("No enumeration was done.")
  if (enumerate.up.to.dim>0) {
    ### -- Read enumerated posterior scores
    enumerated.model.likelihood.scores <- NULL
    for (max.model.dimension in c(0:enumerate.up.to.dim)) { # From 0 to account for the null model
      n.models.this.dimension <- choose(bglims.arguments$V - bglims.arguments$startRJ, max.model.dimension)
      model.scores.this.dimension <- read.table(
        results.file,
        skip = n.lines.until.rjmcmc.output,
        header=FALSE,
        nrows=n.models.this.dimension)
      n.lines.until.rjmcmc.output <- n.lines.until.rjmcmc.output+n.models.this.dimension
      enumerated.model.likelihood.scores <- rbind(enumerated.model.likelihood.scores, model.scores.this.dimension)    
    }
    colnames(enumerated.model.likelihood.scores) <- c("Model", "PosteriorScore")
    enumerated.model.likelihood.scores$Model <- as.character(enumerated.model.likelihood.scores$Model)
    enumerated.posterior.inference <- .ApproxPostProbsByModelEnumeration(enumerated.model.likelihood.scores, model.space.priors, enumerate.up.to.dim)
  }
  
  ### --- Read MCMC output
  n.rows.written <- bglims.arguments$iterations/bglims.arguments$thin
  mcmc.output <- read.table(
    results.file,
    skip = n.lines.until.rjmcmc.output,
    header=TRUE,
    nrows=n.rows.written)
  Lhalf <- round(nrow(mcmc.output)/2)     	 # Burnin is a half	
  mcmc.output <- mcmc.output[(Lhalf+1):nrow(mcmc.output),]   # Remove burnin

  ### --- Summary table
  posterior.summary.table <- matrix(NA,ncol(mcmc.output)+length(model.space.priors),8)
  colnames(posterior.summary.table) = c("PostProb","Median","CrI_Lower","CrI_Upper",
    "Median_Present","CrI_Lower_Present","CrI_Upper_Present","BF")
  rownames(posterior.summary.table) <- c(colnames(mcmc.output),paste("ModelSizePartition",c(1:length(model.space.priors)),sep=""))
  # Prior probabilties - used for Bayes Factors below
  prior.probs <- rep(NA, nrow(posterior.summary.table))
  names(prior.probs) <- rownames(posterior.summary.table)
  for (c in 1:length(model.space.priors)) {
    # Calculate partition specific covariate specific prior probabilities of inclusion
    if ("Rate" %in% names(model.space.priors[[c]]) ) {
      prior.probs[model.space.priors[[c]]$Variables] <- .ModelSpaceSpecProb(length(model.space.priors[[c]]$Variables), model.space.priors[[c]]$Rate)
    } else {
      prior.probs[model.space.priors[[c]]$Variables] <- model.space.priors[[c]]$a/(model.space.priors[[c]]$a + model.space.priors[[c]]$b)
    }
    # Calculate posterior on model dimension in each partition
    model.dim.posterior.c <- apply(as.matrix(mcmc.output[,model.space.priors[[c]]$Variables]),MAR=1,function(x)sum(x!=0)) # Wrap in as.matrix incase model space partition has one variable
    posterior.summary.table[paste("ModelSizePartition",c,sep=""),c("CrI_Lower", "Median", "CrI_Upper")] <- quantile(model.dim.posterior.c,c(0.025, 0.5, 0.975))
  }
  # Fill in the summary table
  for (v in colnames(mcmc.output)) {
    posterior.summary.table[v,c("CrI_Lower", "Median", "CrI_Upper")] <- quantile(mcmc.output[,v],c(0.025, 0.5, 0.975))
    if (v %in% unlist(lapply(model.space.priors, function(x) x$Variables))) {
      posterior.summary.table[v,c("CrI_Lower_Present", "Median_Present", "CrI_Upper_Present")] <- quantile(mcmc.output[,v][mcmc.output[,v]!=0],c(0.025, 0.5, 0.975) )
      posterior.summary.table[v,"PostProb"] <- length( mcmc.output[,v][mcmc.output[,v]!=0] ) / nrow(mcmc.output)      
      if (enumerate.up.to.dim>0) {
        # Replace with enumeration probs
        posterior.summary.table[v,"PostProb"] <- enumerated.posterior.inference$marg.probs[v]
      }
      posterior.summary.table[v,"BF"] <- .BayesFactor(prior.probs[v], posterior.summary.table[v,"PostProb"])
      if (likelihood %in% c("Weibull", "Cox", "CaseCohort_Prentice", "CaseCohort_Barlow", "Logistic", "CLogLog", "RocAUC_Anchoring") ) {
        # Exponentiate quantiles
        posterior.summary.table[v,c("CrI_Lower", "Median", "CrI_Upper","CrI_Lower_Present", "Median_Present","CrI_Upper_Present")] <- exp(posterior.summary.table[v,c("CrI_Lower", "Median", "CrI_Upper","CrI_Lower_Present", "Median_Present","CrI_Upper_Present")])
      }
    }
  }
  
  t2 <- proc.time()["elapsed"]  
  results.processing.time <- t2-t1
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  cat(paste("Posterior output processed in",hrs,"hrs",mins,"mins and",secs,"seconds.\n"))
  
  ########################
  ########################
  ### --- Clean up --- ###
  ########################
  ########################
  
  if (clean.up.bglims.files) {
    unlink(c(arguments.file, data.file, results.file), recursive = FALSE, force = TRUE)
  }
  hrs <-floor( (t2-t1)/(60*60) )
  mins <- floor( (t2-t1-60*60*hrs)/60 )
  secs <- round(t2-t1-hrs*60*60 - mins*60)
  
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

  # Create dummies for any missing slots in order to create S4 object
  if (is.null(confounders)) {confounders <- c("None")}
  if (is.null(beta.prior.partitions)) {beta.prior.partitions <- list("NA")}
  if (is.null(n)) {n <- nrow(data)}
  results <- methods::new(
    "R2BGLiMS_Results",
    likelihood = likelihood,
    n = n,
    posterior.summary.table = posterior.summary.table,
    enumerate.up.to.dim = enumerate.up.to.dim,
    enumerated.posterior.inference = enumerated.posterior.inference,
    n.iterations = n.iter,
    thin = bglims.arguments$thin,
    model.space.priors = model.space.priors,
    beta.prior.partitions = beta.prior.partitions,
    confounders = confounders,
    run.times=run.times,
    n.covariate.blocks.for.jam = 1,
    bglims.arguments=bglims.arguments,
    mcmc.output=mcmc.output
    )
  
  ########################
  ### --- Finished --- ###
  ########################
  
  cat(paste("Finished.\n"))
  return(results)
}
