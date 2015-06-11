#' Writes a data file in the format required for Java MCMC program
#' @export
#' @title Write Java MCMC format data file
#' @name .WriteData
#' @param data.file Desired path for .txt data file to be written to
#' @param likelihood Type of model to fit. Current options are "Logistic" (for binary data), "Weibull" (for survival data), 
#' "Gaussian" (for continuous data), "GaussianMarg" (for analysis of univariate associations from Gaussian linear 
#' regressions) and "GaussianMargConj" (for analysis under a marginal conjugate linear regression model).
#' @param data Matrix of data to write - rows indiviuals, columns variables.
#' @param predictors vector of predictors. Leave as the default NULL to include all variables available in the data.frame will be used.
#' @param outcome.var Which column in data contains the binary outcome for logistic and survival data, or the integer count for
#' Poisson data (default "Disease")
#' @param times.var If survival data or Poisson data, the column in data which contains the follow-up times (default NULL)
#' @param xTx GaussianMarg and GaussianMargConj ONLY: List containing each block's plug-in estimate for X'X.
#' @param z GaussianMarg and GaussianMargConj ONLY: Vector of quantities calculated from the summary statistics.
#' @param sigma2_invGamma_a GaussianMarg and GaussianMargConj ONLY: Inverse-Gamma parameter one for the residual
#' precision. For the conjugate model this parameter is intergrated out, and this may be provided as in the
#' Bottolo and Richardson 2010 notation. For an informative prior for the non-conjugate model
#' (taking into account Java parameterisation) choose N/2.
#' @param sigma2_invGamma_b GaussianMarg and GaussianMargConj ONLY: Inverse-Gamma parameter one for the residual
#' precision. For the conjugate model this parameter is intergrated out, and this may be provided as in the
#' Bottolo and Richardson 2010 notation. For an informative prior for the non-conjugate model
#' (taking into account Java parameterisation) choose N/(2*variance estimate).
#' @param g.prior GaussianMargConj ONLY: Whether to use a g-prior for the beta's - i.e. a multivariate normal 
#' with correlation structure proportional to sigma^2*X'X^-1 or to use independence priors (default = FALSE).
#' @param tau GaussianMargConj ONLY: Value to use for sparsity parameter tau (tau*sigma^2 parameterisation).
#' Default residual.var.n. If modelling this parameter, this value is used to center the Zellner-Siow prior
#' and as an intial value.
#' @param model.tau GaussianMargConj ONLY: Whether to model tau or not (default FALSE). If set to true,
#' then a Zellner-Siow prior is used, centred on the value provide by tau. The value provided in tau is also used as the
#' initial value.
#' @param tau.proposal.sd GaussianMargConj ONLY: When modelling tau an initial SD to use in the adaption
#' of the proposal distribution. Under the conjugate model, tau makes use of the BGLiMS parameter `betaPriorSd'. Thus
#' the parameters is used on a very different scale to normal and as such it is best to specify this seperately.
#' Defaults to 0.05.
#' @param all.model.scores.up.to.dim GaussianMargConj ONLY: When NOT modelling tau, whether to output the posterior scores
#' for every possible model of dimension up to the integer specified. Currenly maximum allowed dimension is 2. Setting
#' to default 0 means this is not carried out. Can be used as an alternative to running the RJMCMC.
#' 
#' @param block.indices If Guassian marginal tests are being analysed, the external xTx data may be divided
#' into blocks, to simplify inversion. This vector should contain the indices of the block break points (default NULL)
#' @param cluster.var If hierarchical data and random intercepts are required, the column in data contains the clustering variable (default NULL)
#' @param initial.model Optionally, an initial model can be provided as a vector of 0's and 1's. Default is NULL
#' and the null model is used. If set to 1, the saturated model is used.
#' @return NA
#' @author Paul Newcombe
.WriteData <- function(
  data.file,
  likelihood,
  data,
  predictors=NULL,
  confounders=NULL,
  outcome.var=NULL,
  times.var=NULL,
  xTx=NULL,
  z=NULL,
  sigma2_invGamma_a=NULL,
  sigma2_invGamma_b=NULL,
  g.prior=FALSE,
  tau=NULL,
  model.tau=FALSE,
  tau.proposal.sd=0.05,
  all.model.scores.up.to.dim=0,
  residual.var=NULL,
  residual.var.n=NULL,  
  cluster.var=NULL,
  beta.priors=NULL,
  model.space.priors,
  initial.model=NULL
) {
	### Pre-processing
  if (!likelihood%in%c("GaussianMarg", "GaussianMargConj")) {
    n.start <- nrow(data)
    if (!is.null(predictors)) {
      data <- data[, c(outcome.var, times.var, cluster.var, predictors)]
    }
    for (v in colnames(data)) {
      data <- data[!is.na(data[,v]), ]
    }
    # Extract disease var
    disease <- data[,outcome.var]
    data <- data[, colnames(data)[!colnames(data)%in%c(outcome.var)] ]
    # Extract survival times var
    if (!is.null(times.var) ) {
      times <- data[,times.var]
      data <- data[, colnames(data)[!colnames(data)%in%times.var] ]
    }
    # Extract cluster var
    if (!is.null(cluster.var) ) {
      clusters <- as.integer(as.factor(as.integer(data[,cluster.var])))	# This makes range 1,..nClusters
      data <- data[, colnames(data)[!colnames(data)%in%cluster.var] ]
      n.clusters <- length(unique(clusters))
    }
    
    V <- ncol(data)
    N <- nrow(data)
    var.names <- colnames(data)
    cat(paste((n.start-N),"observations deleted due to missingness"))
  } else {
    V <- length(z)
    N <- 1
    var.names <- unlist(lapply(xTx,colnames))
    block.sizes <- unlist(lapply(xTx,ncol))
    block.indices <- rep(1,length(xTx)+1)
    for (b in 1:length(xTx)) {
      block.indices[b+1] <- block.indices[b] + block.sizes[b]
    }
  }
  
	### Writing
	# Model
  write(likelihood, file = data.file , ncolumns = 1)    
	# V: Total number of variables
	write(V, file = data.file , ncolumns = 1, append = T)
  # Varnames: Variable names
	write(var.names, file = data.file , ncolumns = V, append = T)
  # VnoRJ: Total number of variables (from the beginning) to be fixed in the model
	write(length(confounders), file = data.file , ncolumns = 1, append = T) 				# from which variable RJ is performed
  # N: Number of individuals
	write(N, file = data.file , ncolumns = 1, append = T)
  # R: Number of clusters
	if (!is.null(cluster.var) ) {
		write( n.clusters, file = data.file , ncolumns = 1, append = T)				
	} else {
		write(0, file = data.file , ncolumns = 1, append = T)		
	}
  # N x V matrix of covariate values
  cat("\n\nWriting a BGLiMS formatted datafile...\n")
  # All Guassian models
  if (length(grep("Gaussian",likelihood))==1) {
    write(sigma2_invGamma_a, file = data.file , ncolumns = 1, append = T)
    write(sigma2_invGamma_b, file = data.file , ncolumns = 1, append = T)    
  }
  # Conjugate models
  if (length(grep("Conj",likelihood))==1) {
    write(as.integer(g.prior), file = data.file , ncolumns = 1, append = T)
    write(tau, file = data.file, ncolumns = 1, append = T)      
    write(as.integer(model.tau), file = data.file , ncolumns = 1, append = T)      
    write(tau.proposal.sd, file = data.file, ncolumns = 1, append = T)      
    write(all.model.scores.up.to.dim, file = data.file , ncolumns = 1, append = T)
  }
  # Covariate data - different if marginal setup
  if (length(grep("Marg",likelihood))==1) { # Write summary data
    write(length(xTx), file = data.file , ncolumns = 1, append = T)
    write(block.indices, file = data.file , ncolumns = length(block.indices), append = T)
    if (likelihood == "GaussianMarg") {
      for (b in 1:length(xTx)) {
        write.table(xTx[[b]], row.names=F, col.names=F, file = data.file, append = T)
      }
    } else if (likelihood == "GaussianMargConj") {
      L_Inv <- list()      
      for (b in 1:length(xTx)) {
        L <- chol(xTx[[b]]) # NB: UPPER triangle. So L'L = X'X
        write.table(L, row.names=F, col.names=F, file = data.file, append = T)
        cat("Taking Cholesky decomposition of block",b,"...\n")
        L_Inv[[b]] <- solve(L) # Take inverse. Check: id <- t(L[[b]]) %*% xTx[[b]] %*% L[[b]]
        cat("...done")
      }      
    }
  } else { # Write IPD Covariate data
    if (likelihood == "GaussianConj") { data <- apply(data,MAR=2,function(x) x-mean(x)) }
    write.table(data, row.names=F, col.names=F, file = data.file , append = T)    
  }
  cat("... finished writing datafile.\n")
  # If R>0 the vector of cluster labels
	if (!is.null(cluster.var) ) {
		write(t(clusters), file = data.file , ncolumns = n.clusters, append = T)
	}
  # Vector of disease labels
  if (likelihood %in% c("Logistic", "Weibull")) {
    write(t(as.integer(disease)), file = data.file , ncolumns = N, append = T)    
  } else if (likelihood %in% c("Gaussian","GaussianConj")) {
    if (likelihood == "GaussianConj") { disease <- disease - mean(disease) }
    write(t(disease), file = data.file , ncolumns = N, append = T)        
  } else if (likelihood %in% c("GaussianMarg")) {
    write(t(z), file = data.file , ncolumns = V, append = T)        
  } else if (likelihood %in% c("GaussianMargConj")) {
    # Multiply z by inverse cholesky decomposition
    L_Inv_z <- z
    for (b in 1:length(xTx)) {
      block.vars <- c(block.indices[b]:(block.indices[b+1]-1))
      L_Inv_z[block.vars] <- L_Inv[[b]] %*% z[block.vars]
    }
    write(t(L_Inv_z), file = data.file , ncolumns = V, append = T)        
  }
  if (!is.null(times.var)) {
    write(t(times), file = data.file , ncolumns = N, append = T)    
  }
  
  ### --- Fixed beta priors
  use.unknown.priors <- TRUE
  if (is.null(beta.priors)) {
    write(0, file = data.file , ncolumns = 1, append = T)    
  } else {
    write(nrow(beta.priors), file = data.file , ncolumns = 1, append = T)
    for (v in 1:nrow(beta.priors)) {
      write( c(beta.priors[v,1],beta.priors[v,2]), file = data.file , ncolumns = 2, append = T)
    }
    if (nrow(beta.priors)==length(predictors)) {
      use.unknown.priors <- FALSE
    }
  }
  
  ### --- Shared unknown beta priors
  if (!use.unknown.priors) {
    # Fixed priors provided for everything, no need for common unknown prior
    write(0, file = data.file , ncolumns = 1, append = T)    
  } else {
    # Use common unknown priors for model space components
    write(length(model.space.priors), file = data.file , ncolumns = 1, append = T)
    if (length(model.space.priors)>1) {
      for (c in 1:(length(model.space.priors)-1)) {
        write(length(model.space.priors[[c]]$Variables), file = data.file , ncolumns = 1, append = T)      
      }
    }    
  }
  
  # Initial model
  if (is.null(initial.model)) {
    write(0, file = data.file , ncolumns = 1, append = T)    
  } else if (length(initial.model)==1) {
    write(initial.model, file = data.file , ncolumns = 1, append = T)    
  } else {
    write(2, file = data.file , ncolumns = 1, append = T)    
    write(t(initial.model), file = data.file , ncolumns = length(initial.model), append = T)    
  }  
}
