#' Writes a data file in the format required for Java MCMC program
#' @export
#' @title Write Java MCMC format data file
#' @name .WriteData
#' @inheritParams R2BGLiMS
#' @param data.file Desired path for .txt data file to be written to
#' @return NA
#' @author Paul Newcombe
.WriteData <- function(
  data.file,
  likelihood,
  data,
  outcome.var=NULL,
  times.var=NULL,
  confounders=NULL,
  predictors=NULL,
  model.space.priors,
  beta.priors=NULL,
  beta.prior.partitions=NULL,
  g.prior=FALSE,
  model.tau=FALSE,
  tau=NULL,
  xtx.ridge.term=0,
  enumerate.up.to.dim=0,
  n=n,
  xTx=NULL, # X.ref * X.ref
  z=NULL, # X'y
  ns.each.ethnicity=NULL,
  initial.model=NULL,
  trait.variance = NULL,
  mrloss.w = 0,
  mrloss.function = "variance",
  mrloss.marginal.causal.effects = NULL,
  mrloss.marginal.causal.effect.ses = NULL
) {
	### Pre-processing
  if (likelihood%in%c("JAM_MCMC", "JAM")) {
    N <- 1
    n.ethnicities <- 1
    if (!is.null(ns.each.ethnicity)) {
      V <- length(z[[1]])
      var.names <- colnames(xTx[[1]])
      n.blocks <- 1
      n.ethnicities <- length(ns.each.ethnicity)
      block.sizes <- V
      block.indices <- c(1,(V+1))
    } else {
      V <- length(z)
      var.names <- unlist(lapply(xTx,colnames))
      n.blocks <- length(xTx)
      block.sizes <- unlist(lapply(xTx,ncol))
      block.indices <- rep(1,length(xTx)+1)
      for (b in 1:length(xTx)) {
        block.indices[b+1] <- block.indices[b] + block.sizes[b]
      }
    }
  } else {
    n.start <- nrow(data)
    if (!is.null(predictors)) {
      data <- data[, c(outcome.var, times.var, predictors)]
    }
    for (v in colnames(data)) {
      data <- data[!is.na(data[,v]), ]
    }
    # Extract outcome var
    outcome <- data[,outcome.var]
    data <- data[, colnames(data)[!colnames(data)%in%c(outcome.var)] ]
    # Extract survival times var from the main data
    if (!is.null(times.var) ) {
      times <- data[,times.var]
      data <- data[, colnames(data)[!colnames(data)%in%times.var] ]
    }

    V <- ncol(data)
    N <- nrow(data)
    var.names <- colnames(data)
    cat(paste((n.start-N),"observations deleted due to missingness\n"))
    if (N==0) stop("All observations have been deleted due to missingness!\n")  
  }
  
	### Writing
  
	# Model
  if (likelihood == "JAM" & !is.null(trait.variance)) {
    write("JAMv2", file = data.file , ncolumns = 1)
  } else {
    write(likelihood, file = data.file , ncolumns = 1)
  }
  
	# V: Total number of variables
	write(paste("totalNumberOfCovariates",format(V,sci=F)), file = data.file , ncolumns = 1, append = T)
	
  # Varnames: Variable names
	write(var.names, file = data.file , ncolumns = V, append = T)
	
  # VnoRJ: Total number of variables (from the beginning) to be fixed in the model
	write(paste("numberOfCovariatesToFixInModel",format(length(confounders),sci=F)), file = data.file , ncolumns = 1, append = T)
	
	# N: Number of individuals
	write(paste("numberOfIndividuals",format(N,sci=F)), file = data.file , ncolumns = 1, append = T)
	
	### --- Fixed beta priors
	if (is.null(beta.priors)) {
	  write(paste("numberOfCovariatesWithInformativePriors",0), file = data.file , ncolumns = 1, append = T)    
	} else {
	  write(paste("numberOfCovariatesWithInformativePriors",nrow(beta.priors)), file = data.file , ncolumns = 1, append = T)
	  write("FixedPriors", file = data.file , ncolumns = 1, append = T)
	  for (v in 1:nrow(beta.priors)) {
	    write( c(beta.priors[v,1],beta.priors[v,2]), file = data.file , ncolumns = 2, append = T)
	  }
	}    
	
	### --- Hierarchical beta priors
	write(paste("numberOfHierarchicalCovariatePriorPartitions",format(length(beta.prior.partitions),sci=F) ), file = data.file , ncolumns = 1, append = T)
	if (length(beta.prior.partitions)>0) {
	  # --- Partion indices
	  partitioned.covariates <- var.names[var.names %in% unlist(lapply(beta.prior.partitions, function(x) x$Variables))] # Make sure order is consistent with the above
	  partition.indices <- rep(NA, length(partitioned.covariates))
	  for (c in 1:length(beta.prior.partitions)) {
	    partition.indices[which(partitioned.covariates %in% beta.prior.partitions[[c]]$Variables)] <- c
	  }
	  write("hierarchicalCovariatePriorPartitionPicker", file = data.file , ncolumns = 1, append = T)
	  write(t(partition.indices), file = data.file , ncolumns = length(partition.indices), append = T)
	  # --- Partition-specific hyper priors
	  write("BetaPriorPartitionVarianceHyperPriors", file = data.file , ncolumns = 1, append = T)
	  for (c in 1:length(beta.prior.partitions)) {
	    if (beta.prior.partitions[[c]]$Family=="Uniform") {
	      hyperparam1 <- beta.prior.partitions[[c]]$UniformA
	      hyperparam2 <- beta.prior.partitions[[c]]$UniformB
	    } else if (beta.prior.partitions[[c]]$Family=="Gamma") {
	      hyperparam1 <- beta.prior.partitions[[c]]$GammaA
	      hyperparam2 <- beta.prior.partitions[[c]]$GammaB
	    }
	    write(paste(beta.prior.partitions[[c]]$Family,
	                format(hyperparam1,sci=F),
	                format(hyperparam2,sci=F),
	                format(beta.prior.partitions[[c]]$Init,sci=F)),
	          file = data.file, append = T)
	  }
	}
	
	# Initial model
	if (is.null(initial.model)) {
	  write("initialModelOption 0", file = data.file , ncolumns = 1, append = T)    
	} else if (length(initial.model)==1) {
	  write("initialModelOption 1", file = data.file , ncolumns = 1, append = T)    
	} else {
	  write("initialModelOption 2", file = data.file , ncolumns = 1, append = T)    
	  write("InitialModel", file = data.file , ncolumns = 1, append = T)
	  write(t(initial.model), file = data.file , ncolumns = length(initial.model), append = T)    
	}
	
	# Conjugate-only modelling options
	if (likelihood %in% c("GaussianConj", "JAM")) {
	  write(paste("useGPrior", as.integer(g.prior)), file = data.file , ncolumns = 1, append = T)
	  write(paste("tau",format(tau,sci=F)), file = data.file, ncolumns = 1, append = T)      
	  write(paste("modelTau",as.integer(model.tau)), file = data.file , ncolumns = 1, append = T)      
	  write(paste("enumerateUpToDim",format(enumerate.up.to.dim,sci=F)), file = data.file , ncolumns = 1, append = T)
	  if (!is.null(trait.variance)) {
	    YtY <- trait.variance*(n-1)
	    write(paste("YtY",format(YtY,sci=F)), file = data.file , ncolumns = 1, append = T)
	    write(paste("nForJamV2Likelihood",format(n,sci=F)), file = data.file , ncolumns = 1, append = T)
	  }
	}

	############################
  # --- Covariate matrix --- #
	############################
	
	cat("Writing data into an input file for BGLiMS...\n")
  # Covariate data - different if marginal setup
  if (likelihood %in% c("JAM", "JAM_MCMC")) { # Write summary data
    write(paste("nEthnicities",format(n.ethnicities,sci=F)), file = data.file , ncolumns = 1, append = T)
    write(paste("nBlocks",format(n.blocks,sci=F)), file = data.file , ncolumns = 1, append = T)
    write("blockIndices", file = data.file , ncolumns = 1, append = T)
    write(block.indices, file = data.file , ncolumns = length(block.indices), append = T)
    if (likelihood == "JAM_MCMC") {
      for (b in 1:length(xTx)) {
        write.table(xTx[[b]], row.names=F, col.names=F, file = data.file, append = T)
      }
    } else if (likelihood == "JAM") {
      Lt_Inv <- list()
      for (b in 1:length(xTx)) { # Could be blocks or ethnicities
        if (xtx.ridge.term!=0) {
          diag(xTx[[b]]) <- diag(xTx[[b]]) + xtx.ridge.term # Add a constrant to diagonal first
        }
        L <- chol(xTx[[b]]) # NB: UPPER triangle. So L'L = X'X (LIKE IN PAPER)
        write.table(L, row.names=F, col.names=F, file = data.file, append = T) # Multiply by L' (like in JAVA_test)
        cat("Taking Cholesky decomposition of block",b,"...\n")
        Lt_Inv[[b]] <- solve(t(L)) # Take TRANSPOSE inverse for below
      }      
    }
    # MR Pleiotropic loss function stuff
    write(mrloss.w, file = data.file , ncolumns = 1, append = T)
    if (mrloss.w != 0) {
      write(mrloss.function, file = data.file , ncolumns = 1, append = T)
      write(t(mrloss.marginal.causal.effects), file = data.file , ncolumns = length(mrloss.marginal.causal.effects), append = T)    
      write(t(mrloss.marginal.causal.effect.ses), file = data.file , ncolumns = length(mrloss.marginal.causal.effect.ses), append = T)    
    }
  } else { 
    # Write IPD Covariate data
    if (likelihood == "GaussianConj") {
      # Mean centre covariates for Gaussian conjugate model
      data <- apply(data,MAR=2,function(x) x-mean(x))
    }
    write.table(data, row.names=F, col.names=F, file = data.file , append = T)    
  }
  
  # Vector of outcomes
  if (likelihood %in% c("Logistic", "CLogLog", "Weibull")) {
    write(t(as.integer(outcome)), file = data.file , ncolumns = N, append = T)    
  } else if (likelihood %in% c("Gaussian","GaussianConj")) {
    if (likelihood == "GaussianConj") { outcome <- outcome - mean(outcome) }
    write(t(outcome), file = data.file , ncolumns = N, append = T)        
  } else if (likelihood %in% c("JAM_MCMC")) {
    write(t(z), file = data.file , ncolumns = V, append = T)        
  } else if (likelihood %in% c("JAM")) {
    # Multiply z by inverse cholesky decomposition
    if (length(ns.each.ethnicity)>1) {
      z_L <- NULL
      for (e in 1:length(xTx)) {
        z_L <- c(z_L, Lt_Inv[[e]] %*% z[[e]])
      }
    } else {
      z_L <- z
      for (b in 1:length(xTx)) {
        block.vars <- c(block.indices[b]:(block.indices[b+1]-1))
        z_L[block.vars] <- Lt_Inv[[b]] %*% z[block.vars]
      }
    }
    write(t(z_L), file = data.file , ncolumns = V, append = T)        
  }
  if (likelihood=="Weibull") {
    write(t(times), file = data.file , ncolumns = N, append = T)    
  }

  ########################
  ### --- Finished --- ###
  ########################  
}
