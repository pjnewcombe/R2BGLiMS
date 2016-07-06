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
  subcohort.var=NULL,
  confounders=NULL,
  predictors=NULL,
  model.space.priors,
  beta.priors=NULL,
  beta.prior.partitions=NULL,
  dirichlet.alphas.for.roc.model=NULL,
  g.prior=FALSE,
  model.tau=FALSE,
  tau=NULL,
  enumerate.up.to.dim=0,
  xTx=NULL,
  z=NULL,
  subcohort.sampling.fraction=NULL,
  casecohort.pseudo.weight=NULL,
  max.fpr=1,
  min.tpr=0,
  initial.model=NULL
) {
	### Pre-processing
  if (likelihood%in%c("Cox", "CaseCohort_Prentice", "CaseCohort_Barlow")) {
    # Re-order rows of data in ascending order of follow-up time
    # Must be done BEFORE extracting individual variables
    data <- data[order(data[,times.var], decreasing=T),]
  }  
  if (!likelihood%in%c("JAM_MCMC", "JAM")) {
    n.start <- nrow(data)
    if (!is.null(predictors)) {
      data <- data[, c(outcome.var, times.var, subcohort.var, predictors)]
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
    # Extract sub-cohort var from the main data
    if (!is.null(subcohort.var) ) {
      subcohort.indicators <- data[,subcohort.var] # NB This is done after the ordering step above
      data <- data[, colnames(data)[!colnames(data)%in%subcohort.var] ]
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
	write(paste("totalNumberOfCovariates",format(V,sci=F)), file = data.file , ncolumns = 1, append = T)
	
  # Varnames: Variable names
	write(var.names, file = data.file , ncolumns = V, append = T)
	
  # VnoRJ: Total number of variables (from the beginning) to be fixed in the model
	write(paste("numberOfCovariatesToFixInModel",format(length(confounders),sci=F)), file = data.file , ncolumns = 1, append = T)
	
	# N: Number of individuals
	write(paste("numberOfIndividuals",format(N,sci=F)), file = data.file , ncolumns = 1, append = T)
	
	### --- Dirichlet alphas (if ROC model)
	if (likelihood %in% c("RocAUC")) {
	  if (length(dirichlet.alphas.for.roc.model)==1) {
	    dirichlet.alphas.for.roc.model <- c(rep(dirichlet.alphas.for.roc.model,V))
	  }
	  for (v in 1:length(dirichlet.alphas.for.roc.model)) {
	    write( dirichlet.alphas.for.roc.model[v], file = data.file , ncolumns = 2, append = T)
	  }
	}
	
	### --- Fixed beta priors
	if (!likelihood %in% c("RocAUC")) {
	  if (is.null(beta.priors)) {
	    write(paste("numberOfCovariatesWithInformativePriors",0), file = data.file , ncolumns = 1, append = T)    
	  } else {
	    write(paste("numberOfCovariatesWithInformativePriors",nrow(beta.priors)), file = data.file , ncolumns = 1, append = T)
	    write("FixedPriors", file = data.file , ncolumns = 1, append = T)
	    for (v in 1:nrow(beta.priors)) {
	      write( c(beta.priors[v,1],beta.priors[v,2]), file = data.file , ncolumns = 2, append = T)
	    }
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
	}
	
	############################
  # --- Covariate matrix --- #
	############################
	
	cat("Writing data into an input file for BGLiMS...\n")
  # Covariate data - different if marginal setup
  if (likelihood %in% c("JAM", "JAM_MCMC")) { # Write summary data
    write(paste("nBlocks",format(length(xTx),sci=F)), file = data.file , ncolumns = 1, append = T)
    write("blockIndices", file = data.file , ncolumns = 1, append = T)
    write(block.indices, file = data.file , ncolumns = length(block.indices), append = T)
    if (likelihood == "JAM_MCMC") {
      for (b in 1:length(xTx)) {
        write.table(xTx[[b]], row.names=F, col.names=F, file = data.file, append = T)
      }
    } else if (likelihood == "JAM") {
      Lt_Inv <- list()
      for (b in 1:length(xTx)) {
        L <- chol(xTx[[b]]) # NB: UPPER triangle. So L'L = X'X (LIKE IN PAPER)
        write.table(L, row.names=F, col.names=F, file = data.file, append = T) # Multiply by L' (like in JAVA_test)
        cat("Taking Cholesky decomposition of block",b,"...\n")
        Lt_Inv[[b]] <- solve(t(L)) # Take TRANSPOSE inverse for below
      }      
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
  if (likelihood %in% c("Logistic", "Weibull", "Cox", "CaseCohort_Prentice", "CaseCohort_Barlow", "RocAUC", "RocAUC_Anchoring")) {
    write(t(as.integer(outcome)), file = data.file , ncolumns = N, append = T)    
  } else if (likelihood %in% c("Gaussian","GaussianConj")) {
    if (likelihood == "GaussianConj") { outcome <- outcome - mean(outcome) }
    write(t(outcome), file = data.file , ncolumns = N, append = T)        
  } else if (likelihood %in% c("JAM_MCMC")) {
    write(t(z), file = data.file , ncolumns = V, append = T)        
  } else if (likelihood %in% c("JAM")) {
    # Multiply z by inverse cholesky decomposition
    z_L <- z
    for (b in 1:length(xTx)) {
      block.vars <- c(block.indices[b]:(block.indices[b+1]-1))
      z_L[block.vars] <- Lt_Inv[[b]] %*% z[block.vars]
    }
    write(t(z_L), file = data.file , ncolumns = V, append = T)        
  }
  if (likelihood=="Weibull") {
    write(t(times), file = data.file , ncolumns = N, append = T)    
  }
  if (likelihood %in% c("CaseCohort_Prentice", "CaseCohort_Barlow") ) {
    write(t(as.integer(subcohort.indicators)), file = data.file , ncolumns = N, append = T)    
    write(casecohort.pseudo.weight, file = data.file , ncolumns = 1, append = T)
  }
  if (likelihood == "CaseCohort_Barlow") {
    write(subcohort.sampling.fraction, file = data.file , ncolumns = 1, append = T)
  }
  if (likelihood %in% c("RocAUC","RocAUC_Anchoring") ) {
    write(max.fpr, file = data.file , ncolumns = N, append = T)    
    write(min.tpr, file = data.file , ncolumns = N, append = T)    
  }
  
  ########################
  ### --- Finished --- ###
  ########################  
}
