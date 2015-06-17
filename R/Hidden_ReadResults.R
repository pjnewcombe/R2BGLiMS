#' Reads in a .txt results file, output by my Java program, and returns as a list
#' @export
#' @title Read Java Results
#' @name .ReadResults
#' @param results.file results file to be read (default is NULL)
#' @param mi.results.files character vector of results files to be read in if analysing multiply imputed datasets (default is NULL)
#' @param burnin.denominator fraction of results (from the end) to keep
#' @param thin.factor keep every nth saved iteration
#' @param first.n.mil.its how many million iterations to read in from the beginning. Useful if you would like to explore whether
#' convergence could have been achieved with less iterations (default is NULL)
#' 
#' @return A Reversible Jump results object is returned. This is a list of two elements: "args" which records various modelling
#' arguments used in the analysis, and "results" - a matrix containing all saved posterior samples from the analysis. Columns
#' correspond to parameters and each row contains values from a particular itertation, with 0 indicating exclusion from the model.
#' 
#' @author Paul Newcombe
.ReadResults <- function(
  results.file=NULL,
  mi.results.files=NULL,
  burnin.denominator=2,
  thin.factor=1,
  first.n.mil.its=NULL) {
  cat("\n\nReading in the BGLiMS results file...\n")  
  if (is.null(results.file)) {
    results.file <- mi.results.files[1]
  }
  resultsRead <- list()	# Initialise list
  
  # Read arguments
  resultsRead$args <- as.list(read.table(results.file, header=T, sep=" ", nrows=1))
  resultsRead$args$Likelihood <- as.character(resultsRead$args$Likelihood)
  resultsRead$args$ModelSpacePriorFamily <- as.character(resultsRead$args$ModelSpacePriorFamily)
  if (!is.null(first.n.mil.its)) {
    resultsRead$args$iterations <- first.n.mil.its*1000000 	# Number of iterations
  }
  n.rows.written <- resultsRead$args$iterations/resultsRead$args$thin
  
  # Get model space prior information FROM THE RESULTS FILE
  modelSpacePriorQuantities <- scan(results.file, skip = 2, nlines = 1)
  if (resultsRead$args$ModelSpacePriorFamily=="Poisson") {
    resultsRead$args$poisson.mu <- c(1:resultsRead$args$nRjComp)    
  } else if (resultsRead$args$ModelSpacePriorFamily=="BetaBinomial") {
    resultsRead$args$beta.binom.a <- c(1:resultsRead$args$nRjComp)
    resultsRead$args$beta.binom.b <- c(1:resultsRead$args$nRjComp)    
  }
  for (c in 1:resultsRead$args$nRjComp) {
    if (resultsRead$args$ModelSpacePriorFamily=="Poisson") {
      resultsRead$args$poisson.mu[c] <- modelSpacePriorQuantities[c]
    } else if (resultsRead$args$ModelSpacePriorFamily=="BetaBinomial") {
      resultsRead$args$beta.binom.a[c] <- modelSpacePriorQuantities[2*(c-1)+1]
      resultsRead$args$beta.binom.b[c] <- modelSpacePriorQuantities[2*(c-1)+2]
    }
  }
  if (resultsRead$args$nRjComp > 1) {
    resultsRead$args$modSpaceSplits <- c(resultsRead$args$nRjComp+1)
    for (c in 1:(resultsRead$args$nRjComp+1) ) {
      if (resultsRead$args$ModelSpacePriorFamily=="Poisson") {
        resultsRead$args$modSpaceSplits[c] <- modelSpacePriorQuantities[resultsRead$args$nRjComp+c]
      } else if (resultsRead$args$ModelSpacePriorFamily=="BetaBinomial") {
        resultsRead$args$modSpaceSplits[c] <- modelSpacePriorQuantities[resultsRead$args$nRjComp*2+c]        
      }
    }		
  }
  
  ### Read in posterior scores for single SNP models
  if (resultsRead$args$allModelScoresUpToDim>0) {
    # Null and single dimension models
    n.model.scores.dim1 <- (resultsRead$args$V - resultsRead$args$startRJ) + 1
    skip.to.main.results=3+n.model.scores.dim1
    resultsRead$model.scores <- read.table(
      results.file,
      skip = 3,
      header=FALSE,
      nrows=n.model.scores.dim1)
    # Two dimension models
    if (resultsRead$args$allModelScoresUpToDim>=2) {
      n.model.scores.dim2 <- choose(resultsRead$args$V - resultsRead$args$startRJ,2)
      model.scores.dim2 <- read.table(
        results.file,
        skip = 3+n.model.scores.dim1,
        header=FALSE,
        nrows=n.model.scores.dim2)
      resultsRead$model.scores <- rbind(resultsRead$model.scores, model.scores.dim2)
      skip.to.main.results <- skip.to.main.results + n.model.scores.dim2
    }
    if (resultsRead$args$allModelScoresUpToDim>=3) {
      n.model.scores.dim3 <- choose(resultsRead$args$V - resultsRead$args$startRJ,3)
      model.scores.dim3 <- read.table(
        results.file,
        skip = 3+n.model.scores.dim1+n.model.scores.dim2,
        header=FALSE,
        nrows=n.model.scores.dim3)
      resultsRead$model.scores <- rbind(resultsRead$model.scores, model.scores.dim3)
      skip.to.main.results <- skip.to.main.results + n.model.scores.dim3
    }
    if (resultsRead$args$allModelScoresUpToDim>=4) {
      n.model.scores.dim4 <- choose(resultsRead$args$V - resultsRead$args$startRJ,4)
      model.scores.dim4 <- read.table(
        results.file,
        skip = 3+n.model.scores.dim1+n.model.scores.dim2+n.model.scores.dim3,
        header=FALSE,
        nrows=n.model.scores.dim4)
      resultsRead$model.scores <- rbind(resultsRead$model.scores, model.scores.dim4)
      skip.to.main.results <- skip.to.main.results + n.model.scores.dim4
    }
    if (resultsRead$args$allModelScoresUpToDim>=5) {
      n.model.scores.dim5 <- choose(resultsRead$args$V - resultsRead$args$startRJ,5)
      model.scores.dim5 <- read.table(
        results.file,
        skip = 3+n.model.scores.dim1+n.model.scores.dim2+n.model.scores.dim3+n.model.scores.dim4,
        header=FALSE,
        nrows=n.model.scores.dim5)
      resultsRead$model.scores <- rbind(resultsRead$model.scores, model.scores.dim5)
      skip.to.main.results <- skip.to.main.results + n.model.scores.dim5
    }
    colnames(resultsRead$model.scores) <- c("Model", "PosteriorScore")
  } else {
    skip.to.main.results=3    
  }
  
  if (is.null(mi.results.files)) {
    resultsRead$results <- read.table(
      results.file,
      skip = skip.to.main.results,
      header=TRUE,
      nrows=n.rows.written)
    # Remove burnin  and thin if asked
    Lhalf <- round(nrow(resultsRead$results)/burnin.denominator)   		 # Burnin is a half	
    resultsRead$results <- resultsRead$results[(Lhalf+1):nrow(resultsRead$results),]
    if(thin.factor > 0) {
      keepIts <- seq(1,nrow(resultsRead$results),by = thin.factor)
      resultsRead$results <- resultsRead$results[keepIts,]
    }
  } else {
    # Reading in and combining results from different multiple imputation chains
    resultsRead$args$n.mi.chains <- length(mi.results.files)
    resultsRead$results <- NULL
    for (results.file.ch in mi.results.files) {
      results.ch <- read.table(
        results.file.ch,
        skip = skip.to.main.results,
        header=TRUE,
        nrows=n.rows.written)
      # Remove burnin  and thin if asked
      Lhalf <- round(nrow(results.ch)/burnin.denominator)     	 # Burnin is a half	
      results.ch <- results.ch[(Lhalf+1):nrow(results.ch),]
      if(thin.factor > 0) {
        keepIts <- seq(1,nrow(results.ch),by = thin.factor)
        results.ch <- results.ch[keepIts,]
      }      
      # Merge
      resultsRead$results <- rbind(resultsRead$results, results.ch)
    }
  }
  
  # Correct number of saved iterations and return list
  cat("... finished reading BGLiMS results.\n")
  return(resultsRead)
}
