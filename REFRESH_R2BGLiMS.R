library(roxygen2)
package.name <- 'R2BGLiMS'
package.location <- '/Users/paul/Documents/Work Projects/R Packages'

#######################################################
# --- Copy in latest Java code and arguments file --- #
#######################################################
setwd(paste(package.location,"/",package.name,sep=""))
system("rm DESCRIPTION")
system("rm -rf inst/BGLiMS/*")
system( paste("cp -rf '/Users/paul/Documents/Work Projects/RJ Java Package/Java Code/BGLiMS/dist/'* inst/BGLiMS/", sep="") )
system( paste("cp '/Users/paul/Documents/Work Projects/RJ Java Package/Test/Arguments/Arguments_Package.txt' inst/BGLiMS/", sep="") )  

####################################
# --- Write `DESCRIPTION' file --- #
####################################
write.dcf( 
  x = list(
    Package = package.name, 
    Type = "Package",
    Title = "R interface to the java package BGLiMS (Bayesian Generalised Linear Model Selection).",
    Version = paste("0.1",format(Sys.time(), "%d-%m-%Y"),sep="-"),
    Date=format(Sys.time(), "%d-%m-%Y"),
    Author = "Paul J Newcombe <paul.newcombe@mrc-bsu.cam.ac.uk>",
    Maintainer = "Paul J Newcombe <paul.newcombe@mrc-bsu.cam.ac.uk>",
#    Imports = "MendelianRandomization, jsonlite", # jsonlite for an error Frank Dudbridge saw.
    Imports = "jsonlite", # jsonlite for an error Frank Dudbridge saw.
    Description = "Provides an R interface to a Java package for fitting 
    Bayesian GLMs. Approximate posterior samples are drawn using an MCMC 
    sampler with a (Reversible Jump) Metropolis-Hastings acceptance ratio. 
    Currently, Logistic and Weibull, and Linear regression models are 
    available, as well as JAM for genetic summary data fine-mapping, 
    JAMPred for polygenic risk scores, and JAMMR for summary data Mendelian
    Randomization.",
    License = "GPL (>=2)"),
  file = file.path(package.location,package.name,"DESCRIPTION"),
  append = FALSE
)

#####################################################
# --- Run Roxygen to generate NAMESPACE and Rds --- #
#####################################################
roxygenise(package.dir=file.path(package.location,package.name),clean=TRUE)

######################
# --- Re-install --- #
######################
system(command=paste("R CMD INSTALL '",file.path(package.location,package.name),"'",sep=""))
