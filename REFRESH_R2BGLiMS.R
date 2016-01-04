package.name <- 'R2BGLiMS'
package.location <- '/Users/pauln/Dropbox/Work Projects/R Packages'
github.username <-'pjnewcombe'
commit.message <- NULL
#commit.message <- "Unnessary specification of N for the marginal model removed."

### --- Java code and arguments file
setwd(paste(package.location,"/",package.name,sep=""))
system("rm DESCRIPTION")
system("rm -rf inst/BGLiMS/*")
system( paste("cp -rf '/Users/pauln/Dropbox/Work Projects/RJ Java Package/Java Code/BGLiMS/dist/'* inst/BGLiMS/", sep="") )
system( paste("cp '/Users/pauln/Dropbox/Work Projects/RJ Java Package/Test/Arguments/Arguments_Package.txt' inst/BGLiMS/", sep="") )  

### --- Load Pmisc and refresh the package
library(Pmisc)
RefreshPackage(package.name=package.name,commit.message=commit.message)

do.not.run <- FALSE
if (do.not.run) {
  library(devtools)
  install_github("pjnewcombe/R2BGLiMS")
}
