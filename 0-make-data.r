# Download data from GEO and save cleaned expressionSet objects as rda files in
# the data directory

# Check for necessary packages --------------------------------------------
packages <- c("Biobase", "plyr", "devtools")
package.paths <- find.package(packages, quiet = TRUE)

if(length(package.paths) < length(packages)) {
  status <- sapply(packages, function(x) sum(grepl(x, package.paths)))
  cat("\nThe following packages are missing:\n", names(status)[status == 0])
  
  stop("\n\nInstall the packages listed above before proceeding", call. = FALSE)
}


# Download data -----------------------------------------------------------
source("1-download-data.r")


# Clean data and export rda files -----------------------------------------
source("2-reformat-data.r")


# Check integrity of rda files --------------------------------------------
exp.verified <- devtools:::md5("data/eset-exp.rda") == "3f9afdefd2469691f9d54f916d33089a"
meth.verified <- devtools:::md5("data/eset-meth.rda") == "800538004a90b6e3b9c11b3497fc7db9"
micro.verified <- devtools:::md5("data/eset-micro.rda") == "604b2f6dcb5de62cd80e4bdd8c5fc658"


# Remove raw geo data -----------------------------------------------------
if(all(exp.verified, meth.verified, micro.verified)) {
  file.remove(c("data/geo-raw-data.rda", dir("data/geo", full.names = TRUE)))  
} else {
  stop("rda files generated from GEO downloads differ from Aaron's.", call. = FALSE)
}


