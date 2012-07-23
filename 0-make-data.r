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

