# For each gene calculate correlations between its mRNA levels and the 
# methylation levels of nearby loci. This script should be batch processed and
# 'tissue' should be supplied as a command line argument equal to CRBLM, PONS,
# FCTX or TCTX.
###########################################################################

library(Biobase)
source("functions/cor_cross_eset.r")
options(stringsAsFactors = FALSE)

# Read in command line arguments ------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Evaluate the expressions containing each argument
if(length(args) == 0){
  stop("tissue argument not supplied.")
} else {
  for(i in seq_along(args)){
    eval(parse(text = args[[i]]))
  }
}

tissue <- match.arg(tissue, c("PONS", "CRBLM", "FCTX", "TCTX"))


# Load adjusted data ------------------------------------------------------
load("data/crossdata-correlations-meth-exp/eset-exp-adjusted.rda")
load("data/crossdata-correlations-meth-exp/eset-meth-adjusted.rda")


# Subset ExpressionSets ---------------------------------------------------

# Identify individuals present in both datasets
common.ids <- intersect(eset.exp[, eset.exp$tissue == tissue]$individual,
  eset.meth[, eset.meth$tissue == tissue]$individual)

cat(length(common.ids), "individual IDs present in both datasets.\n")

# Region-specific eSet objects
exp.sub <- eset.exp[, eset.exp$individual %in% common.ids & 
                      eset.exp$tissue == tissue]

meth.sub <- eset.meth[, eset.meth$individual %in% common.ids &
                      eset.meth$tissue == tissue]

# Swap geo ids for individual ids
sampleNames(exp.sub) <- exp.sub$individual
sampleNames(meth.sub) <- meth.sub$individual


# Construct cross-platform probe look-up table ----------------------------
if( exists("probe.key") ) {
  
  if(!all(names(probe.key) %in% c("target", "meth", "exp"))) {
    stop("Supplied probe.key data.frame must contain three columns named:",
         "target, meth and exp", call. = FALSE)
  }
  
  cat("Supplied probe.key contains", length(unique(probe.key$meth)), 
      "methylation probes and", length(unique(probe.key$exp)), 
      "mRNA probes with", length(unique(probe.key$target)), "common targets\n")
} else {
  # Identify meth and exp probes with common targets 
  meth.key <- with(fData(meth.sub), data.frame(target = symbol, meth = id))
  exp.key <- with(fData(exp.sub), data.frame(target = symbol, exp = id))
  probe.key <- merge(meth.key, exp.key, by = "target", sort = FALSE)
  
  cat("Identified", length(unique(probe.key$meth)), "methylation probes and",
      length(unique(probe.key$exp)), "mRNA probes with",
      length(unique(probe.key$target)), "common targets\n")
}

# Correlate expression/methylation residuals ------------------------------
cors.df <- cor_cross_eset(meth.sub, exp.sub, probe.key)

write.csv(cors.df, file = file.path("results",
  paste("crossdata-meth-exp-correlations-", tolower(tissue), ".csv", sep = "")), row.names = FALSE)
