# For each gene calculate correlations between its mRNA levels and the 
# methylation levels of nearby loci. This script should be batch processed and
# 'tissue' should be supplied as a command line argument equal to CRBLM, PONS,
# FCTX or TCTX.
###########################################################################

library(Biobase)

library(foreach)
library(doMC)


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

# Register multicore backend ----------------------------------------------
n.cores <- multicore:::detectCores()
registerDoMC(n.cores)


# Load data ---------------------------------------------------------------
source("4a-load-exp-data.r")
source("4b-load-meth-data.r")


# Subset ExpressionSets ---------------------------------------------------

# Identify individuals present in both datasets
common.ids <- intersect(eset.exp[, eset.exp$tissue == region]$individual,
  eset.meth[, eset.meth$tissue == region]$individual)

# Region-specific eSet objects
exp.sub <- eset.exp[, eset.exp$individual %in% common.ids & 
                      eset.exp$tissue == region]

meth.sub <- eset.meth[, eset.meth$individual %in% common.ids &
                      eset.meth$tissue == region]

# Swap geo ids for individual ids
sampleNames(exp.sub) <- exp.sub$individual
sampleNames(meth.sub) <- meth.sub$individual

# Reorder second dataset to match first
meth.sub <- meth.sub[, sampleNames(exp.sub)]
stopifnot(sampleNames(exp.sub) == sampleNames(meth.sub))


# Extract expression matrices ---------------------------------------------

exp.mat <- exprs(exp.sub)
meth.mat <- exprs(meth.sub)

# Remove probes with NA values
if(!"na.count" %in% names(fData(meth.sub))) {
  na.filter <- apply(exprs(meth.sub), 1, function(x) sum(is.na(x))) == 0
  cat("Removed", sum(!na.filter), "methylation probes with NA values.\n")
  
  meth.mat <- meth.mat[na.filter,]
  meth.sub <- meth.sub[na.filter,]
}

# Identify meth and exp probes with common targets 
genes <- intersect(fData(exp.sub)$symbol, fData(meth.sub)$symbol)


# Correct expression for technical variables ------------------------------

exp.lm <- lm(t(exp.mat) ~ exp.sub$gender + exp.sub$age + exp.sub$tissuebank)
meth.lm <- lm(t(meth.mat) ~ meth.sub$gender + meth.sub$age + meth.sub$tissuebank)

exp.mat <- exp.lm$residuals
meth.mat <- meth.lm$residuals


# Correlate expression/methylation residuals ------------------------------

cors.df <- foreach(g = genes, .combine = "rbind") %dopar% {

  ep <- rownames(fData(exp.sub)[fData(exp.sub)$symbol == g,])
  mp <- rownames(fData(meth.sub)[fData(meth.sub)$symbol == g,])
  
  out <- cor(meth.mat[, mp], exp.mat[, ep])
  
  if(length(mp) == 1 | length(ep) == 1) {
    out <- data.frame(Var1 = mp, Var2 = ep, Freq = as.numeric(out))
  } else {
    out <- as.data.frame.table(out)  
  }
  
  data.frame(symbol = g, out)
}

names(cors.df) <- c("symbol", "meth", "exp", "r")
cors.df$pvalue <- WGCNA::corPvalueStudent(cors.df$r, length(common.ids))

write.csv(cors.df, file = file.path("results",
  paste("meth-exp-correlations-nonnormalized-", tolower(region), ".csv", sep = "")),
  row.names = FALSE)
