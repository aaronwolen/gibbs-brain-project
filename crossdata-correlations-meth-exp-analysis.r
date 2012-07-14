# For each gene calculate correlations between its mRNA levels and the 
# methylation levels of nearby loci. This script should be batch processed and
# 'tissue' should be supplied as a command line argument equal to CRBLM, PONS,
# FCTX or TCTX.
###########################################################################

library(Biobase)

library(foreach)
library(doMC)

source("functions/lm_eset.r")

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

# Reorder second dataset to match first
meth.sub <- meth.sub[, sampleNames(exp.sub)]
stopifnot(sampleNames(exp.sub) == sampleNames(meth.sub))


# Extract expression matrices ---------------------------------------------

exp.mat <- exprs(exp.sub)
meth.mat <- exprs(meth.sub)

# Identify meth and exp probes with common targets 
genes <- intersect(fData(exp.sub)$symbol, fData(meth.sub)$symbol)

# Identify probes corresponding to the common targets
meth.probes <- subset(fData(meth.sub), symbol %in% genes)$id
exp.probes <- subset(fData(exp.sub), symbol %in% genes)$id

cat("Identified", length(meth.probes), "methylation probes and ", length(exp.probes),
    "corresponding to the", length(genes), "present in both datasets.")

# Correct expression for technical variables ------------------------------

# Setting lm_eset's use argument = "pairwise.complete" causes it to loop through
# probes so missing values are be handled on a case by case basis. The result
# is a matrix of residuals.

meth.mat <- lm_eset(~ age + tissuebank, eset = meth.sub, 
  probesets = meth.probes, use = "pairwise.complete", resid.only = TRUE)
meth.mat <- t(meth.mat)
  
exp.mat <- lm_eset(~ age + tissuebank, eset = exp.sub, 
  probesets = exp.probes, use = "pairwise.complete", resid.only = TRUE)
exp.mat <- t(exp.mat)

# Correlate expression/methylation residuals ------------------------------

cors.df <- foreach(g = genes, .combine = "rbind") %dopar% {

  ep <- rownames(fData(exp.sub)[fData(exp.sub)$symbol == g,])
  mp <- rownames(fData(meth.sub)[fData(meth.sub)$symbol == g,])
  
  out <- cor(meth.mat[, mp], exp.mat[, ep], use = "pairwise.complete")
  
  if(length(mp) == 1 | length(ep) == 1) {
    out <- data.frame(Var1 = mp, Var2 = ep, Freq = as.numeric(out))
  } else {
    out <- as.data.frame.table(out)  
  }
  
  data.frame(symbol = g, out)
}

# Rename and reoder
names(cors.df) <- c("symbol", "meth", "exp", "r", "n")
cors.df <- cors.df[, c("symbol", "meth", "exp", "n", "r")]

# Calculate p-values
cors.df$pvalue <- WGCNA::corPvalueStudent(cors.df$r, cors.df$n)

# Calculate q-values
cors.df$qvalue <- p.adjust(cors.df$pvalue, method = "BH")

write.csv(cors.df, file = file.path("results",
  paste("crossdata-meth-exp-correlations-", tolower(tissue), ".csv", sep = "")),
  row.names = FALSE)
