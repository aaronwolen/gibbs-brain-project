# For each gene calculate correlations between its mRNA levels and the 
# methylation levels of nearby loci.
###########################################################################

library(Biobase)
source("functions/cor_cross_eset.r")
options(stringsAsFactors = FALSE)


# Load adjusted data ------------------------------------------------------
load("data/crossdata-correlations-meth-exp/eset-exp-adjusted.rda")
load("data/crossdata-correlations-meth-exp/eset-meth-adjusted.rda")


# Construct cross-platform probe look-up table ----------------------------

# Identify meth and exp probes with common targets 
meth.key <- with(fData(eset.meth), data.frame(target = symbol, meth = id))
exp.key <- with(fData(eset.exp), data.frame(target = symbol, exp = id))
probe.key <- merge(meth.key, exp.key, by = "target", sort = FALSE)


# Correlate expression/methylation residuals for each tissue --------------

cors.df <- list()

for(t in levels(eset.meth$tissue)) {

  # Tissue-specific eSet objects
  meth.sub <- eset.meth[, eset.meth$tissue == t]
  exp.sub <- eset.exp[, eset.exp$tissue == t]

  # Swap geo ids for individual ids
  sampleNames(meth.sub) <- as.character(meth.sub$individual)
  sampleNames(exp.sub) <- as.character(exp.sub$individual)
  
  # Identify samples present in both datasets
  common.ids <- intersect(sampleNames(meth.sub), sampleNames(exp.sub))
  
  meth.sub <- meth.sub[, common.ids]
  exp.sub <- exp.sub[, common.ids]

  cors.df[[t]] <- cor_cross_eset(meth.sub, exp.sub, probe.key)
}

# Combine results
cors.df <- reshape2::melt(cors.df, id.vars = colnames(cors.df[[1]]))
cors.df <- cbind(tissue = cors.df$L1, subset(cors.df, select = -L1))


# Export results ----------------------------------------------------------
write.csv(cors.df, "results/crossdata-meth-exp-correlations.csv", 
  row.names = FALSE)