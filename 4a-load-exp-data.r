library(Biobase)

# Load data ---------------------------------------------------------------

load("data/eset-exp.rda")


# Sample filters ----------------------------------------------------------

# Remove samples with an excessive number of negative probes.
# (exp-neg-probes-per-individual.pdf)

# Count negative probes per sample
eset.exp$neg.count <- apply(exprs(eset.exp), 2, function(x) sum(x < 0))

neg.filter <- eset.exp$neg.count == 0
cat("Removed", sum(!neg.filter), "samples with negative probes.\n")
eset.exp <- eset.exp[, neg.filter]
rm(neg.filter)



# Log transformation ------------------------------------------------------

cat("Log transformed.\n")
exprs(eset.exp) <- log2(exprs(eset.exp))


# Quantile normalization --------------------------------------------------

cat("Quantile normalized.\n")

exprs(eset.exp) <- matrix(broman::normalize(exprs(eset.exp)), 
  ncol = ncol(eset.exp), 
  dimnames = list(featureNames(eset.exp), sampleNames(eset.exp)))
