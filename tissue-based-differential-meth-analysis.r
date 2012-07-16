# Within regions, which probes are differentially methylated bentween gender?
############################################################################
library(genefilter)
source("functions/lm_eset.r")

# Nice html table for wiki
print_xtable <- function(x) {
  print(xtable::xtable(x), type = "html", include.colnames = FALSE)
}

# Load data ---------------------------------------------------------------
source("4b-load-meth-data.r")


# Filter probe-sets -------------------------------------------------------

# Identify highly variable probes
fData(eset.meth)$sd <- rowSds(exprs(eset.meth), na.rm = TRUE)

par(mfrow = c(2, 1))
hist(fData(eset.meth)$sd)
boxplot(fData(eset.meth)$sd, horizontal = TRUE)

test.probes <- featureNames(eset.meth)[fData(eset.meth)$sd > 0.1]

# Differential methylation - ttest ----------------------------------------

tissues <- levels(eset.meth$tissue)

# All combinations of tissue pairs
pairs <- t(combn(tissues, 2))

t.results <- list()

for(i in seq_len(nrow(pairs))) {
  pair.label <- paste(pairs[i, ], collapse = "-")
  
  # Subset for current tissue pair
  meth.sub <- eset.meth[, eset.meth$tissue %in% pairs[i,]]
  meth.sub$tissue <- factor(meth.sub$tissue)
  
  t.results[[pair.label]] <- rowttests(meth.sub, "tissue")
  t.results[[pair.label]] <- data.frame(pair = pair.label, 
    id = rownames(t.results[[i]]), t.results[[i]])
}

t.results <- data.frame(do.call("rbind", t.results), row.names = NULL)


# Correct for multiple testing --------------------------------------------

t.results$q.value <- unlist(
  with(t.results, tapply(p.value, pair, p.adjust, method = "BH")))

print_xtable(with(t.results, table(pair[p.value < .05])))
print_xtable(with(t.results, table(pair[q.value < .1])))

# Annotate t-test results -------------------------------------------------

t.results <- merge(
  fData(eset.meth)[, c("id", "symbol", "chr", "position", "strand", "description")],
  t.results, by = "id")


# Export t-test results ---------------------------------------------------

# Reorder
t.results <- t.results[order(t.results$p.value),]

write.csv(t.results, "results/tissue-based-differential-methylation.csv",
  row.names = FALSE)

