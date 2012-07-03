# Within regions, which probes are differentially methylated bentween gender?
############################################################################
library(genefilter)
source("functions/lm_eset.r")

# Load data ---------------------------------------------------------------
source("4b-load-meth-data.r")


# Remove sex chromosomes --------------------------------------------------
eset.meth <- eset.meth[fData(eset.meth)$chr %in% 1:22, ]


# Filter probe-sets -------------------------------------------------------

# Identify highly variable probes
fData(eset.meth)$sd <- rowSds(exprs(eset.meth), na.rm = TRUE)

par(mfrow = c(2, 1))
hist(fData(eset.meth)$sd)
boxplot(fData(eset.meth)$sd, horizontal = TRUE)

test.probes <- featureNames(eset.meth)[fData(eset.meth)$sd > 0.1]


# Correct expression for technical variables ------------------------------

# Setting lm_eset's use argument = "pairwise.complete" causes it to loop through
# probes so missing values are be handled on a case by case basis. The result
# is a matrix of residuals.

adj.exp <- list()
for(t in levels(eset.meth$tissue)) {
  adj.exp[[t]] <- lm_eset(~ age + tissuebank, eset = eset.meth[,eset.meth$tissue == t], 
  probesets = test.probes, use = "pairwise.complete", resid.only = TRUE)
}

adj.exp <- do.call("cbind", adj.exp)

# ExpressionSet to hold lm adjusted values
eset.adj <- eset.meth[test.probes, ]
exprs(eset.adj) <- adj.exp

# Remove columns that are completely NA due to samples missing phenotypes
eset.adj <- eset.adj[, !sampleNames(eset.adj) %in% 
  names(which(apply(exprs(eset.adj), 2, function(x) sum(!is.na(x))) == 0))]

# Differential methylation - ttest ----------------------------------------

t.results <- list()

for(t in levels(eset.adj$tissue)) {
  t.results[[t]] <- rowttests(eset.adj[, eset.adj$tissue == t], "gender")
  t.results[[t]] <- data.frame(tissue = t, id = rownames(t.results[[t]]),
                               t.results[[t]])
}

t.results <- data.frame(do.call("rbind", t.results), row.names = NULL)


# Perform t-test individiually on probes with missing values --------------

na.ts <- which(is.na(t.results$statistic))

for(i in na.ts) {
  t <- t.results$tissue[i]
  id <- t.results$id[i]
  eset.sub <- eset.adj[, eset.adj$tissue == t]
  
  t.data <- data.frame(exp = exprs(eset.sub)[id,], gender = eset.sub$gender)
  t.data <- t.data[!(is.na(t.data$exp) | is.na(t.data$gender)),]
  
  t.out <- rowttests(matrix(t.data$exp, nrow = 1), fac = t.data$gender)
  t.results[i, c("statistic", "dm", "p.value")] <- t.out
}



# Correct for multiple testing --------------------------------------------

t.results$q.value <- unlist(
  with(t.results, tapply(p.value, tissue, p.adjust, method = "BH")))

with(t.results, table(tissue[p.value < .05]))
with(t.results, table(tissue[q.value < .1]))

# Annotate t-test results -------------------------------------------------

t.results <- merge(
  fData(eset.meth)[, c("id", "symbol", "chr", "position", "strand", "product")],
  t.results, by = "id")


# Export t-test results ---------------------------------------------------

# Reorder
t.results <- t.results[order(t.results$p.value),]

write.csv(t.results, "results/gender-based-differential-methylation.csv",
  row.names = FALSE)

