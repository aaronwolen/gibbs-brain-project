# Adjust data for correlation analysis between methylation and mRNA levels
###########################################################################

library(Biobase)
source("functions/lm_eset.r")

# Load data ---------------------------------------------------------------
source("4a-load-exp-data.r")
source("4b-load-meth-data.r")


# Adjust expression data --------------------------------------------------
tissues <- levels(eset.exp$tissue)

exp.mat <- list()

for(t in tissues) {
  exp.mat[[t]] <- lm_eset(~age + tissuebank, 
    eset = eset.exp[, eset.exp$tissue == t], 
    use = "pairwise.complete", resid.only = TRUE)
}

exp.mat <- do.call("cbind", exp.mat)

stopifnot(all(rownames(exp.mat) == rownames(exprs(eset.exp))))
stopifnot(all(colnames(exp.mat) == colnames(exprs(eset.exp))))

exprs(eset.exp) <- exp.mat

eset.exp@annotation <- c(eset.exp@annotation,
  "Adjusted data for covariates age and tissuebank")


# Adjust methylation data -------------------------------------------------
tissues <- levels(eset.meth$tissue)

meth.mat <- list()

for(t in tissues) {
  meth.mat[[t]] <- lm_eset(~age + tissuebank, 
    eset = eset.meth[, eset.meth$tissue == t], 
    use = "pairwise.complete", resid.only = TRUE)
}

meth.mat <- do.call("cbind", meth.mat)

stopifnot(all(rownames(meth.mat) == rownames(exprs(eset.meth))))
stopifnot(all(colnames(meth.mat) == colnames(exprs(eset.meth))))

exprs(eset.meth) <- meth.mat

eset.meth@annotation <- c(eset.meth@annotation,
  "Adjusted data for covariates age and tissuebank")


# Save eSet objects as compressed rda files -------------------------------
data.dir <- "data/crossdata-correlations-meth-exp"
dir.create(data.dir)

save(list = "eset.exp", file = file.path(data.dir, "eset-exp-adjusted.rda"))
save(list = "eset.meth", file = file.path(data.dir, "eset-meth-adjusted.rda"))

