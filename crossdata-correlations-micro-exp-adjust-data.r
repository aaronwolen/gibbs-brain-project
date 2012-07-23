# Adjust data for microRNA/mRNA expression correlation analysis.
###########################################################################

library(Biobase)
source("functions/lm_eset.r")

# Load data ---------------------------------------------------------------
source("4a-load-exp-data.r")
source("4c-load-micro-data.r")


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


# Adjust microRNA data ----------------------------------------------------
tissues <- levels(eset.micro$tissue)

micro.mat <- list()

for(t in tissues) {
  micro.mat[[t]] <- lm_eset(~age + tissuebank, 
    eset = eset.micro[, eset.micro$tissue == t], 
    use = "pairwise.complete", resid.only = TRUE)
}

micro.mat <- do.call("cbind", micro.mat)

stopifnot(all(rownames(micro.mat) == rownames(exprs(eset.micro))))
stopifnot(all(colnames(micro.mat) == colnames(exprs(eset.micro))))

exprs(eset.micro) <- micro.mat

eset.micro@annotation <- c(eset.micro@annotation,
  "Adjusted data for covariates age and tissuebank")


# Save eSet objects as compressed rda files -------------------------------
data.dir <- "data/crossdata-correlations-micro-exp"
dir.create(data.dir)

save(list = "eset.exp", file = file.path(data.dir, "eset-exp-adjusted.rda"))
save(list = "eset.micro", file = file.path(data.dir, "eset-micro-adjusted.rda"))

