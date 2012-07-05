# Extract phenodata from each of the geo datasets and combine into a single 
# data.frame containing sample phenotypic data to be added to each expressionSet
# object's pData slot.

# Use appropriate columns from pData to constrcut MIAME objects

library(Biobase)
library(plyr)
library(ggplot2)
source("functions/cleaning-functions.r")
source("functions/load_gists.r")

# Data' -------------------------------------------------------------------
load("data/GSE15745_raw-data.rda")

# sample phenotype data from publication
sample.pdata <- read.csv("data/pm20485568_supptable_sampleinfo.csv")

# Notes -------------------------------------------------------------------

# title (renamed to individual) uniquely identifiers each individiual that
# donated tissue.

# geo_accession (renamed to sample) uniquely identifies a specific tissue from
# each individual.


# Identify platforms ------------------------------------------------------

# GPL6104 Illumina humanRef-8 v2.0 expression beadchip
# GPL8178	Illumina Human v1 MicroRNA expression beadchip
# GPL8490	Illumina HumanMethylation27 BeadChip

platforms <- data.frame(
  label = names(geo.data),
  gpl = unlist(lapply(geo.data, slot, name = "annotation")), 
  row.names = NULL)

# Look-up table for geo platform
gpl <- c("exp" = "GPL6104", "meth" = "GPL8490", "micro" = "GPL8178")

geo.exp <- geo.data[platforms$label[platforms$gpl == gpl["exp"]]]
geo.meth <- geo.data[platforms$label[platforms$gpl == gpl["meth"]]]
geo.micro <- geo.data[platforms$label[platforms$gpl == gpl["micro"]]]


# Construct pdata data.frame and MIAME objects ----------------------------

source("2a-reformat-pdata.r")
source("2b-reformat-fdata.r")


# Expression - construct eSet object --------------------------------------

# Combine into a single expression matrix
exp.mat <- do.call("cbind", lapply(geo.exp, exprs))

# Ensure fdata matches order of expression matrix
fdata.exp <- fdata.exp[rownames(exp.mat),]

# Ensure pdata matches order of expression matrix
pdata.exp <- pdata.exp[colnames(exp.mat),]

eset.exp <- ExpressionSet(exp.mat, 
  phenoData = AnnotatedDataFrame(pdata.exp), 
  featureData = AnnotatedDataFrame(fdata.exp, varMetadata = fmeta.exp),
  experimentData = miame.exp)



# Methylation - construct eSet object -------------------------------------

# Combine into a single expression matrix
meth.mat <- do.call("cbind", lapply(geo.meth, exprs))

# Ensure fdata matches order of expression matrix
fdata.meth <- fdata.meth[rownames(meth.mat),]

# Ensure pdata matches order of expression matrix
pdata.meth <- pdata.meth[colnames(meth.mat),]

eset.meth <- ExpressionSet(meth.mat, 
  phenoData = AnnotatedDataFrame(pdata.meth), 
  featureData = AnnotatedDataFrame(fdata.meth, varMetadata = fmeta.meth),
  experimentData = miame.meth)



# microRNA - construct eSet object --------------------------------------

# Combine into a single microRNA matrix
micro.mat <- do.call("cbind", lapply(geo.micro, exprs))

# Ensure fdata matches order of microRNA matrix
fdata.micro <- fdata.micro[rownames(micro.mat),]

# Ensure pdata matches order of microRNA matrix
pdata.micro <- pdata.micro[colnames(micro.mat),]

eset.micro <- ExpressionSet(micro.mat, 
  phenoData = AnnotatedDataFrame(pdata.micro), 
  featureData = AnnotatedDataFrame(fdata.micro, varMetadata = fmeta.micro),
  ExperimentData = miame.micro)



# Reorder eSet objects ----------------------------------------------------

# Order all data-set samples by tissue and then individual id
reorder_eset <- function(eset) {
  tissues <- as.character(eset$tissue)
  samples <- as.numeric(as.character(eset$individual))
  return(eset[, order(tissues, samples)])
}

eset.exp <- reorder_eset(eset.exp)
eset.meth <- reorder_eset(eset.meth)
eset.micro <- reorder_eset(eset.micro)

# Save eSet objects as compressed rda files -------------------------------

save(list = "eset.exp", file = "data/eset-exp.rda", compress = TRUE)
save(list = "eset.meth", file = "data/eset-meth.rda", compress = TRUE)
save(list = "eset.micro", file = "data/eset-micro.rda", compress = TRUE)
