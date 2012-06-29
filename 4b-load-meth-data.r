library(Biobase)

# Load methylation data ---------------------------------------------------

load("data/eset-meth.rda")


# Filter samples based on qc analysis -------------------------------------

# Remove samples with an excessive number of missing probes.
# (meth-qc_NA-probes-histogram.pdf)

# Count NAs per sample
eset.meth$na.count <- apply(exprs(eset.meth), 2, function(x) sum(is.na(x)))

na.filter <- eset.meth$na.count < 50
cat("Removed", sum(!na.filter),
    "samples with an excessive number of missing probes.\n")
eset.meth <- eset.meth[, na.filter]
rm(na.filter)                    


# Remove outlier samples with median intensities > 0.2
# (meth-qc_histogram.pdf, meth-qc_boxplot.pdf)
sample.mdn <- rowMedians(t(exprs(eset.meth)), na.rm = TRUE)
mdn.filter <- sample.mdn < 0.2

cat("Removed", sum(!mdn.filter),
    "samples with excessively high median intensities.\n")

eset.meth <- eset.meth[, mdn.filter]
rm(list = c("sample.mdn", "mdn.filter"))
