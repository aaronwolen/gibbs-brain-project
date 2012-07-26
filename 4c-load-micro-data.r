library(Biobase)

# Load microRNA data ------------------------------------------------------

load("data/eset-micro.rda")

# Filter samples based on qc analysis -------------------------------------

# Remove samples squared median deviations 2 SD's above the mean
ref.comps <- read.csv("figures/qc-microrna/microrna-median-ref-comparisons.csv",
                     stringsAsFactors = FALSE)

outlier <- sampleNames(eset.micro) %in% ref.comps$sample[ref.comps$outlier]
cat("Removed", sum(outlier),
    "samples deemed outliers based on their deviations from the tissue-median.\n")

eset.micro <- eset.micro[, !outlier]
rm(list = c("ref.comps", "outlier"))                   


# Log transformation ------------------------------------------------------

cat("Log transformed.\n")
exprs(eset.micro) <- log2(exprs(eset.micro))
