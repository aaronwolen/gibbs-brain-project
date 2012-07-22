library(Biobase)

# Load microRNA data ------------------------------------------------------

load("data/eset-micro.rda")

# Filter samples based on qc analysis -------------------------------------

# Remove samples median-ref correlations 2 SD's below the mean
ref.cors <- read.csv("figures/qc-microrna/microrna-median-ref-correlations.csv",
                     stringsAsFactors = FALSE)

outlier <- sampleNames(eset.micro) %in% ref.cors$sample[ref.cors$outlier]
cat("Removed", sum(outlier),
    "samples deemed outliers based on their correlation with tissue-median.\n")

eset.micro <- eset.micro[, !outlier]
rm(list = c("ref.cors", "outlier"))                   


# Log transformation ------------------------------------------------------

cat("Log transformed.\n")
exprs(eset.micro) <- log2(exprs(eset.micro))
