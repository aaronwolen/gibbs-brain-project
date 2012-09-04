# Within each tissue, calculate correlations between each microRNA and its
# putative, based on various target prediction databases.
###########################################################################

library(Biobase)
library(RmiR.hsa)
library(org.Hs.eg.db)
source("functions/cor_cross_eset.r")
options(stringsAsFactors = FALSE)

source("functions/eval-args.r")

# Load adjusted data ------------------------------------------------------
load("data/crossdata-correlations-micro-exp/eset-exp-adjusted.rda")
load("data/crossdata-correlations-micro-exp/eset-micro-adjusted.rda")


# Create look-up table for microRNAs and their targets --------------------

if( !exists("mirdb") ) {
  mirdb <- "tarbase"
}

# Check validity of mirdb
mirdb <- match.arg(mirdb, dbListTables(RmiR.hsa_dbconn()))

cat("Performing correlations between microRNAs and targets predicted by", 
  mirdb, "\n")

# Extract mirdb database table
mir.key <- dbGetQuery(RmiR.hsa_dbconn(), 
  paste("SELECT * FROM", mirdb, "where mature_miRNA like '%'"))
mir.key$gene_id <- as.character(mir.key$gene_id)
mir.key <- subset(mir.key, select = c(mature_miRNA, gene_id))

# Add microRNA probes
mir.key <- merge(mir.key, fData(eset.micro)[, c("id", "symbol")],
  by.x = "mature_miRNA", by.y = "symbol", sort = FALSE)
names(mir.key)[ncol(mir.key)] <- "micro"

# Add mRNA probes based on entrez mappings
mir.key <- merge(mir.key, fData(eset.exp)[, c("id", "entrez")],
  by.x = "gene_id", by.y = "entrez", sort = FALSE)
names(mir.key)[ncol(mir.key)] <- "exp"

mir.key <- subset(mir.key, select = -gene_id)
names(mir.key)[1] <- "microrna"

# Correlate microRNA/mRNA residuals for each tissue -----------------------

cors.df <- list()

for(t in levels(eset.micro$tissue)) {
  
  # Tissue-specific eSet objects
  micro.sub <- eset.micro[, eset.micro$tissue == t]
  exp.sub <- eset.exp[, eset.exp$tissue == t]
  
  # Swap geo ids for individual ids
  sampleNames(micro.sub) <- as.character(micro.sub$individual)
  sampleNames(exp.sub) <- as.character(exp.sub$individual)
  
  # Identify individuals present in both datasets
  common.ids <- intersect(sampleNames(micro.sub), sampleNames(exp.sub))
  
  micro.sub <- micro.sub[, common.ids]  
  exp.sub <- exp.sub[, common.ids]
  
  cors.df[[t]] <- cor_cross_eset(micro.sub, exp.sub, mir.key)
}

# Combine results
cors.df <- reshape2::melt(cors.df, id.vars = colnames(cors.df[[1]]))
cors.df <- cbind(tissue = cors.df$L1, subset(cors.df, select = -L1))


# Export results ----------------------------------------------------------
save(list = "cors.df", file = file.path("results",
  paste("crossdata-micro-exp-correlations-", mirdb, ".rda", sep = "")))

