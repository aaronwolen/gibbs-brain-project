# Calculate cross-tissue correlations for each CpG probe.
###########################################################################

library(Biobase)
library(reshape2)

# Load adjusted data ------------------------------------------------------
load("data/crosstissue-correlations-meth/eset-meth-adjusted.rda")


# Remove sex chromosomes --------------------------------------------------
eset.meth <- eset.meth[fData(eset.meth)$chr %in% 1:22, ]


# Filter probe-sets -------------------------------------------------------

# Identify highly variable probes
fData(eset.meth)$sd <- genefilter::rowSds(exprs(eset.meth), na.rm = TRUE)

par(mfrow = c(2, 1))
hist(fData(eset.meth)$sd)
boxplot(fData(eset.meth)$sd, horizontal = TRUE)

probes <- featureNames(eset.meth)[fData(eset.meth)$sd > 0.05]


# Loop through each probe and calculate correlations ----------------------

pb <- txtProgressBar(max = length(probes), style = 3)

cors.df <- list()

for(p in probes) {
  setTxtProgressBar(pb, value = which(probes == p))
  # Assemble and reshape methylation data for current probe
  p.df <- as.data.frame.table(exprs(eset.meth)[p,], stringsAsFactors = F)
  p.df <- cbind(p.df, pData(eset.meth)[p.df$Var1, c("individual", "tissue")])
  p.df <- dcast(p.df, individual ~ tissue, value.var = "Freq")
  
  # Calculate correlations and reshape results
  p.cor <- cor(as.matrix(p.df[, -1]), use = "pairwise")
  p.cor[lower.tri(p.cor, diag = TRUE)] <- NA
  p.cor <- as.data.frame.table(p.cor, stringsAsFactors = F)
  p.cor <- p.cor[!is.na(p.cor$Freq),]
  
  # Count number of complete cases for each tissue pair
  p.cor$n <- apply(p.cor, 1, function(x)
    sum(!is.na(p.df[, x["Var1"]]) & !is.na(p.df[, x["Var2"]])))
 
  cors.df[[p]] <- p.cor
}


# Reformat correlation results --------------------------------------------

# Melt list
cors.df <- melt(cors.df, id.vars = names(cors.df[[1]]))

# Combine tissue columns into a single variable
cors.df <- transform(cors.df, tissue.pair = paste(Var1, Var2, sep = "-"))
cors.df <- subset(cors.df, select = c(-Var1, -Var2))

# Rename and reoder
names(cors.df) <- c("r", "n", "id", "tissue.pair")
cors.df <- cors.df[, c("id", "tissue.pair", "r", "n")]


# Calculate p-values and FDR ----------------------------------------------

# Calculate p-values
cors.df$pvalue <- WGCNA::corPvalueStudent(cors.df$r, cors.df$n)

# Calculate q-values
cors.df <- ddply(cors.df, .(tissue.pair), transform,
  qvalue = p.adjust(pvalue, method = "BH"))


# Export results ----------------------------------------------------------
write.csv(cors.df, file = file.path("results",
  "crosstissue-correlations-meth.csv"), row.names = FALSE)
