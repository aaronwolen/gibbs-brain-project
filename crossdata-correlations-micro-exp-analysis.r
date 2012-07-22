# For each microRNA calculate correlations with all genes annotated as putative
# targets by mirBase. This script should be batch processed and 'tissue' should
# be supplied as a command line argument equal to CRBLM, PONS, FCTX or TCTX.
###########################################################################

library(Biobase)
library(RmiR.hsa)

library(foreach)
library(doMC)

source("functions/lm_eset.r")

# Read in command line arguments ------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Evaluate the expressions containing each argument
if(length(args) == 0){
  stop("tissue argument not supplied.")
} else {
  for(i in seq_along(args)){
    eval(parse(text = args[[i]]))
  }
}

tissue <- match.arg(tissue, c("PONS", "CRBLM", "FCTX", "TCTX"))

# Register multicore backend ----------------------------------------------
n.cores <- multicore:::detectCores()
registerDoMC(n.cores)


# Load data ---------------------------------------------------------------
source("4a-load-exp-data.r")
source("4c-load-micro-data.r")



# Identify  miroRNA targets -----------------------------------------------

# Extract mirbase database
mirbase <- dbGetQuery(RmiR.hsa_dbconn(), 
  "SELECT * FROM mirbase where mature_miRNA like '%'")

# Identify microRNA ID's in eSet and mirbase
mirs <- intersect(fData(eset.micro)$symbol, mirbase$mature_miRNA)

mirbase <- subset(mirbase, mature_miRNA %in% mirs)


# Subset ExpressionSets ---------------------------------------------------

# Identify individuals present in both datasets
common.ids <- intersect(eset.exp[, eset.exp$tissue == tissue]$individual,
  eset.micro[, eset.micro$tissue == tissue]$individual)

cat(length(common.ids), "individual IDs present in both datasets.\n")

# Region-specific eSet objects
exp.sub <- eset.exp[, eset.exp$individual %in% common.ids & 
                      eset.exp$tissue == tissue]

micro.sub <- eset.micro[, eset.micro$individual %in% common.ids &
                      eset.micro$tissue == tissue]

# Swap geo ids for individual ids
sampleNames(exp.sub) <- exp.sub$individual
sampleNames(micro.sub) <- micro.sub$individual

# Reorder second dataset to match first
micro.sub <- micro.sub[, sampleNames(exp.sub)]
stopifnot(sampleNames(exp.sub) == sampleNames(micro.sub))


# Extract expression matrices ---------------------------------------------

exp.mat <- exprs(exp.sub)
micro.mat <- exprs(micro.sub)

# Identify microRNA targets and exp probes with common entrez IDs 
entrezs <- intersect(fData(exp.sub)$entrez, mirbase$gene_id)

# Identify probes corresponding to the microRNA and putative targets
micro.probes <- subset(fData(micro.sub), symbol %in% mirs)$id
exp.probes <- subset(fData(exp.sub), entrez %in% entrezs)$id

# Correct expression for technical variables ------------------------------

# Setting lm_eset's use argument = "pairwise.complete" causes it to loop through
# probes so missing values are be handled on a case by case basis. The result
# is a matrix of residuals.

# micro.mat <- lm_eset(~ age + tissuebank, eset = micro.sub, 
#   probesets = micro.probes, use = "pairwise.complete", resid.only = TRUE)
micro.mat <- t(micro.mat)
  
exp.mat <- lm_eset(~ age + tissuebank, eset = exp.sub, 
  probesets = exp.probes, use = "pairwise.complete", resid.only = TRUE)
exp.mat <- t(exp.mat)

# Correlate expression/methylation residuals ------------------------------

cors.df <- foreach(m = mirs, .combine = "rbind") %dopar% {

  # MicroRNA probe
  mp <- subset(fData(eset.micro), symbol == m)$id
  
  # Putative targets of m
  m.genes <- subset(mirbase, mature_miRNA == m)$gene_id
  
  # Probes for putative targets
  ep <- subset(fData(exp.sub), entrez %in% m.genes)$id
  symbols <- subset(fData(exp.sub), entrez %in% m.genes)$symbol

  mp.mat <- matrix(micro.mat[, mp], ncol = 1, 
    dimnames = list(rownames(micro.mat), mp))
  
  out <- cor(mp.mat, exp.mat[, ep], use = "pairwise.complete")
  
  if(length(ep) == 1) {
    out <- data.frame(Var1 = mp, Var2 = ep, Freq = as.numeric(out))
  } else {
    out <- as.data.frame.table(out)  
  }
  
  # Count number of values used in each correlation
  out$n <- apply(as.matrix(out), 1, function(x)
    sum(!is.na(micro.mat[, x[1]]) & !is.na(exp.mat[, x[2]])))
  
  data.frame(mir = m, symbol = symbols, out)
}

# Rename and reoder
names(cors.df) <- c("micorna", "symbol", "micro", "exp", "r", "n")
cors.df <- cors.df[, c("micorna", "symbol", "micro", "exp", "n", "r")]

# Calculate p-values
cors.df$pvalue <- WGCNA::corPvalueStudent(cors.df$r, cors.df$n)

# Calculate q-values
cors.df$qvalue <- p.adjust(cors.df$pvalue, method = "BH")

write.csv(cors.df, file = file.path("results",
  paste("crossdata-micro-exp-correlations-", tolower(tissue), 
  ".csv", sep = "")), row.names = FALSE)
