library(affy)
library(ggplot2)
library(ggaffy)

source("functions/gg_save.r")


# Graphing variables ------------------------------------------------------

log.label <- expression(log[2] * " expression")

# Shrink text so labels are legible along x-axis
axis.text.x <- theme_text(size = 7, colour = "grey50", angle = 90, hjust = 1)


# Load data ---------------------------------------------------------------

load("data/eset-exp.rda")


# Sample filters ----------------------------------------------------------

# Should samples with negative expression values be removed?

# Count negative probes per sample
eset.exp$neg.count <- apply(exprs(eset.exp), 2, function(x) sum(x < 0))

# A small number of samples have excessive negative values
table(eset.exp$neg.count)

ggplot(pData(eset.exp)) + 
  aes(x = individual, y = neg.count, fill = batch) +
  geom_bar() + facet_wrap(~tissue, ncol = 1) +
  ylab("Negative probe count") +
  scale_fill_brewer(palette = "Set1") +
  opts(legend.position = "top", axis.text.x = axis.text.x)

gg_save("exp-neg-probes-per-individual.pdf", 
  path = "figures/qc", height = 9, width = 13)



# Log transform -----------------------------------------------------------

exprs(eset.exp) <- log2(exprs(eset.exp))


# Distributions of unnormalized data --------------------------------------

# Quantiles plot
ggaffy_quantiles(eset.exp,
  x = "individual", probs = c(0.01, 0.1, 0.25),
  fill = "batch", shape = "gender", facet.row = "tissue") +
  ylab(log.label) + opts(axis.text.x = axis.text.x)
gg_save("exp-quantiles.pdf", path = "figures/qc", width = 13, height = 9)

# Boxplot
ggaffy_box(eset.exp, x = "individual", fill = "batch", facet.row = "tissue") +
  ylab(log.label) + opts(axis.text.x = axis.text.x)
gg_save("exp-boxplot.pdf", path = "figures/qc", width = 13, height = 9)

ggaffy_box(eset.exp, x = "individual", method = "rle",
  fill = "batch", facet.row = "tissue") + opts(axis.text.x = axis.text.x)
gg_save("exp-boxplot-rle.pdf", path = "figures/qc", width = 13, height = 9)

# Histogram
ggaffy_hist(eset.exp, facet.row = "tissue", color = "batch") + xlab(log.label)
gg_save("exp-histogram.pdf", path = "figures/qc", width = 13, height = 9)
  

