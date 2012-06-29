library(affy)
library(ggplot2)
library(ggaffy)

source("functions/gg_save.r")


# Graphing variables ------------------------------------------------------

# Shrink text so labels are legible along x-axis
axis.text.x <- theme_text(size = 7, colour = "grey50", angle = 90, hjust = 1)

# Load data ---------------------------------------------------------------

load("data/eset-meth.rda")


# Sample filters ----------------------------------------------------------

# Should samples with high NA counts be removed?

# Count NAs per sample
eset.meth$na.count <- apply(exprs(eset.meth), 2, function(x) sum(is.na(x)))

qplot(eset.meth$na.count) + 
  ylab("Sample count") + xlab("NA methylation probe count")
gg_save("meth-NA-probes-histogram.pdf", path = "figures/qc")

ggplot(pData(eset.meth)) + 
  aes(x = individual, y = na.count, fill = batch) +
  geom_bar() + facet_wrap(~tissue, ncol = 1) +
  ylab("NA methylation probe count") +
  scale_fill_brewer(palette = "Set1") +
  opts(legend.position = "top", axis.text.x = axis.text.x)

gg_save("meth-NA-probes-per-individual.pdf", 
  height = 9, width = 13, path = "figures/qc")


# Count NAs per probe
# A few red probes stand out as having a high number of missing values
fData(eset.meth)$na.count <- apply(exprs(eset.meth), 1, function(x) sum(is.na(x)))

ggplot(fData(eset.meth)) + aes(x = color_channel, y = na.count) + geom_point()
head(plyr::arrange(fData(eset.meth), plyr::desc(na.count)), 10)



# Distributions of unnormalized data --------------------------------------

# Quantiles plot
ggaffy_quantiles(eset.meth,
  x = "individual", probs = c(0.01, 0.1, 0.25),
  fill = "batch", shape = "gender", facet.row = "tissue") +
  opts(axis.text.x = axis.text.x)
gg_save("meth-quantiles.pdf", path = "figures/qc", width = 13, height = 9)

# Boxplot
ggaffy_box(eset.meth, x = "individual", fill = "batch", facet.row = "tissue") +
  opts(axis.text.x = axis.text.x)
gg_save("meth-boxplot.pdf", path = "figures/qc", width = 13, height = 9)

ggaffy_box(eset.meth, x = "individual", method = "rle",
  fill = "batch", facet.row = "tissue") + opts(axis.text.x = axis.text.x)
gg_save("meth-boxplot-rle.pdf", path = "figures/qc", width = 13, height = 9)

# Histogram
ggaffy_hist(eset.meth, facet.row = "tissue", color = "batch", adjust = 0.5)
gg_save("meth-histogram.pdf", path = "figures/qc", width = 9, height = 13)


