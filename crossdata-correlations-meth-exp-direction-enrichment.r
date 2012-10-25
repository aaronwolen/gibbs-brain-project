# Are positive methylation/mrna correlations more likely to involve CpG islands
# than what would be expected by chance?

library(Biobase)
library(plyr)
library(ggplot2)
library(reshape2)

# Load and prepare data ---------------------------------------------------

load("data/eset-meth.rda")

cors.df <- read.csv("results/crossdata-meth-exp-correlations.csv",
  stringsAsFactors = FALSE)

# Classify correlations
cors.df$direction <- factor(ifelse(cors.df$r > 0, "positive", "negative"))

# Annotate probes in correlation results
cors.df <- data.frame(cors.df,
  fData(eset.meth)[cors.df$meth, 
  c("chr", "position", "strand", "cpg.island", "distance.to.tss", "description")],
  row.names = NULL)



# Enrichment analysis -----------------------------------------------------
sig.cors <- subset(cors.df, qvalue < 0.1)

ddply(sig.cors, .(tissue), summarise,
      positive.cor = sum(r > 0),
      negative.cor = sum(r < 0))

counts.df <- ddply(sig.cors, .(tissue, cpg.island), summarise,
      positive = sum(r > 0),
      negative = sum(r < 0))
counts.df

# Counts
counts.table <- with(sig.cors, table(cpg.island, direction))
addmargins(counts.table)

# Proportions
round(addmargins(prop.table(counts.table)), d = 2)

# Chi-squared test
chisq.test(counts.table)
fisher.test(counts.table, alternative = "greater")


# Visualize enrichment ----------------------------------------------------
plot.df <- melt(counts.df, measure.vars = c("positive", "negative"),
                variable.name = "direction", value.name = "count")
ggplot(plot.df) + 
  aes(x = cpg.island, y = count, fill = direction) + 
  geom_bar(position = position_fill()) +
  facet_grid(. ~ tissue, scales = "free") +
  scale_y_continuous(labels = scales::percent_format())