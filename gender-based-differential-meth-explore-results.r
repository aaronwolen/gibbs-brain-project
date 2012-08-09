# Explore and visualize results from the gender-based differential methylation
# analysis performed in gender-based-differential-meth.r
###########################################################################

library(plyr)
library(reshape2)
library(ggplot2)

# ggpairs function
devtools::source_gist("3298388")  

# Load data ---------------------------------------------------------------

source("4b-load-meth-data.r")
t.results <- read.csv("results/gender-based-differential-methylation.csv")


# Figure variables --------------------------------------------------------

fig.dir <- "figures/gender-based-meth-differences"
dir.create(fig.dir)

fig.prefix <- "meth-gender-"

source("functions/save-figures.r")


# Volcano plot ------------------------------------------------------------

vp <- ggplot(t.results) +
  aes(x = dm, y = -log10(p.value)) +
  geom_point(size = 1) +
  facet_wrap(~tissue) +
  xlab("Mean gender difference") + 
  ylab(expression(-log[10] * " p-value"))

# Remove TLE1 outlier
vp %+% subset(t.results, id != "cg15915418")
save_plot("volcano-plot.png")


# Pairwise tissue comparisons of mean gender differences ------------------
t.dm <- dcast(t.results, id ~ tissue, value.var = "dm")

ggpairs(data = subset(t.dm, id != "cg15915418"), 
  facet.vars = c("CRBLM", "FCTX", "PONS", "TCTX"), alpha = 0.5, size = 0.75) + 
  stat_smooth(method = "lm") + 
  xlab("Mean gender difference")  + ylab("Mean gender difference")
              

save_plot("dm-tissue-pairs.png", width = 5.5)

