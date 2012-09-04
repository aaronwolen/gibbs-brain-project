### <% mirdb %> predicted interactions
* results: [`crossdata-micro-exp-correlations-<% mirdb %>.rda`](https://github.com/aaronwolen/gibbs-brain-project/blob/master/results/crossdata-micro-exp-correlations-<% mirdb %>.rda)


```{r load-<% mirdb %>-results, cache=FALSE, message=FALSE}
load("results/crossdata-micro-exp-correlations-<% mirdb %>.rda")


# Annotate probes
cors.df <- data.frame(cors.df,
  fData(eset.exp)[cors.df$exp, c("symbol", "entrez")],
  fData(eset.micro)[cors.df$micro, c("chr", "position", "strand")],
  row.names = NULL)
  
print(head(cors.df))
```


```{r <% mirdb %>-results-table, results='asis'}
results.table <- ddply(cors.df, .(tissue), summarise, 
  'p-value < 0.05' = sum(pvalue < 0.05),
  'q-value < 0.05' = sum(adj.pvalue < 0.05))

print_xtable(results.table)
```

#### Correlation and p-value distributions
Histograms of the Pearson correlation values (left) and corresponding p-values (right).

```{r <% mirdb %>-histograms, message=FALSE}
p.hist <- ggplot(cors.df) + 
  aes(x = pvalue) + geom_histogram() +
  facet_grid(tissue~.)

r.hist <- ggplot(cors.df) + 
  aes(x = r) + geom_histogram() +
  geom_vline(xintercept = 0, color = "white") +
  facet_grid(tissue~.)

ggaffy:::ggaffy_layout(nrow = 1, ncol = 2, plot.list = list(r.hist, p.hist))
```


#### Frequency of significant microRNA correleations
A subset of microRNAs are correlated with very large lists of genes. The following figure summarizes the number of significant correlations associated with individual microRNAs. The tables below provide for each tissue the ten microRNAs associated with the largest number of significant correlations. 

```{r <% mirdb %>-count-correlations, fig.width=5, fig.height=3}
sig.cors <- subset(cors.df, adj.pvalue < .05)
sig.cors$tissue <- factor(sig.cors$tissue)
sig.cors <- arrange(sig.cors, adj.pvalue)

# Count significant cors associated with each microRNA                  
mir.counts <- with(sig.cors, tapply(microrna, tissue, table))                  
                  
mir.counts <- lapply(mir.counts, as.data.frame.table)
mir.counts <- reshape2::melt(mir.counts, id.vars = c("Var1", "Freq"))
names(mir.counts) <- c("microrna", "count", "tissue")

mir.counts <- subset(mir.counts, count > 0)
                  
mir.counts <- arrange(mir.counts, tissue, desc(count))

ggplot(mir.counts) + aes(x = tissue, y = count) + geom_boxplot() +
  ylab("microRNA correlations count (q-value < 0.05)")
```

```{r <% mirdb %>-count-correlations-tables, results='asis'}
tissues <- c("CRBLM" = "Cerebellum", "FCTX" = "Frontal cortex", 
   "TCTX" = "Temporal cortex", "PONS" = "Pons")

for(t in unique(mir.counts$tissue)) {
  cur.t <- tissues[t]
  tab <- subset(mir.counts, tissue == t, select = -tissue)
  
  if(nrow(tab) > 10) tab <- tab[1:10,]
  
  cat("\n#####", tissues[t], "\n")
  print_xtable(tab)
}
```


#### Positive vs negative correlations
In all tissues the number of positive correlations exceeded that of negative correlations. 

```{r <% mirdb %>-count-correlations-direction, fig.width=5, fig.height=3}
ggplot(data.frame(sig.cors, 
  dir = ifelse(sig.cors$r > 0, "positive", "negative"))) + 
  geom_bar(aes(x = tissue, fill = dir), position = position_dodge()) +
  ylab("Correlation count (q-value < 0.05)") + 
  guides(fill = guide_legend(title = "correlation\ndirection"))
```

#### Scatter plots
Scatter plots for the 50 most significant microRNA/mRNA correlations are provided [here](https://github.com/aaronwolen/gibbs-brain-project/blob/master/figures/crossdata-correlations-micro-exp/<% mirdb %>-cors-micro-exp-top-probes.pdf).

```{r <% mirdb %>-sig-probes-scatterplots, include=FALSE}
# Produce visualizations for individual probes with largest correlations
nplots <- ifelse(nrow(sig.cors) > 50, 50, nrow(sig.cors))

# Phenotype data for plots
pdata.cols <- c("geo", "individual", "gender")
pdata <- rbind(pData(eset.exp)[, pdata.cols], pData(eset.micro)[, pdata.cols])

# Vector of microRNA/mRNA probe ids
sig.cors$id <- with(sig.cors, paste(exp, micro, sep = "-"))

probe.plots <- list(); i <- 0

while(length(probe.plots) < nplots) {  
  
  i <- i + 1
  cur.id  <- sig.cors$id[i]
  
  # Skip current correlation is its been previously plotted
  if(cur.id %in% names(probe.plots)) {
    next()
  }
  
  m.probe <- as.character(sig.cors$micro[i])
  e.probe <- as.character(sig.cors$exp[i])

  probe.df <- data.frame(sample = sampleNames(eset.exp),
    tissue = eset.exp$tissue,
    expression = exprs(eset.exp)[e.probe,], assay = "mRNA", 
    row.names = NULL)
  
  probe.df <- rbind(probe.df,
    data.frame(sample = sampleNames(eset.micro),
      tissue = eset.micro$tissue,
      expression = exprs(eset.micro)[m.probe,], assay = "microRNA",
      row.names = NULL))
   
  probe.df <- merge(probe.df, pdata, by.x = "sample", by.y = "geo") 
    
  probe.df <- reshape2::dcast(probe.df, 
    individual + gender + tissue ~ assay, value.var = "expression")
  
  probe.plots[[cur.id]] <- ggplot(probe.df) +
    aes(x = microRNA, y = mRNA) +
    geom_point(aes(color = gender), alpha = 0.5) +
    geom_text(aes(label = individual), color = "white", size = 0.5) +
    geom_rug(alpha = 0.5) +
    facet_wrap(~tissue, ncol = 2)  +
    stat_smooth(method = "lm") +
    xlab(paste("microRNA expression", " (", m.probe, ")", sep = "")) +
    ylab(paste("mRNA expression", " (",  e.probe, ")", sep = "")) +
    opts(title = paste(sig.cors$symbol[i], "/", sig.cors$microrna[i], 
      " (", tolower(sig.cors$tissue[i]), 
      " cor = ", round(sig.cors$r[i], 2), ")", sep = ""))
}

pdf(file.path(fig.dir, paste(fig.prefix, "<% mirdb %>-top-probes.pdf", sep = "")),
    width = 7, height = 6)
temp <- lapply(probe.plots, print)
dev.off()
```
