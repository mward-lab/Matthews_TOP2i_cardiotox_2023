# ---
#   title: "Sources of variation in the data"
# date: 2016-07-07
# output:
#   html_document:
#   toc: true
# toc_float: true
# ---
#
#   **Last updated:** `r Sys.Date()`
#
# **Code version:** `r system("git log -1 --format='%H'", intern = TRUE)`
#
# ```{r chunk-options}
# library("knitr")
# opts_chunk$set(cache = FALSE, fig.pos = "center", fig.width = 7, fig.height = 7)
# ```
#
# ## Setup
#
# ```{r packages, message=FALSE}
library("dplyr")
library("tidyr")
library("edgeR")
library("DT")
library("gplots")
library("ggplot2")
library("cowplot")
```

Input data.

# ```{r input-counts}
counts <- read.delim("../data/counts-filtered.txt", check.names = FALSE,
                     row.names = 1)
info <- read.delim("../data/experiment-info-filtered.txt",
                   stringsAsFactors = FALSE, row.names = 1)
stopifnot(colnames(counts) == rownames(info))
```

Calculate counts per million.

# ```{r calc-cpm}
counts_cpm <- cpm(counts, log = TRUE)

##min:  Count per million  is dat_cpm
counts_cpm <- dat_cpm
## PCA

# ```{r calc-pca}
pca <- prcomp(t(counts_cpm), scale. = TRUE)
variances <- pca$sdev^2
explained <- variances / sum(variances)
plot(pca, main = "Variance per PC")


Calculate the relationship between each recorded covariate and the top 6 PCs.

```{r pc-covariate-correlation}
p_comps <- 1:6
pc_cov_cor <- matrix(nrow = ncol(anno), ncol = length(p_comps),
                     dimnames = list(colnames(anno), colnames(pca$x)[p_comps]))
for (pc in p_comps) {
  for (covariate in 1:ncol(anno)) {
    lm_result <- lm(pca$x[, pc] ~ anno[, covariate])
    r2 <- summary(lm_result)$r.squared
    pc_cov_cor[covariate, pc] <- r2
  }
}
datatable(pc_cov_cor)
```

# ```{r pc-covariate-correlation-heatmap, fig.width=8}
heatmap.2(pc_cov_cor, trace = "none", margins = c(4, 8))


```{r long}
pc_cov_cor_2 <- as.data.frame(pc_cov_cor)
pc_cov_cor_2$variable <- rownames(pc_cov_cor)
pc_cov_cor_2 <- gather(pc_cov_cor_2, key = "pc", value = "cor", -variable)
head(pc_cov_cor_2)
```

```{r pca-heatmap}
d_heatmap <- pc_cov_cor_2
d_heatmap$variable <- factor(d_heatmap$variable,
                             levels = c("rin", "master_mix",
                                        "extraction", "arrival", "infection",
                                        "individual", "status", "treatment"),
                             labels = c("RNA quality", "Library prep batch",
                                        "RNA extraction batch", "Arrival batch",
                                        "Infection batch", "Individual",
                                        "Susceptibility status", "Treatment"))
pca_heat <- ggplot(d_heatmap, aes(x = pc, y = variable)) +
  geom_tile(aes(fill = cor), colour = "white") +
  scale_fill_gradient(low = "white", high = "red", limits = c(0, 1)) +
  labs(x = "Principal Component", y = "",
       title = "Correlation between principal components and experimental variables")
pca_heat
```


```{r pca-data}
pca_data <- cbind(info, pca$x[, p_comps])
pca_data$status <- factor(pca_data$status, levels = c("contact", "tb"),
                          labels = c("resistant", "susceptible"))
pca_data$treatment <- factor(pca_data$treatment, levels = c("none", "infected"))
```

PC1 versus PC2.

```{r pc1-pc2}
pc1v2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = treatment)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 1, round(explained[1] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 2, round(explained[2] * 100, 2)))
pc1v2
```

PC3 versus PC4.

```{r pc3-pc4}
pc3v4 <- ggplot(pca_data, aes(x = PC3, y = PC4, color = infection)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 3, round(explained[3] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 4, round(explained[4] * 100, 2)))
pc3v4
```

PC5 versus PC6.


pc5v6 <- ggplot(pca_data, aes(x = PC5, y = PC6, color = infection)) +
  geom_text(aes(label = individual)) +
  labs(x = sprintf("PC%d (%.2f%%)", 5, round(explained[5] * 100, 2)),
       y = sprintf("PC%d (%.2f%%)", 6, round(explained[6] * 100, 2)))
pc5v6


## Effect of treatment

PC1 is correlated with treatment.

```{r}
pc_cov_cor_2 %>% filter(variable == "treatment", cor > 0.5)
pc1_treatment <- ggplot(pca_data, aes(x = treatment, y = PC1)) +
  geom_boxplot() +
  labs(title = "PC1 captures the effect of treatment")
pc1_treatment
```

## Any technical batch effects?

There are multiple factors with a strong correlation with at least one of the first 6 PCs.


pc_cov_cor_2 %>% filter(cor > 0.5)

Most of the major correlations are with biological variables: treatment or individual.
It is expected that the largest source of variation in gene expression is due to the transcriptional changes induced by the immune response to MTB.
Also, we expect there to be differences between individuals, and this will be explicitly modeled for in the differential expression analysis.
Thus the main concern is the variable infection.
The infection experiments were done as the whole blood samples were obtained, so they are the most confounded with susceptibility status.

```{r}
pc_cov_cor_2 %>% filter(variable == "infection", cor > 0.5)
table(info$infection[info$treatment == "none"],
      info$status[info$treatment == "none"])
```

The infection dates with susceptible samples are 2013-04-05, 2014-06-10, and 2014-11-24.
We clearly see there is variation in PC3 by the date of the infection experiment.

```{r pc3-infection}
pc3_infection <- ggplot(pca_data, aes(x = infection, y = PC3, fill = treatment)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "PC3 varies by date of infection")
pc3_infection
```

Reassuringly, however, this variation does not correlate strongly with susceptibilty status.

```{r pc3-status}
pc_cov_cor_2 %>% filter(variable == "status", pc == "PC3")
pc3_status <- ggplot(pca_data, aes(x = reorder(individual, PC3), y = PC3, color = status)) +
  geom_point() +
  facet_wrap(~treatment) +
  labs(title = "PC3 does not correlate with susceptibility status",
       x = "Individual") +
  theme(axis.text.x = element_text(angle = 90))
pc3_status
```

The situation is the same with PC5.
It varies with infection date.

```{r pc5-infection}
pc5_infection <- pc3_infection %+% aes(y = PC5) +
  labs(title = "PC5 varies by date of infection")
pc5_infection
```

But does not correlate strongly with susceptibilty status.

```{r pc5-status}
pc_cov_cor_2 %>% filter(variable == "status", pc == "PC5")
pc5_status <- pc3_status %+% aes(x = reorder(individual, PC5), y = PC5) +
  labs(title = "PC5 does not correlate with susceptibility status",
       x = "Individual")
pc5_status
```

## Multipanel figures

PCA

#```{r pca-mulit, fig.width=14, fig.height=14}
plot_grid(pc1v2 + theme(legend.position = "none"),
          pc3v4 + theme(legend.position = "none"),
          pc5v6 + theme(legend.position = "none"),
          pca_heat,
          labels = LETTERS[1:4])


Batch effects

#```{r effect-of-infection, fig.width=14, fig.height=14}
plot_grid(pc3_infection + theme(legend.position = "bottom"),
          pc5_infection  + theme(legend.position = "bottom"),
          pc3_status + theme(legend.position = "bottom"),
          pc5_status + theme(legend.position = "bottom"),
          labels = LETTERS[1:4])


## Session information

```{r info}
sessionInfo()
```
