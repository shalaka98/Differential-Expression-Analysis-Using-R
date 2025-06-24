<h1 align="center"> Visualization, Differential Expression Analysis and Downstream Analysis of Glioma Transcriptomics Count Data </h1>

The whole analysis will use a [count dataset of glioblastoma transcriptomic samples](https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv).
# Installing gplots for visualization
```{r, echo=TRUE, message=FALSE, warning=FALSE, results='hide'}
install.packages("tidyverse")
library(tidyverse)
```
# Generating Matrix from a CSV file
```{r}
gene_data <- read.csv('https://raw.githubusercontent.com/HackBio-Internship/public_datasets/main/Cancer2024/glioblastoma.csv', row.names=1)
mat <- as.matrix(gene_data)

log_data_matrix <- log2(mat+1) #log transformation of matrix as data range is too broad
head(log_data_matrix)
```

# Preparation for Differential Expression Analysis

# Calculate fold change
```{r, echo=TRUE, message=FALSE, warning=FALSE, results='hide'}
# group samples using index positions: separate treatment and control data
group1 <-  gene_data[, 1:5]
group2 <-  gene_data[, 6:10]

# Calculate the mean for control and treatment groups
group1_mean <- rowMeans(group1)
group2_mean <- rowMeans(group2)

# Calculate log2 fold change
log2_fold_change = log2(group1_mean/group2_mean)
log2_fold_change
```
# Calculate the pvalues of fold changes

```{r}
# Calculate the pvalues of fold changes
pvalues <- apply(gene_data, 1, function(row) {
  t.test(row[1:5], row[6:10])$p.value
})
pvalues

padj <- p.adjust(pvalues, method = "BH")

```

# Filtering out Significant Genes (padj <0.05 and fold change cutoff 2)
```{r}

gene_data$padj <- padj
gene_data$log2FC <- log2_fold_change

up_genes <- gene_data %>% filter(padj < 0.05, log2FC > 1)
down_genes <- gene_data %>% filter(padj < 0.05, log2FC < -1)


write.csv(gene_data, "background_genes.csv", row.names = TRUE)
write.csv(up_genes, "DE_significantly_up_genes.csv", row.names = TRUE)
write.csv(down_genes, "DE_significantly_down_genes.csv", row.names = TRUE)


```

### Using the default parameters in [ShinyGO 0.80](http://bioinformatics.sdstate.edu/go/), for 10 significantly differential expressed genes, we found the following top 5 enriched biological pathways, visualizated using R
Top 3 Enriched Pathways includes: Glutathione Derivative Metabolic process, Glutathione Derivative Metabolic Biosynthesis Process, the Linoleic Acid Metabolic Process

# Visualization of Biological Process Pathways
Dendrograms cluster genes/samples based on similarity.
Colors indicate high/low expression.
Creates a heatmap with row scaling.
```{r, fig.width=10, fig.height=10}
col_palette <- colorRampPalette(brewer.pal(11, "RdBu"))(100)
heatmap.2(x=mat, col = col_palette, 
          density.info = 'none', dendrogram = 'both',
          scale = 'row', trace = 'none')
```



