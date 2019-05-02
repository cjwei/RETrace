library(swne)
library(umap)
library(NNLM)
library(dplyr)

## Read counts
counts <- read.table("/media/Scratch_SSD/cjwei/20190409_Methyl-Brain/5_plotPCA/20190409_Brain_SC-only.methRate.csv",
                     sep = ",", header = T, row.names = 1)
counts <- as.matrix(counts)

## Filter out features with too many NAs
na.frac <- apply(counts, 1, function(x) sum(is.na(x))/length(x))
counts <- counts[na.frac < 0.8,]
dim(counts)

## Run NMF
nmf.res <- NNLM::nnmf(t(counts), k = 12, n.threads = 8, max.iter = 500, check.k = F)
cpg.loadings <- nmf.res$H
cell.scores <- t(nmf.res$W)

# ## Mean imputation for any remaining missing values
# cell.scores <- t(apply(cell.scores, 1, function(x) {
#    x[is.na(x)] <- mean(na.omit(x))
#    x
#  }))

## Heatmap of cell scores
ggHeat(cell.scores, clustering = "both", labRow = T, labCol = F)

## Run UMAP on cell scores for visualization
umap.res <- umap(t(cell.scores))
umap.emb <- umap.res$layout

## Pull out some labels (you can change this depending on what you want to extract)
labels <- sapply(colnames(cell.scores), ExtractField, field = 2, delim = "_")
#labels <- sapply(colnames(cell.scores), ExtractField, field = 2, delim = "\\.")
#labels <- sapply(labels, ExtractField, field = 1, delim = "_")

## Plot results
PlotDims(umap.emb, sample.groups = labels, x.lab = "umap1", y.lab = "umap2")
