library(swne)
library(umap)
library(NNLM)
library(dplyr)
library(Rtsne)

## Read counts
#counts <- read.table("/media/Scratch_SSD/cjwei/20190409_HCT116/Methyl/4_methRate/SC-filtered_only.methRate.csv",
#                     sep = ",", header = T, row.names = 1, check.names = F)
counts <- read.table("/media/Scratch_SSD/cjwei/20190409_HCT116/Methyl/4_methRate/cellLine-only.methRate.csv",
                     sep = ",", header = T, row.names = 1, check.names = F)
counts <- as.matrix(counts)

## Filter out features with too many NAs
na.frac <- apply(counts, 1, function(x) sum(is.na(x))/length(x))
counts <- counts[na.frac < 0.9,]
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

## Run tSNE on cell scores for visualization
tsne <- Rtsne(t(cell.scores), perplexity = 5) #Choose perplexity according to <https://distill.pub/2016/misread-tsne/>
row.names(tsne$Y) <- colnames(cell.scores)

## Pull out some labels (you can change this depending on what you want to extract)
#labels <- sapply(colnames(cell.scores), ExtractField, field = 1:2, delim = "_")
labels <- sapply(colnames(cell.scores), ExtractField, field = 1, delim = "_")


## Plot results
#pdf(file = "/media/Scratch_SSD/cjwei/20190409_HCT116/Methyl/4_methRate/SC-filtered_only.methRate.pdf")
pdf(file = "/media/Scratch_SSD/cjwei/20190409_HCT116/Methyl/4_methRate/cellLine-only.methRate.pdf")

PlotDims(umap.emb, sample.groups = labels, x.lab = "umap1", y.lab = "umap2", main.title = "HCT116 Single Cell Methylation (UMAP)", pt.size = 3, do.label = T, show.legend = F)
PlotDims(tsne$Y, sample.groups = labels, x.lab = "tsne1", y.lab = "tsne2", main.title = "HCT116 Single Cell Methylation (tSNE)", pt.size = 3, do.label = T, show.legend = F)

dev.off()
