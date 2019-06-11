library(swne)
library(umap)
library(NNLM)
library(dplyr)
library(Rtsne)

## Read counts
counts <- read.table("/media/Scratch_SSD/cjwei/20190409_Methyl-HCT116/2_RETrace-methDict/HCT116_SC-only.SC-filtered.200targets.300000CpG.methRate.csv",
                     sep = ",", header = T, row.names = 1)
counts <- as.matrix(counts)

## Filter out features with too many NAs
na.frac <- apply(counts, 1, function(x) sum(is.na(x))/length(x))
counts <- counts[na.frac < 0.9,]
#counts.cortex <- counts[,grepl("Cortex", colnames(counts))] #This contains CpG methRate for cortex samples only
#counts.NeuNpos <- counts[,grepl("Cortex1.NeuNpos|Cortex2.NeuNpos", colnames((counts)))] #This contains CpG meth rate for cortex NeuNpos samples only
dim(counts)

## Run NMF
nmf.res <- NNLM::nnmf(t(counts), k = 12, n.threads = 8, max.iter = 500, check.k = F)
cpg.loadings <- nmf.res$H
cell.scores <- t(nmf.res$W)

#nmf.res.cortex <- NNLM::nnmf(t(counts.cortex), k = 12, n.threads = 8, max.iter = 500, check.k = F)
#cpg.loadings.cortex <- nmf.res.cortex$H
#cell.scores.cortex <- t(nmf.res.cortex$W)

#nmf.res.NeuNpos <- NNLM::nnmf(t(counts.NeuNpos), k = 12, n.threads = 8, max.iter = 500, check.k = F)
#cpg.loadings.NeuNpos <- nmf.res.NeuNpos$H
#cell.scores.NeuNpos <- t(nmf.res.NeuNpos$W)

# ## Mean imputation for any remaining missing values
# cell.scores <- t(apply(cell.scores, 1, function(x) {
#    x[is.na(x)] <- mean(na.omit(x))
#    x
#  }))

## Heatmap of cell scores
ggHeat(cell.scores, clustering = "both", labRow = T, labCol = F)
#ggHeat(cell.scores.cortex, clustering = "both", labRow = T, labCol = F)
#ggHeat(cell.scores.NeuNpos, clustering = "both", labRow = T, labCol = F)

## Run UMAP on cell scores for visualization
umap.res <- umap(t(cell.scores))
umap.emb <- umap.res$layout

#umap.res.cortex <- umap(t(cell.scores.cortex))
#umap.emb.cortex <- umap.res.cortex$layout

#umap.res.NeuNpos <- umap(t(cell.scores.NeuNpos))
#umap.emb.NeuNpos <- umap.res.NeuNpos$layout

## Run tSNE on cell scores for visualization
tsne <- Rtsne(t(cell.scores), perplexity = 5) #Choose perplexity according to <https://distill.pub/2016/misread-tsne/>
row.names(tsne$Y) <- colnames(cell.scores)

#tsne.cortex <- Rtsne(t(cell.scores.cortex), perplexity = 5)
#row.names(tsne.cortex$Y) <- colnames(cell.scores.cortex)

#tsne.NeuNpos <- Rtsne(t(cell.scores.NeuNpos), perplexity = 5)
#row.names(tsne.NeuNpos$Y) <- colnames(cell.scores.NeuNpos)

## Pull out some labels (you can change this depending on what you want to extract)
labels <- sapply(colnames(cell.scores), ExtractField, field = 2, delim = "_")
#labels.cortex <- sapply(colnames(cell.scores.cortex), ExtractField, field = 2, delim = "_")
#labels.NeuNpos <- sapply(colnames(cell.scores.NeuNpos), ExtractField, field = 2, delim = "_")
#labels <- sapply(colnames(cell.scores), ExtractField, field = 2, delim = "\\.")
#labels <- sapply(labels, ExtractField, field = 1, delim = "_")

## Plot results
pdf(file = "/media/Scratch_SSD/cjwei/20190409_Methyl-HCT116/2_RETrace-methDict/20190409_HCT116_SC-only.min_CpG1.methRate.pdf")

PlotDims(umap.emb, sample.groups = labels, x.lab = "umap1", y.lab = "umap2", main.title = "HCT116 Single Cell Methylation (UMAP)")
PlotDims(tsne$Y, sample.groups = labels, x.lab = "tsne1", y.lab = "tsne2", main.title = "HCT116 Single Cell Methylation (tSNE)")

#PlotDims(umap.emb.cortex, sample.groups = labels.cortex, x.lab = "umap1", y.lab = "umap2", main.title = "Cortex only (UMAP)")
#PlotDims(tsne.cortex$Y, sample.groups = labels.cortex, x.lab = "tsne1", y.lab = "tsne2", main.title = "Cortex only (tSNE)")

#PlotDims(umap.emb.NeuNpos, sample.groups = labels.NeuNpos, x.lab = "umap1", y.lab = "umap2", main.title = "NeuNpos only (UMAP)")
#PlotDims(tsne.NeuNpos$Y, sample.groups = labels.NeuNpos, x.lab = "tsne1", y.lab = "tsne2", main.title = "NeuNpos only (tSNE)")
#We want to label all SC uniquely to get a clearer view of what is happening
#labels.NeuNpos.SC <- sapply(colnames(cell.scores.NeuNpos), ExtractField, field = 4, delim = "\\.")
#PlotDims(umap.emb.NeuNpos, sample.groups = labels.NeuNpos.SC, x.lab = "umap1", y.lab = "umap2", show.legend = F, main.title = "NeuNpos only (UMAP); All SC labeled")
#PlotDims(tsne.NeuNpos$Y, sample.groups = labels.NeuNpos.SC, x.lab = "tsne1", y.lab = "tsne2", show.legend = F, main.title = "NeuNpos only (tSNE); All SC labeled")

dev.off()