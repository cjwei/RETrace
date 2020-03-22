# RETrace: simultaneous retrospective lineage tracing and methylation profiling of single cells
This repo contains multiple tools to perform retrospective lineage tracing and call cell type information from the same single cell using microsatellite and methylation, respectively.  Below is a summary of the analysis tools developed.  Paper can be found on Genome Research, doi:10.1101/gr.255851.119
## Microsatellite Analysis
- **HipSTR_allelotype**: Run HipSTR for microsatellite calling and allelotyping
- **Custom_allelotype**: Run scripts for custom microsatellite calling and allelotyping (based off of LobSTR)
- **merge_allelotype**: Merge multiple alleleDict pickle files from either Custom or HipSTR allelotyping
- **buildPhylo**: Build phylogenetic tree given allelotype of single cells
- **iterPhylo**: Iteratively build phylogenetic tree given allelotype of single cells
- **evalPhylo**: Evaluate accuracy of phylogeny given exVivo tree data
- **viewPhylo**: Utilize ete3 to view phylogenetic tree derived from buildPhylo
## Methylation Analysis
- **importMethyl**: Import methylation calls derived from methylpy TSV files for samples/cell type ref (make sure to run methylpy prior to import)
- **refPD**: Calculate and plot pairwise dissimilarity between single cells and reference cell types
- **pairwise_methRate**: Calculate methRate across Ensembl Regulatory Build windows
- **combinedPhylo**: Build phylogenty by combining MS and Methyl distances between single cell pairs
