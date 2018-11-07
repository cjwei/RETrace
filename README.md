# RETrace
This repo contains multiple tools to perform retrospective lineage tracing and call cell type information from the same single cell using microsatellite and methylation, respectively.  Below is a summary of the analysis tools developed:
## Microsatellite Analysis
- **msCount**: This directory contains the scripts necessary to call microsatellite counts from each single cell.  This uses a custom-built script for microsatellite calling
- **HipSTR**: Contains scripts to analyze HipSTR output of msCounts/allelotype.
## Methylation Analysis
- **methylCov.py**: This script will calculate basic CpG coverage statistics from mapped RRBS reads
