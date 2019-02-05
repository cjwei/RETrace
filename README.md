# RETrace
This repo contains multiple tools to perform retrospective lineage tracing and call cell type information from the same single cell using microsatellite and methylation, respectively.  Below is a summary of the analysis tools developed:
## Microsatellite Analysis
- **msCount**: This directory contains the scripts necessary to call microsatellite counts from each single cell.  This uses a custom-built script for microsatellite calling
- **HipSTR**: Contains scripts to analyze HipSTR output of msCounts/allelotype.

## Methylation Analysis
- **methylStats.py**: This script will calculate basic CpG coverage statistics from mapped RRBS reads
- **methylTrace.py**: Expanding upon the methylation calling from methylStats.py, we want to ca
lculate the pairwise dissimilarity between all files inputted.

#### Step 1: Convert all bam files to pileup (if pileup not detected)
This is done using the following command:
```r
samtools mpileup -f reference.fa input.bam >output.pileup
```

#### Step 2: Perform methylation calling using the pileup file
Contained within the `parsePileup` function, we can make methylation calls for each relevant C in CpG context

#### Step 3: Calculate the pairwise dissimilarity matrix
We will calculate the pairwise dissimilarity matrix as described in Hui et al, 2018.  For each pair of single cell bam files, we calculate the following:
1. Determine shared CpGs between each pair of cells
2. Filter out only CpGs with 0% or 100% methylation (and above the read cov cutoff [default: 1])
3. In order to calculate dissimilarity, we will assign CpG sites with the same methylation call to 0 and different to 100 then take the average over all shared sites
4. Output pairwise dissimilarity into csv file `output.pd.csv`

#### Step 4 (Optional): Calculate CpG capture statistics
This will output the following statistics into the `output.methStats.txt` file:
- Total CpG/CHG/CHH
- Unique CpG covered by >=1 read + Mean Cov
- Unique CpG covered by >=5 reads + Mean Cov
- Unique CpG covered by >=10 reads + Mean Cov
- Number CGI covered + Mean Cov
