We can download the locations of CpG islands in hg19 genome from UCSC genome browser using the following commands: (from <https://www.biostars.org/p/236141/>)
wget -qO- http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cpgIslandExt.txt.gz \
   | gunzip -c \
   | awk 'BEGIN{ OFS="\t"; }{ print $2, $3, $4, $5$6, substr($0, index($0, $7)); }' \
   | sort-bed - \
   > cpgIslandExt.hg19.bed
