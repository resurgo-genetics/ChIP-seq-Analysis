## installationScript03.R
## Version 4.0 Wednesday, December 10th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION
# This program installs all of the R packages necessary for the Kepler workflows or the other R Scripts. 
# The default TxDb and org.eg.db packages installed are those for hg19, mm9, and mm10. 
# See: http://bioconductor.org/packages/3.0/data/annotation/ if using other species or assemblys. 

# Run Speed: Usually less than 10 minutes.

## DIRECTIONS:
# Run the program.

## Example use of species or assembly other than hg19, mm9, mm10:

# For hg18 installation, run:
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Hsapiens.UCSC.hg18.knownGene")
# You will not need to install a different org db, since org.Hs.eg is installed by default and is for all assemblys of human.

# For C elegans, run:
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Celegans.UCSC.ce6.ensGene")
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Ce.eg.db")


########################### BEGIN CODING BLOCK ###########################
source("http://bioconductor.org/biocLite.R")
biocLite()
packages <- c("GenomicRanges", "rtracklayer", "GOstats", "GenomicAlignments", "gage", "pathview", "org.Hs.eg.db", "TxDb.Hsapiens.UCSC.hg19.knownGene", 
"TxDb.Mmusculus.UCSC.mm10.knownGene", "TxDb.Mmusculus.UCSC.mm9.knownGene", "org.Mm.eg.db", "Rsamtools")
biocLite(packages)
########################### END CODING BLOCK ############################