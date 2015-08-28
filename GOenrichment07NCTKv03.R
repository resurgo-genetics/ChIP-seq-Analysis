## GOenrichment07NCTKv03.R
## Version 7.3 Thursday, December 10th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program which performs a Gene Ontology Enrichment Analysis.
# It takes the annotatePeaks output file or any .txt list of gene symbols with the header "Gene_Symbol" as input. 

## Run Speed:

## DIRECTIONS: 
# Change the appropriate parameters 
# Run the program.
# Output is four .txt files in the outputDirectory. One of the combined GO categories, and one for each
# individual GO category, BP (Biological Processes), MF (Molecular Function), CC (Cellular Component).

## Example Input:
# inputFile <- "/home/user/inputData/GATA1_5000up_2000down_allData.txt"
# outputDirectory <- "/home/user/outputData/GOenrichment" 
# annotationDB <- "org.Hs.eg.db" 
# runName <- "Test1"
# pValueCutoff <- 0.05 

## See this example when using different species than the human default.
# org.Hs.eg.db is a specific annotation package for annotating human genes only.
# A list of all annotation packages for every species can be found at the bioconductor website: http://www.bioconductor.org/packages/3.0/data/annotation/
# You will have to install the correct organism annotation Db. Example is for Mouse (mus musculus).
# Install by running in R:
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
# Then be sure to change the annotationDB parameter to "org.Mm.eg.db", example:
# annotationDB <- "org.Mm.eg.db"

##########################   SET PARAMETERS HERE   ####################################

## The input file and output directory must be complete file paths with no trailing /
## The input file is intended to be the output of the annotatePeaks Rscript/workflow

inputFile <- "/home/tylerk/massInputData/GATA1_5000up_2000down_allData.txt" #An annotate peaks .txt output file or a .txt list of gene symbols with the header "Gene_Symbol"
outputDirectory <- "/home/tylerk/massTestResults/GOenrichment" # Directory where output files will be located
annotationDB <- "org.Hs.eg.db" # Database used to convert gene symbols to entrez ids
runName <- "massTest1" # runName will be used as a prefix in the output file names
pValueCutoff <- 0.05 # Highest p value from the GO categories to include in the output file

##############################  END PARAMETER BLOCK ##################################

## Load required libraries
library(annotationDB, character.only=TRUE)
library(GOstats)
assign("annotdb", get(annotationDB))

symbolList <- read.delim(inputFile) # read in input file 
symbolList <- unique(as.character(symbolList$Gene_Symbol))

#symbolList <- as.character(symbolList$Gene_Symbol) # get list of gene symbols from the input file

## convert the input gene symbols into entrez IDs
geneList <- select(annotdb, keys=symbolList, columns="ENTREZID", keytype="SYMBOL") 
geneUniverse <- keys(annotdb, keytype="ENTREZID")
GOstatsGenes <- as.character(geneList$ENTREZID)

## perform GO enrichment analysis on BP terms
paramsBP <- new("GOHyperGParams", geneIds=GOstatsGenes, universeGeneIds=geneUniverse,
annotation=annotationDB, ontology="BP", pvalueCutoff=pValueCutoff, conditional=TRUE,
testDirection="over")

## perform GO enrichment analysis on CC terms
paramsCC <- new("GOHyperGParams", geneIds=GOstatsGenes, universeGeneIds=geneUniverse,
annotation=annotationDB, ontology="CC", pvalueCutoff=pValueCutoff, conditional=TRUE,
testDirection="over")

## perform GO enrichment analysis on MF terms
paramsMF <- new("GOHyperGParams", geneIds=GOstatsGenes, universeGeneIds=geneUniverse,
annotation=annotationDB, ontology="MF", pvalueCutoff=pValueCutoff, conditional=TRUE,
testDirection="over")

## format GO enrichment analysis for output
BP_GOdata <- hyperGTest(paramsBP)
BP_Summ <- summary(BP_GOdata)
BP_Summ <- data.frame(BP_Summ[1], BP_Summ[7], BP_Summ[2:6])
CC_GOdata <- hyperGTest(paramsCC)
CC_Summ <- summary(CC_GOdata)
CC_Summ <- data.frame(CC_Summ[1], CC_Summ[7], CC_Summ[2:6])
MF_GOdata <- hyperGTest(paramsMF)
MF_Summ <- summary(MF_GOdata)
MF_Summ <- data.frame(MF_Summ[1], MF_Summ[7], MF_Summ[2:6])

## write GO analysis outputs to file
outputFile <- paste(outputDirectory, "/", runName, sep="")

write.table(BP_Summ, file=paste(outputFile,"_BP.txt",sep=""), quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=TRUE)
write.table(CC_Summ, file=paste(outputFile,"_CC.txt",sep=""), quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=TRUE)
write.table(MF_Summ, file=paste(outputFile,"_MF.txt",sep=""), quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=TRUE)

## combine all GO analysis outputs into one file sorted by p values
colnames(MF_Summ)[1] <- "GOID"
colnames(BP_Summ)[1] <- "GOID"
colnames(CC_Summ)[1] <- "GOID"
temp <- merge(MF_Summ, BP_Summ, all=T)
ALL_Summ <- merge(temp, CC_Summ, all=T)
ALL_Summ <- ALL_Summ[with(ALL_Summ,order(Pvalue)),]
## write data to file.
write.table(ALL_Summ, file=paste(outputFile,"_ALL.txt",sep=""), quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=TRUE)
