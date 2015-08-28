## gagePathway12NCTKv01.R
## Version 12.1, Monday, December 8th, 2014.
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
## This is a program which takes the output file from annotatePeaks as input,
## and calculates which pathways are enriched in the gene set 
## Output varies depending on how many pathways are found within the pValueCutoff.
## It includes a _full_pathway_list.txt file of all the pathways found, their pathway code, name, p and q values,
## .png image files and a matching .xml file for all pathways that meet the pValuecutoff for visalization of up and down regulation of pathway components. 


## Run Speed: 

## USE OF ALTERNATIVE INPUT FILES
## 1) This program can also function with an inputFile that has just two columns: entrezID and score
## 2) File must have entrezID and score as the column names
## 3) File must be tab-delimited and have unix line endings
## 4) EntrezID is the entrez id number for the gene; score should be a non-log2 value indicating enrichment or read height or other measure

## Example Input:
# inputFile <- "/home/user/inputData/GATA1_5000up_2000down_allData.txt"
# outputDirectory <- "/home/user/outputData/gagePathway/"
# runName <- "Test1" 
# pValueCutoff <- 0.2
# KEGGspeciesCode <- "hsa"

##########################   SET PARAMETERS HERE   #####################################

## The input file and output directory must be complete file paths, with no trailing /
inputFile <- "/home/tylerk/massInputData/GATA1_5000up_2000down_allData.txt" #Usually a .txt file from the annotatePeaks output, see above for alternate possible input files
outputDirectory <- "/home/tylerk/massTestResults/gagePathway" #Directory where output files will be located
runName <- "Test2" #runName that will be used as a prefix in all output file names
pValueCutoff <- 0.2 #Numeric highest p value from the pathway analysis to include in the output file
KEGGspeciesCode <- "hsa" #Species code for the KEGG database, "hsa" is human, "mmu" is mouse, See http://www.genome.jp/kegg/catalog/org_list.html for complete list.

#############################  END OF PARAMETER BLOCK ###################################

library(gage)
library(pathview)
#outputDirectory <- paste(outputDirectory, "/", runName, "_Pathway_data", sep="")
#system(paste("mkdir", outputDirectory, sep=" "))
setwd(outputDirectory)

## gets the list of entrez genes and their associated peak scores from the input file
geneData <- unique(read.delim(inputFile)[,c("entrezID", "score")])
geneList <- unique(geneData$entrezID)
geneList <- geneList[!is.na(geneList)]
dataMatrix <- matrix(data=NA, nrow=length(geneList), ncol=1) 
rownames(dataMatrix) <- geneList
colnames(dataMatrix) <- "score"

## for each unique entrez gene, finds the highest peak score associated with that gene
for (i in 1:length(geneList)) {
	gene <- as.character(geneList[i])
	#print(gene)
	score <- max(geneData$score[which(geneData$entrezID == gene)])
	dataMatrix[gene,] <- score
}


data(kegg.gs) # load in the KEGG pathway data
gageOutput <- gage(dataMatrix, gsets = kegg.gs, samp=1) # determines which pathways are enriched in the data set
pathwayData <- gageOutput$greater[,c("p.val", "q.val")]
pathCode <- substr(rownames(pathwayData), 1, 8)
pathName <- substr(rownames(pathwayData), 9, length(rownames(pathwayData)))
pathwayData <- cbind(pathCode, pathName ,pathwayData)
colnames(pathwayData) <- c("Pathway Code", "Pathway Name", "p Values", "q value")
write.table(pathwayData, paste(runName, "full_pathway_list.txt", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

## get enriched pathways with a p value less than the cutoff value
sel <- gageOutput$greater[,"p.val"] < pValueCutoff & !is.na(gageOutput$greater[,"p.val"])
pathIds <- substr(rownames(gageOutput$greater)[sel], 1, 8)

## Use pathview to create output images
runPathview <- function(pid){
	print(pid)
	pathview(gene.data=dataMatrix, pathway.id=pid, kegg.native=TRUE, species=KEGGspeciesCode, plot.col.key = FALSE)
	system(paste("rename 's/^/", runName, "_/' ", pid, "*", sep=""))
	}
pathList <- sapply(pathIds, function(pid) runPathview(pid))

