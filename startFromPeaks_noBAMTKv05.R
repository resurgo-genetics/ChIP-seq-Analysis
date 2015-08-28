## startFromPeaks_noBAMTKv05.R
## Version 5.0 Wednesday, December 10th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program to change run the full pipeline script starting from Mapped Reads without the original mapped read file using just the peak file.

## Note: This program must be run from the command line and not using RStudio.

## Run Speed: Long, dependant on maxDremeRunTime parameter.

## DIRECTIONS:
# Change the appropriate parameters 
# Run the program.
# Output is one .jpeg file into the outputDirectory with the name 

## EXAMPLES:

## To Run:
# Type the file name in a terminal window: Rscript startFromPeaks_noBAMv05.R 

## Sample Input for Parameter Block:
# MACSpeakFile <- "/home/user/inputData/GATA1_peaks.bed"
# outputDirectory <- "/home/user/Results/startPeaksNoBAM"
# keplerBinLocation <- "/home/user/CWF1/keplerBin"
# genomeAssembly <- "hg19"
# transcriptionDB <- "TxDb.Hsapiens.UCSC.hg19.knownGene"
# annotationDB <- "org.Hs.eg.db"
# runName <- "Test1"
# extendReads <- "yes" 
# fragmentSize <- 500 
# useControl <- "no"
# controlFile <- "NA"
# broadCutOff <- 10000 
# upstreamRange <- 5000 
# downstreamRange <- 2000 
# GOpValueCutoff <- 0.05 
# pathway_pValueCutoff <- 0.2 
# KEGGspeciesCode <- "hsa"
# genomeDirectory <- "/home/user/keplerBin/genomes/hg19"
# DREMEeValueCutoff <- 0.001
# maxMotifsToGenerate <- 10
# maxDremeRuntime <- 60 
# titlePrefix <- "My Test Graph!" 
# useCustomYrange <- "yes"
# yMin <- 0
# yMax <- 45
# numberOfExamples <- 10 
# rangeAroundPeak <- 5000 
# nonPeakRegionSize <- 5000 

## See these examples when using different species or genome assembly, all annotation packages can be found at the bioconductor website: http://www.bioconductor.org/packages/3.0/data/annotation/.

## Example Use of Different TxDB. For Mmusculus assembly mm10:
# Visit Bioconductor website link above, find correct TxDB package in list, click link and 
# Install by running in R:
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
# Then be sure to change transcriptDB parameter to "TxDb.Mmusculus.UCSC.mm10.ensGene" and the genomeAssembly to mm10, example:
# transcriptionDB <- "TxDb.Mmusculus.UCSC.mm10.ensGene" 
# genomeAssembly <- "hg19"

## You will also have to install the correct organism annotation Db. 
# Install by running in R:
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Mm.eg.db")
# Then be sure to change the annotationDB parameter to "org.Mm.eg.db", example:
# annotationDB <- "org.Mm.eg.db"

## The genomeDirectory parameter for genomic data must also be changed. hg19 and mm18 are included by default but for any others they must be downloaded ftp://ftp.ncbi.nih.gov/genomes/ and added to the keplerBin.
# Change the genome directory to the new species, example for mm18:
# genomeDirectory <- "/home/user/CWF1/keplerBin/genomes/mm18"

## The KEGGspeciesCode parameter in pathway options must also be changed, "hsa"=human, "mmu"=mouse, See http://www.genome.jp/kegg/catalog/org_list.html for complete list.
# Change the KEGGspeciesCode parameter to the new species, example for mouse:
# KEGGspeciesCode <- "mmu"

## IMPORTANT HOMER SPECIES AND ASSEMBLY GENOME INFORMATION:
# The correct species and genome assembly must be installed within Homer for this program to work. 
# Note: This depends on Homer, Python and Perl being installed correctly. Homer must be in the path.

# To install this:

# Run this perl command to see potential lists of genomes, be sure to change /home/user/Homer to the correct path where you have installed Homer:
# perl /home/user/Homer//configureHomer.pl -list

# From this list get the correct genome code, hg19 = human genome assembly 19, mm10 = mus musculus 10. 
# Run this perl command to install the genome, example for mm10, be sure to change /home/user/Homer to the correct path where you have installed Homer:
# perl /home/user/Homer//configureHomer.pl -install mm10

# The genomeAssembly parameter within the parameter block  must correctly indicate the species & assembly. 
# Change the genome assembly to the new species, example for mm10:
# genomeAssembly <- "mm10"

##########################   SET PARAMETERS HERE   ####################################

## The MACSpeakFile, outputDirectory and keplerBinLocation must be complete file paths, with no trailing /
## Note: This workflow assumes you don't have the sequence read data, it uses a .bed file, NOT a .sam or .bam file.

MACSpeakFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" # Should be a .bed file of peak data.
outputDirectory <- "/home/tylerk/massTestResults/startPeaksNoBAM" # The directory where results will be outputted to
keplerBinLocation <- "/home/tylerk/CWF1/keplerBin" # The location of the keplerBin where the python scripts are located
genomeAssembly <- "hg19" # Assembly of genome being used, see above example for more info
transcriptionDB <- "TxDb.Hsapiens.UCSC.hg19.knownGene" # TxDb package being used, see above example for more info
annotationDB <- "org.Hs.eg.db" # The database used for annotation, see above example for more info

runName <- "massTest1" # runName will be used as a prefix in the output file names

extendReads <- "yes" # If set to yes, reads will be extended, otherwise they will be left unaltered
fragmentSize <- 500 # Numeric estimated fragment size after sonication

## MACS options
useControl <- "no" # Specify "yes" or "no" to use or not use a control
controlFile <- "NA" # Complete file path of control file, example: "/home/user/inputData/GATA1control_peaks.bed", no trailing /

## peakStats options
broadCutOff <- 10000 # Numeric number of bases a peak must have to be considered a broad peak

## annotate peaks options
upstreamRange <- 5000 # Numeric maximum distance upstream of TSS to look for a peak
downstreamRange <- 2000 # Numeric maximum distance downstream of TSS to look for a peak

## GO options
GOpValueCutoff <- 0.05 # Numeric highest p value from the GO categories to include in the output file

## Pathway options
pathway_pValueCutoff <- 0.2 # Numeric highest p value from the pathway analysis to include in the output file
KEGGspeciesCode <- "hsa" # Species code from KEGG database, "hsa"= human, "mmu" = mouse, more at: http://www.genome.jp/kegg/catalog/org_list.html

## get peak sequences options
genomeDirectory <- "/home/tylerk/CWF1/keplerBin/genomes/hg19" # The directory of your genome file, within the kepler bin. hg19=human genome 19, Change if using different species or assembly!

## DREME options
DREMEeValueCutoff <- 0.001 # Numeric e value cut off for motif generation
maxMotifsToGenerate <- 10 # Numeric maximum number of motifs to generate
maxDremeRuntime <- 240  # Numeric run time in seconds. This is the parameter that will affect total program run time the most. At least 1 hour is recommended (3600). Set to 0 for no limit. 

## peak example options
titlePrefix <- "My Test Graph!" # The title prefix will be added to the title of the graph creates


useCustomYrange <- "yes" # Set to "yes" to restrict y axis of the graph to custom settings. Set to "no" for automatic y axis to fit input data.
yMin <- 0 # The minimum range of the y axis if useCustomYrange is set to "yes"
yMax <- 45 # The maximum range of the y axis if useCustomYrange is set to "yes"

numberOfExamples <- 10 # Numeric number of peak example and non-peak example coverage plots to create
rangeAroundPeak <- 5000 # Numeric number of bases up and downstream of the peak to include in the plot
nonPeakRegionSize <- 5000 # Numeric total number of bases to include in each non-peak example coverage plot

#################################   END PARAMETER BLOCK  ########################################

########################################## peakStats ########################################## 

genomeSizeFile <- paste(outputDirectory, "/", runName, "_chromosomeSizes.txt", sep="")
shellCommand <- paste("sh ", keplerBinLocation, "/fetchChromSizes.sh", sep="")
system(paste(shellCommand, genomeAssembly, ">", genomeSizeFile, sep=" "))
outputFile <- paste(outputDirectory, "/", runName, "_peakStats.txt", sep="")
dataIn <- read.table(MACSpeakFile)

totalPeaks <- nrow(dataIn) # count number of peaks
dataIn$peakLength <- dataIn$V3 - dataIn$V2 +1 #TK added +1
smallestPeak <- min (dataIn$peakLength) # find size of the smallest peak
largestPeak <- max (dataIn$peakLength) # find size of the largest peak
meanLength <- mean(dataIn$peakLength) # find the mean peak length
SD <- sd(dataIn$peakLength) # find the standard deviation of the peak length
medianLength <- median(dataIn$peakLength) # find the median peak length

## calculates the total genome size from the genomeSize input file if necessary
if (!is.na(suppressWarnings(as.numeric(genomeSizeFile)))) {
	totalBases <- as.numeric(genomeSizeFile)
} else {
	lengths <- read.table(genomeSizeFile)
	trimmedLengths <- lengths[!grepl(".*_.*", lengths$V1), ]
	totalBases <- sum(as.numeric(lengths$V2))
}
coveragePercent <- (sum(dataIn$peakLength)/totalBases) * 100 # calculates % of genome covered in peaks

## determines how many peaks are classified as broad and what percentage of the total peaks this is
broad <- dataIn[dataIn$peakLength > broadCutOff, ]
narrow <- dataIn[dataIn$peakLength <= broadCutOff, ]
broadPercent <- (nrow(broad)/nrow(dataIn)) * 100

## set up the data for the output file
totalPeaks_str <- paste("Total number of peaks:", totalPeaks, sep="\t")
smallestPeak_str <- paste("Min peak length (bp):", smallestPeak, sep="\t")
largestPeak_str <- paste("Max peak length (bp):", largestPeak, sep="\t")
meanLength_str <- paste("Mean peak length (bp):", round(meanLength, digits=3), sep="\t")
SD_str <- paste("Peak length standard deviation:", round(SD, digits=3), sep="\t")
medianLength_str <- paste("Median peak length (bp):", medianLength, sep="\t")
broad_str <- paste("number of peaks longer than ", broadCutOff/1000, " kb:\t", nrow(broad), sep="")
broadPercent_str <- paste("% of peaks larger than ", broadCutOff/1000, " kb:\t", round(broadPercent, digits=3), sep="")
coveragePercent_str <- paste("% of genome with peaks:", round(coveragePercent, digits=3), sep="\t")

## calculates the number of peaks on each chromosome and sets up output string
chrDistro_str <- paste("Number of peaks per chromosome:", "\n")
for (i in 1:22){
	chrNum <- paste("chr", i, sep="")
	temp <- nrow(dataIn[(dataIn$V1 == chrNum), ])
	chrDistro_str <- paste(chrDistro_str, chrNum, temp, "\n", sep="\t")
}

## adds the correct chromosome symbol to the non-numeric chromosomes	
chrNum <- "chrX"
temp <- nrow(dataIn[(dataIn$V1 == chrNum), ])
chrDistro_str <- paste(chrDistro_str,chrNum, temp, "\n", sep="\t")
chrNum <- "chrY"
temp <- nrow(dataIn[(dataIn$V1 == chrNum), ])
chrDistro_str <- paste(chrDistro_str,chrNum, temp, "\n", sep="\t")
chrNum <- "chrM"
temp <- nrow(dataIn[(dataIn$V1 == chrNum), ])
chrDistro_str <- paste(chrDistro_str,chrNum, temp, "\n", sep="\t")

## writes the output string to file
stringOutput <- paste(totalPeaks_str, smallestPeak_str, largestPeak_str, medianLength_str, meanLength_str, SD_str, coveragePercent_str, broad_str, broadPercent_str, chrDistro_str, sep="\n")

writeLines(stringOutput, outputFile)

########################################## distanceToTSS ########################################## 

## takes in a set of peak data in .bed format
## uses the specified TxDb to calculate distances between each peak and its nearest TSS
## see the bioconductor website for a list of available TxDbs

## Load required libraries
library(transcriptionDB, character.only=TRUE)
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)
library(GenomicAlignments)
exampleDirectory <- paste(outputDirectory, "peak_examples", sep="/")
system(paste("mkdir ", exampleDirectory, sep=""))
assign("txdb", get(transcriptionDB))

## function pull the TSS from the correct stands of the TxDb data
grabTSS <- function(x){
	if(x[3]=="-"){
		x[length(x)+1] <- as.integer(x[5])
	}
	else{
		x[length(x)+1] <- as.integer(x[4])
	}
	#str(x)
	return(x)
}

## extract data from the TxDb
keys <- keys(txdb, keytype="TXID")
columns=c("TXSTART", "TXEND","TXSTRAND", "TXCHROM")
TSSdata <- select(txdb, keys=keys, columns=columns, keytype="TXID")
TSSdata <- data.frame(t(apply(TSSdata, 1, grabTSS)))
TSSdata$V6 <- as.numeric(levels(TSSdata$V6))[TSSdata$V6]

## create Granges object for the TSSs
TSSranges <- GRanges(seqnames = Rle(TSSdata$TXCHROM),
ranges = IRanges(as.numeric(TSSdata$V6), as.numeric(TSSdata$V6)),
strand = Rle(TSSdata$TXSTRAND))

peaks <- import.bed(con=MACSpeakFile, asRangedData=F) # read in peaks as a Grange object

distances <- as.data.frame(distanceToNearest(peaks, TSSranges)) # calculate distance between each peak and its nearest TSS

## determine number of peaks which fall into each distance interval
overlap <- nrow(distances[distances$distance==0,])
under1kb <- nrow(distances[distances$distance<=1000 & distances$distance>0,])
under2kb <- nrow(distances[distances$distance<=2000 & distances$distance>1000,])
under5kb <- nrow(distances[distances$distance<=5000 & distances$distance>2000,])
under10kb <- nrow(distances[distances$distance<=10000 & distances$distance>5000,])
under20kb <- nrow(distances[distances$distance<=20000 & distances$distance>10000,])
over20 <- nrow(distances[distances$distance>20000,])

## format data for output file
overlap <- paste("overlapping TSS", overlap, sep="\t")
within1Kb <- paste("within 1Kb of TSS", under1kb, sep="\t")
within2Kb <- paste("between 1Kb and 2Kb", under2kb, sep="\t")
within5Kb <- paste("between 2Kb and 5Kb", under5kb, sep="\t")
within10Kb <- paste("between 5Kb and 10Kb", under10kb, sep="\t")
within20Kb <- paste("between 20Kb and 10Kb", under20kb, sep="\t")
over20Kb <- paste("over 20Kb from TSS", over20, sep="\t")
outString <- "Number of peaks within set intervals of the TSS"
outString <- paste(outString, overlap, within1Kb, within2Kb, within5Kb, within10Kb, within20Kb, over20Kb, sep="\n")

## print output to file
outputFile <- paste(outputDirectory, "/", runName, "_distanceToTSSdata.txt", sep="")
writeLines(outString, outputFile)

########################################## annotatePeaks ########################################## 

library(transcriptionDB, character.only=TRUE)
assign("txdb", get(transcriptionDB))
library(annotationDB, character.only=TRUE)
assign("annotdb", get(annotationDB))
library(GenomicRanges)
library(rtracklayer)

peaks <- import.bed(con=MACSpeakFile, asRangedData=F)
start(peaks) <- start(peaks)-1

#function for extending the TSS to the ranges specified by upstreamRange and downstreamRange
extendTSS <- function(x){
	if(x[5]=="-"){
		x[length(x)+1] <- as.integer(x[7])
		x[length(x)+1] <- as.integer(x[7])-downstreamRange
		x[length(x)+1] <- as.integer(x[7])+upstreamRange
	}
	else{
		x[length(x)+1] <- as.integer(x[6])
		x[length(x)+1] <- as.integer(x[6])-upstreamRange
		x[length(x)+1] <- as.integer(x[6])+downstreamRange
	}
	return(x)
}

## get information from the TxDb and extend the TSSs
keys <- keys(txdb, keytype="TXID")
columns=c("TXSTART", "TXEND","TXSTRAND", "TXCHROM", "GENEID", "TXNAME")
TSSdata <- select(txdb, keys=keys, columns=columns, keytype="TXID")
TSSdata <- data.frame(t(apply(TSSdata, 1, extendTSS)))
TSSdata$V8 <- as.numeric(levels(TSSdata$V8))[TSSdata$V8]
TSSdata$V9 <- as.numeric(levels(TSSdata$V9))[TSSdata$V9]
TSSdata$V10 <- as.numeric(levels(TSSdata$V10))[TSSdata$V10]

## create Granges object for the TSS ranges
TSSranges <- GRanges(seqnames = Rle(TSSdata$TXCHROM),
ranges = IRanges(as.numeric(TSSdata$V9), as.numeric(TSSdata$V10)),
strand = Rle(TSSdata$TXSTRAND), entrezID=TSSdata$GENEID, UCSCID=TSSdata$TXNAME, TSS=TSSdata$V8)

## find overlaps between peak data and TSS ranges
peakData <- data.frame(as.character(seqnames(peaks)), ranges(peaks), peaks@elementMetadata) 
TSSdata <- data.frame(TSSranges@elementMetadata)
overlap <- findOverlaps(peaks,TSSranges)

## extract overlapping peaks and associated data
outputData <- data.frame(as.character(seqnames(peaks))[queryHits(overlap)], ranges(peaks)[queryHits(overlap)], 
peaks@elementMetadata[queryHits(overlap),],ranges(TSSranges)[subjectHits(overlap)], TSSranges@elementMetadata[subjectHits(overlap),], 
strand(TSSranges)[subjectHits(overlap)])
colnames(outputData)[1] <- "chromosome"
colnames(outputData)[13] <- "strand"
colnames(outputData)[which(colnames(outputData)=="start")] <- "peak start"
colnames(outputData)[which(colnames(outputData)=="end")] <- "peak end"
colnames(outputData)[which(colnames(outputData)=="name")] <- "peak name"
colnames(outputData)[which(colnames(outputData)=="UCSCID")] <- "UCSC_ID"
colnames(outputData)[which(colnames(outputData)=="TXSTART")] <- "TSS"
colnames(outputData)[which(colnames(outputData)=="GENEID")] <- "entrez ID"
colnames(outputData)[which(colnames(outputData)=="start.1")] <- "TSS_region_start"
colnames(outputData)[which(colnames(outputData)=="end.1")] <- "TSS_region_end"

## remove extra data fields
drops <- c("width", "TXCHROM", "width.1")
outputData <- outputData[,!(names(outputData) %in% drops)]

## get additional data from the annotationDB
keys <- keys(annotdb, keytype="UCSCKG")
columns=c("SYMBOL","GENENAME")
annotData <- select(annotdb, keys=keys, columns=columns, keytype="UCSCKG")

## add new annotation data to the previous data list
colnames(annotData)[which(colnames(annotData)=="UCSCKG")] <- "UCSC_ID"
colnames(annotData)[which(colnames(annotData)=="SYMBOL")] <- "Gene_Symbol"
colnames(annotData)[which(colnames(annotData)=="GENENAME")] <- "Full Gene Name"
outputData <- merge(outputData, annotData, by="UCSC_ID")
outputData <- outputData[c(2:6, 10, 7:8, 11, 9, 1, 12:13)]

## write data tables to file
symbolList <- unique(outputData$Gene_Symbol)
annotationFile <- paste(outputDirectory, "/", runName, "_", upstreamRange, "up_", downstreamRange, "down_allData.txt", sep="")
symbolFile <- paste(outputDirectory, "/", runName, "_", upstreamRange, "up_", downstreamRange, "down_geneSymbols.txt", sep="")
write.table(outputData, file=annotationFile,append=FALSE, quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=TRUE)
write.table(symbolList, file=symbolFile,append=FALSE, quote=FALSE, sep="\t", eol = "\n", row.name=FALSE, col.name=FALSE)

########################################## GO analysis ########################################## 

## Load required libraries
library(annotationDB, character.only=TRUE)
library(GOstats)
assign("annotdb", get(annotationDB))

symbolList <- read.delim(symbolFile) # read in input file 
symbolList <- as.character(symbolList[[1]])

#symbolList <- as.character(symbolList$Gene_Symbol) # get list of gene symbols from the input file

## convert the input gene symbols into entrez IDs
geneList <- select(annotdb, keys=symbolList, columns="ENTREZID", keytype="SYMBOL") 
geneUniverse <- keys(annotdb, keytype="ENTREZID")
GOstatsGenes <- as.character(geneList$ENTREZID)

## perform GO enrichment analysis on BP terms
paramsBP <- new("GOHyperGParams", geneIds=GOstatsGenes, universeGeneIds=geneUniverse,
annotation=annotationDB, ontology="BP", pvalueCutoff=GOpValueCutoff, conditional=TRUE,
testDirection="over")

## perform GO enrichment analysis on CC terms
paramsCC <- new("GOHyperGParams", geneIds=GOstatsGenes, universeGeneIds=geneUniverse,
annotation=annotationDB, ontology="CC", pvalueCutoff=GOpValueCutoff, conditional=TRUE,
testDirection="over")

## perform GO enrichment analysis on MF terms
paramsMF <- new("GOHyperGParams", geneIds=GOstatsGenes, universeGeneIds=geneUniverse,
annotation=annotationDB, ontology="MF", pvalueCutoff=GOpValueCutoff, conditional=TRUE,
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

########################################## gagePathwayAnalysis ########################################## 

library(gage)
library(pathview)
pathwayDirectory <- paste(outputDirectory, "pathway_images", sep="/")
system(paste("mkdir ", pathwayDirectory, sep = ""))
setwd(pathwayDirectory)

## gets the list of entrez genes and their associated peak scores from the input file
geneData <- unique(read.delim(annotationFile)[,c("entrezID", "score")])
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
write.table(pathwayData, paste(runName, "full_pathway_list", sep="_"), sep="\t", row.names=FALSE, quote = FALSE)

## get enriched pathways with a p value less than the cutoff value
sel <- gageOutput$greater[,"p.val"] < pathway_pValueCutoff & !is.na(gageOutput$greater[,"p.val"])
pathIds <- substr(rownames(gageOutput$greater)[sel], 1, 8)

## Use pathview to create output images
runPathview <- function(pid){
	print(pid)
	pathview(gene.data=dataMatrix, pathway.id=pid, species=KEGGspeciesCode)
	system(paste("rename 's/^/", runName, "_/' ", pid, "*", sep=""))
	}
pathList <- sapply(pathIds, function(pid) runPathview(pid))

########################################## getPeakSequences ########################################## 

pythonCommand <- paste("python ", keplerBinLocation, "/getPeakSequences.py", sep="")
peakSequenceFile=paste(outputDirectory, "/", runName, "_peakSequences.fa", sep="")
system(paste(pythonCommand, MACSpeakFile, genomeDirectory, peakSequenceFile, sep=" "))

########################################## runDREME ################################################## 

dremeOutputDirectory <- paste(outputDirectory, "/", runName, "_DREME_output", sep="")
system(paste("dreme -oc", dremeOutputDirectory, "-p", peakSequenceFile, "-t", maxDremeRuntime, "-m", maxMotifsToGenerate, "-e", DREMEeValueCutoff, "-desc -l", sep=" "))


