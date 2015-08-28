## distancetoTSS02NCTK06.R
## Version 6.0 Dec 2 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This program takes in a set of peak data in .bed format
# uses the specified TxDb to calculate distances between each peak and its nearest TSS

## Note: If you are not using the default TxDb: See the bioconductor website 
# for a list of available TxDbs and install the package to R 
# http://bioconductor.org/packages/3.0/data/annotation/

## Run Speed: Usually less than 30 seconds.

## DIRECTIONS:
# Change the appropriate parameters (readFile, outputDirectory, chromosome, startPosition, endPosition, runName, titlePrefix).
# Run the program.
# Output is one .txt file into the outputDirectory

## Example Input:
# inputFile <- "/home/user/inputData/GATA1_peaks.bed"
# outputDirectory <- "/home/user/outputData"
# runName <- "test01"
# transcriptionDB <- "TxDb.Hsapiens.UCSC.hg19.knownGene"

## Example Use of Different TxDB for Mmusculus assembly mm10:
# Visit Bioconductor website link above, find correct TxDB package in list, click link and 
# Install by running in R:
# source("http://bioconductor.org/biocLite.R")
# biocLite("TxDb.Mmusculus.UCSC.mm10.ensGene")
# Then be sure to change transcriptDB parameter to "TxDb.Mmusculus.UCSC.mm10.ensGene"

##########################   Set Parameters Here ####################################

## The input file and output directory must be complete file paths
inputFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" #Input file must be a set of peak data in .bed file format, no trailing /
outputDirectory <- "/home/tylerk/massTestResults/distanceToTSS" #The output directory, no trailing /

runName <- "distanceToTSStk1" # runName will be used as a prefix in the output file names

transcriptionDB <- "TxDb.Hsapiens.UCSC.hg19.knownGene"  ##Specify the correct TxDb package for your species of interest and genome assembly.

##################### End of Parameter Block ###############################

## Load required libraries
library(transcriptionDB, character.only=TRUE)
library(GenomicRanges)
library(rtracklayer)
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

peaks <- import.bed(con=inputFile, asRangedData=F) # read in peaks as a Grange object

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
