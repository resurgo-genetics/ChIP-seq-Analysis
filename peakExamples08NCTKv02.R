## peakExamples08NCTKv02.R
## Version 8.2 Thursday, December 11th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This program takes in both a set of peak data in .bed format and mapped sequence reads in .bam format
# Creates a specified number of example peak plots
# also creates several plots of non-peak regions as a comparison

## Note: This program must be run from the command line and not using RStudio.

## Run Speed: Usually less than 30 seconds.

## DIRECTIONS:
# Change the appropriate parameters 
# Run the program.
# Output is in the form .jpeg files into the outputDirectory, amount is user-specified by numberofExamples

## Example Input:
# peakFile <- "/home/user/inputData/GATA1_peaks.bed"
# readFile <-  "/home/user/inputData/GATA1_mappedReads_sorted.bam"
# outputDirectory <- "/home/user/outputData"
# runName <- "test01"
# useCustomYrange <- "yes"
# yMin <- 0
# yMax <- 45
# numberOfExamples <- 10 
# rangeAroundPeak <- 10000 
# nonPeakRegionSize <- 5000 

##########################   SET PARAMETERS HERE  ####################################

## The input file and output directory must be complete file paths, with no trailing /
peakFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" #Set of peak Data in .bed file format
readFile <- "/home/tylerk/massInputData/GATA1_mappedReads_sorted.bam" #Set of sorted mapped Reads in .bam file format
outputDirectory <- "/home/tylerk/massTestResults/peakExamples" #Directory where output files will be located

runName <- "runTest" #runName will be used as a prefix in the output file names
titlePrefix <- "Title Test" #The title prefix will be added to the title of the graph creates

## If it is set to yes it will use the axis limits specified, if set to any other value
## then the axis will be created automatically to fit the input data 
useCustomYrange <- "yes" # Set to "yes" to specify custom Y axis restrictions, set to "no" for automatic
yMin <- 0  # Numeric minimum of y axis
yMax <- 45  # Numeric maximum of y axis

numberOfExamples <- 10 # Number of peak example and non-peak example coverage plots to create
rangeAroundPeak <- 10000 # Number of bases up and downstream of the peak to include in the plot
nonPeakRegionSize <- 5000 # Total number of bases to include in each non-peak example coverage plot

############################ END OF PARAMETER BLOCK  ################################

library(GenomicRanges)
library(GenomicAlignments)

## function to create the plots with the automatic y axis range
plotCovDefault <- function(mycov=cov, mychr=1, mypos=c(1,1000), mymain="Coverage", ...) {
	op <- par(mar=c(8,5,6,1))
	plot(as.numeric(mycov[[mychr]][mypos[1]:mypos[2]]), type="l",
	lwd=1, col="blue", ylab="Tag Count", main=mymain, xlab="", xaxt="n", ...)
	axis(1, las = 2, at=seq(1,mypos[2]-mypos[1], length.out= 10),
	labels=as.integer(seq(mypos[1], mypos[2], length.out= 10)))
	mtext(paste("Chromosome position (", mychr, ")", sep=""), side=1, line=6)
	par(op)
}

## function to create the same plot but with a fixed y axis
plotCovCustom <- function(mycov=cov, mychr=1, mypos=c(1,1000), mymain="Coverage", ...) {
	op <- par(mar=c(8,5,6,1))
	plot(as.numeric(mycov[[mychr]][mypos[1]:mypos[2]]), type="l",
	lwd=1, col="blue",ylim=c(yMin, yMax), ylab="Tag Count", main=mymain, xlab="", xaxt="n", ...)
	axis(1, las = 2, at=seq(1,mypos[2]-mypos[1], length.out= 10),
	labels=as.integer(seq(mypos[1], mypos[2], length.out= 10)))
	mtext(paste("Chromosome position (", mychr, ")", sep=""), side=1, line=6)
	par(op)
}

peaks <- read.table(peakFile)
## sort peaks by score and take top scoring examples +- the specified range around them
peaks <- peaks[order(peaks$V5, decreasing=T),]
POI <- peaks[1:numberOfExamples,]
range <- rangeAroundPeak
POI$V2 <- POI$V2-range
POI$V3 <- POI$V3+range

## read in sequences from the .bam file which overlap the peak range
which <- GRanges(POI$V1, IRanges(POI$V2, POI$V3))
param <- ScanBamParam(which=which)
reads <- readGAlignments(readFile, param=param)

## removes extra chr length data from the Granges object
seqlevels(reads) <- seqlevels(reads)[seqlevels(reads) %in% as.character(POI$V1)]
cov <- coverage(reads) # get the total coverage at each base in the specified range

if (useCustomYrange!="yes") { 
	## checks to make sure the graph does not extend past the end of its chromosome
	for (i in 1:numberOfExamples){
		jpeg(file=paste(outputDirectory, "/", runName, "_peak",i,".jpeg",sep=""))
		rangeMax <- sum(cov[[as.character(POI$V1[i])]]@lengths)
		if (POI$V2[i] < 0) {
			POI$V2[i] <- 0
		}
		if (POI$V3[i] > rangeMax) {
			POI$V3[i] <- rangeMax
		}
		## create the plots using automatic y axis scales for each plot
		plot <- plotCovDefault(mycov=cov, mychr=as.character(POI$V1[i]), mypos=c(POI$V2[i],POI$V3[i]), mymain=paste(titlePrefix, POI$V1[i],POI$V2[i],"-",POI$V3[i], sep=" "))
		dev.off()
	}
} else { 
	## checks to make sure the graph does not extend past the end of its chromosome
	for (i in 1:numberOfExamples){
		jpeg(file=paste(outputDirectory, "/", runName, "_peak",i,".jpeg",sep=""))
		rangeMax <- sum(cov[[as.character(POI$V1[i])]]@lengths)
		if (POI$V2[i] < 0) {
			POI$V2[i] <- 0
		}
		if (POI$V3[i] > rangeMax) {
			POI$V3[i] <- rangeMax
		}
		## create the plots using the set y axis
		plot <- plotCovCustom(mycov=cov, mychr=as.character(POI$V1[i]), mypos=c(POI$V2[i],POI$V3[i]), mymain=paste(titlePrefix, POI$V1[i],POI$V2[i],"-",POI$V3[i], sep=" "))
		dev.off()
	}

}

###################################################################
## creates a list of locations which are between peaks and likely have very few reads
## the regions with the least (but not zero) sequence reads are plotted
## This script assumes that the peak file is fairly large (several thousand peaks)
## and it may not work correctly if the file is very small

## pick a set of peaks evenly spaced across the peak file
peaks <- peaks[order(peaks$V1), ]
peakChecks <- numberOfExamples*10
negativeList <- data.frame(chr=(character()), start=(integer()), end=(integer()))
cut <- floor(nrow(peaks)/peakChecks)
regionSize <- nonPeakRegionSize/2

## Select a region in between each chosen peak and the previous peak in the input file
for (i in 1:peakChecks){
	x <- peaks[i*cut, ]
	y <- peaks[(i*cut)-1, ]
	if (x[1] == y[1]){
		midpoint = floor((y$V3 + x$V2)/2)
		newRow <- x[1:3]
		newRow[2] <- midpoint - regionSize
		newRow[3] <- midpoint + regionSize
		negativeList <- rbind(negativeList, newRow)
	}

}

## read in sequences from the .bam file which overlap the range
which <- GRanges(negativeList$V1, IRanges(negativeList$V2, negativeList$V3))
param <- ScanBamParam(which=which)
reads <- readGAlignments(readFile, param=param)	
seqlevels(reads) <- seqlevels(reads)[seqlevels(reads) %in% as.character(negativeList$V1)]
cov <- coverage(reads) # get the total coverage at each base in the specified range
totalReads <- c()

## checks to make sure the graph does not extend past the end of its chromosome
for (i in 1:nrow(negativeList)){
	rangeMax <- sum(cov[[as.character(negativeList$V1[i])]]@lengths)
	if (negativeList$V2[i] < 0) {
		negativeList$V2[i] <- 0
	}
	if (negativeList$V3[i] > rangeMax) {
		negativeList$V3[i] <- rangeMax
	}
	totalReads <- append(totalReads, sum(cov[[as.character(negativeList$V1[i])]][negativeList$V2[i]:negativeList$V3[i]]))
}

## Sorts possible negative regions by total coverage, removing any regions with no sequence reads
## the regions with the lowest total sequence reads are graphed
negativeList[,"reads"] <- totalReads
negativeList <- negativeList[order(negativeList$reads), ]
negativeList <- negativeList[!(negativeList$reads==0), ]
negativeList <- negativeList[1:numberOfExamples, ]
if (useCustomYrange!="yes") {	# create the plots using automatic y axis scales for each plot
	for (i in 1:numberOfExamples){
		jpeg(file=paste(outputDirectory, "/", runName, "_NonPeakRegion",i,".jpeg",sep=""))
		plot <- plotCovDefault(mycov=cov, mychr=as.character(negativeList$V1[i]), mypos=c(negativeList$V2[i],negativeList$V3[i]), 
		mymain=paste(titlePrefix, negativeList$V1[i],negativeList$V2[i],"-",negativeList$V3[i], sep=" "))
		dev.off()
	}
} else {
	for (i in 1:numberOfExamples){ # create the plots using the set y axis
		jpeg(file=paste(outputDirectory, "/", runName, "_NonPeakRegion",i,".jpeg",sep=""))
		plot <- plotCovCustom(mycov=cov, mychr=as.character(negativeList$V1[i]), mypos=c(negativeList$V2[i],negativeList$V3[i]), 
		mymain=paste(titlePrefix, negativeList$V1[i],negativeList$V2[i],"-",negativeList$V3[i], sep=" "))
		dev.off()
	}

}
