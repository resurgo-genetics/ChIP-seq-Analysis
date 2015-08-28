## graphRange05NCTKv02.R
## Version 5.2 Thursday, December 11th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program to create a graph of the read coverage of a specified interval
# Version 4.0 Dec 2 2014
# Created By: Tyler Kolisnik, Nathan Cormier and Mark Bieda.

## Run Speed: Usually less than 30 seconds.

## DIRECTIONS:
# Change the appropriate parameters (readFile, outputDirectory, chromosome, startPosition, endPosition, runName, titlePrefix).
# Run the program.
# Output is one .jpeg file into the outputDirectory with the name 

## Example Input:
# readFile <- "/home/user/inputData/GATA1_mappedReads_sorted.bam"
# outputDirectory <- "/home/user/outputData"
# chromosome <- "chr21"
# startPosition <- 9825881 
# endPosition <- 982755
# runName <- "test01"
# titlePrefix <- "My GATA1 Graph"
# useCustomYrange <- "yes"
# yMin <- 5
# yMax <- 25

##################### Parameters #############################
readFile <- "/home/tylerk/massInputData/GATA1_mappedReads_sorted.bam" # Must be a sorted, indexed .bam file with associated .bai file, no trailing /
outputDirectory <- "/home/tylerk/massTestResults/graphRange" # Must be an existing directory, no trailing /

chromosome <- "chr21" # Must be written in the same form as found in the .bam file, chr# is the standard

startPosition <- 9825881 #Numeric Start Position of specified interval within chromosome, Cannot exceed 536870912 (limitation of readGAlignments)
endPosition <- 9827255 #Numeric Start Position of specified interval within chromosome, Cannot exceed 536870912 (limitation of readGAlignments)

runName <- "test01" # runName will be used as a prefix in the output file names
titlePrefix <- "Test Title" # The Title of the graph created


## If it is set to yes it will use the axis limits specified, if set to any other value
## then the axis will be created automatically to fit the input data 
useCustomYrange <- "yes"  #Set to "yes" to specify custom Y axis restrictions, set to "no" for automatic
yMin <- 5 #Numeric minimum of y Axis
yMax <- 25 #Numeric maximum of y Axis

##################### End of Parameter Block #############################


## Load required libraries
library(GenomicRanges)
library(GenomicAlignments)

## Function to create the plots with the automatic y axis range
plotCovDefault <- function(mycov=cov, mychr=1, mypos=c(1,1000), mymain="Coverage", ...) {
	op <- par(mar=c(8,5,6,1))
	plot(as.numeric(mycov[[mychr]][mypos[1]:mypos[2]]), type="l",
	lwd=1, col="blue", ylab="Tag Count", main=mymain, xlab="", xaxt="n", ...)
	axis(1, las = 2, at=seq(1,mypos[2]-mypos[1], length.out= 10),
	labels=as.integer(seq(mypos[1], mypos[2], length.out= 10)))
	mtext(paste("Chromosome position (", chromosome, ")", sep=""), side=1, line=6)
	par(op)
}

## Function to create the same plot but with a fixed y axis
plotCovCustom <- function(mycov=cov, mychr=1, mypos=c(1,1000), mymain="Coverage", ...) {
	op <- par(mar=c(8,5,6,1))
	plot(as.numeric(mycov[[mychr]][mypos[1]:mypos[2]]), type="l",
	lwd=1, col="blue",ylim=c(yMin, yMax), ylab="Tag Count", main=mymain, xlab="", xaxt="n", ...)
	axis(1, las = 2, at=seq(1,mypos[2]-mypos[1], length.out= 10),
	labels=as.integer(seq(mypos[1], mypos[2], length.out= 10)))
	mtext(paste("Chromosome position (", chromosome, ")", sep=""), side=1, line=6)
	par(op)
}

chr <- chromosome
start <- as.numeric(startPosition)
end <- as.numeric(endPosition)

## Reads in all sequence reads from the bam file in the specified interval
which <- GRanges(chr, IRanges(start, end))
param <- ScanBamParam(which=which)
reads <- readGAlignments(readFile, param=param)

## Removes extra chr length data from the Granges object
seqlevels(reads) <- seqlevels(reads)[seqlevels(reads) %in% as.character(chr)]
cov <- coverage(reads) # get the total coverage at each base in the specified range
rangeMax <- sum(cov[[1]]@lengths)
if (end > rangeMax) {
	end <- rangeMax
}
## Creates and outputs the coverage graph
jpeg(file=paste(outputDirectory, "/", runName, "_", chr, "_",start,"-",end,".jpeg",sep=""))
if (useCustomYrange!="yes") {
	plot <- plotCovDefault(mycov=cov, mychr=as.character(chr), mypos=c(start,end), 
	mymain=paste(titlePrefix, " ", chr, ":",as.character(start),"-",as.character(end), sep=""))
} else {
	plot <- plotCovCustom(mycov=cov, mychr=as.character(chr), mypos=c(start,end), 
	mymain=paste(titlePrefix, " ", chr, ":",as.character(start),"-",as.character(end), sep=""))
}
dev.off()
