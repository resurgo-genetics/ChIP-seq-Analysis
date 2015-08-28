## peakStatsTKv07.R
## Version 7.0 Wednesday, December 10th, 2014
## Created By: Tyler Kolisnik, Nathan Cormier and Mark Bieda.

## DESCRIPTION:
# This program takes a set of peak data in .bed format and 
# calculates several statistics about the input peak file

## Run Speed: Usually under 1 minute.

## DIRECTIONS: 
# Change the appropriate parameters within the parameter block in the initial part of this script. 
# Run the program.
# Output is a single .txt file with the name format: runName_peakStats.txt and will be located in the outputDirectory.

## Example Input:
# peakFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" 
# outputDirectory <- "/home/tylerk/massTestResults/peakStats" 
# genomeSize <- "/home/tylerk/massInputData/chrSizes_hg19.txt" 
# runName <- "massTest2" 
# broadCutOff <- 10000 

## Example Partial Contents of Chromosome Size .txt File, NOTE: Must not contain a header:
# chr1  249250621
# chr2	243199373
# chr3	198022430
# chr4	191154276
# chr5	180915260
# chr6	171115067
# chr7	159138663
# chrX	155270560

## Example use of genomeSize using actual size of genome:
# genomeSize <- 30000000
##########################   SET PARAMETERS HERE   ####################################

## The input file and output directory must be complete file paths with no trailing /
peakFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" # Should be a .bed file of peak data.
outputDirectory <- "/home/tylerk/massTestResults/peakStats" # The directory where all output files will be located
genomeSize <- "/home/tylerk/massInputData/chrSizes_hg19.txt" ## Either enter a path for a .txt file of chromosome sizes or Numerically specify the size of the genome (actual number of bases) for the species you are studying 
runName <- "massTest2" # runName will be used as a prefix in the output file names
broadCutOff <- 10000 # Numeric number of bases a peak must have to be considered a broad peak

##########################  END OF PARAMETER BLOCK   ####################################

outputFile <- paste(outputDirectory, "/", runName, "_peakStats.txt", sep="")
dataIn <- read.table(peakFile)

totalPeaks <- nrow(dataIn) # count number of peaks
dataIn$peakLength <- dataIn$V3 - dataIn$V2
smallestPeak <- min (dataIn$peakLength) # find size of the smallest peak
largestPeak <- max (dataIn$peakLength) # find size of the largest peak
meanLength <- mean(dataIn$peakLength) # find the mean peak length
SD <- sd(dataIn$peakLength) # find the standard deviation of the peak length
medianLength <- median(dataIn$peakLength) # find the median peak length

## calculates the total genome size from the genomeSize input file if necessary
if (is.numeric(genomeSize)) {  
	totalBases <- genomeSize
} else {
	lengths <- read.table(genomeSize)
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
