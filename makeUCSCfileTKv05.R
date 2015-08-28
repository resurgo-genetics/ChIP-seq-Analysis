## makeUCSCfileTKv05.R
## Version 5.0 Wednesday, December 10th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, Mark Bieda.

## DESCRIPTION: 
# This program creates a custom UCSC browser track out of a sorted .bam file. 

## Run Speed: Usually less than 5 minutes.

## DIRECTIONS:
# Note: This program will NOT work with RStudio.
# Change the appropriate parameters (inputFile, outputDirectory, runName, trackName, trackResolution, maxFileSize)
# Run the program.

## Output is one file in the outputDirectory with the name format: runName_UCSCtrack.gz

## Example Input:
# inputFile <- "/home/tylerk/massInputData/GATA1_mappedReads_sorted.bam" #sorted bam file
# outputDirectory <- "/home/tylerk/massTestResults/makeUCSC"
# runName <- "Test1" # runName will be used as a prefix in the output file names
# trackName <- "testTrack"
# trackResolution <- 10
# maxFileSize <- 1e8

## Example Run from command line:
# Rscript makeUCSCfile.R
 
###################################### SET PARAMETERS HERE ##########################################
## Input Options:
#The input file and output directory must be complete file paths, with no trailing /

inputFile <- "/home/tylerk/massInputData/GATA1_mappedReads_sorted.bam" # Input must be a sorted bam file
outputDirectory <- "/home/tylerk/massTestResults/makeUCSC" # The directory where output files will be located
runName <- "Test4" # runName will be used as a prefix in the output file names

## UCSC file options:
trackName <- "GATA1 Test Track " # Specify a name for the UCSC track to be created, this name will appear once it is uploaded to the UCSC genome browser.
trackResolution <- 10 # Numeric resolution of the test track
maxFileSize <- 1e8 # Numeric maximum filesize of the test track

###################################### END PARAMETER BLOCK ##########################################

########################################## createUCSCfile ########################################### 
tagDirectory <- paste(outputDirectory, "tagDirectory", sep="/")
outputFileName <- paste(outputDirectory, "/", runName, "_UCSCtrack", sep="")

system(paste("makeTagDirectory", tagDirectory, inputFile, sep=" "))
system(paste("makeUCSCfile", tagDirectory, "-o", outputFileName, "-name", trackName, "-res", trackResolution, "-color 0,0,0 -fsize", maxFileSize, sep=" "))

