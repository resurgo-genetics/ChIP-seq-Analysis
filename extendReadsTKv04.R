## extendReadsTKv04.R
## Version 4.0 Thursday, December 11th, 2014
## Created by: Nathan Cormier, Tyler Kolisnik, and Mark Bieda

## DESCRIPTION:
# This is a program that takes in a .bed mapped sequence read file and creates a new file which extends the reads
# in a specified size in a direction appropriate for the strand.

## Run Speed: Usually less than 2 minutes, depending on inputFile size. 

## DIRECTIONS:
# Change the appropriate parameters (inputFile, outputDirectory, keplerBinLocation, runName, fragmentSize )
# Run the program.
# Outputs one .bed file with the name format: runName_extendedReads_fragmentSize.bed into the outputDirectory.

## Example Input:
# inputFile <- "/home/user/inputData/GATA1_mappedReads.bed" 
# outputDirectory <- "/home/user/outputData/extendReads"
# keplerBinLocation <- "/home/user/CWF1/keplerBin"
# runName <- "Test1" 
# fragmentSize <- 500 

##########################   SET PARAMETERS HERE ####################################

## The input file and output directory must be complete file paths, with no trailing /

inputFile <- "/home/tylerk/massInputData/GATA1_mappedReads.bed" #This should be a .bed file of mapped reads
outputDirectory <- "/home/tylerk/massTestResults/extendReads" #This is the directory where output files will be located
keplerBinLocation <- "/home/tylerk/CWF1/keplerBin" #The location of the Kepler Bin
runName <- "Test1" # runName will be used as a prefix in the output file names
fragmentSize <- 500 # Numeric estimation of fragment size after sonication

##########################   END PARAMETER BLOCK ####################################

##########################     EXTEND READS     #####################################

	pythonCommand <- paste("python ", keplerBinLocation, "/tagExtender2.py", sep="") 
	extendedReads <- paste(outputDirectory, "/", runName, "_extendedReads_", fragmentSize, ".bed", sep="")
	system(paste(pythonCommand,  inputFile, fragmentSize, extendedReads, sep=" "))

print("extendReads completed successfully")
