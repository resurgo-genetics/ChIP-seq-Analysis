## indexBAMTKv07.R
## Version 7.0 Monday, December 8, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program to index BAM files.
# Output is a sorted .bam file and a .bai (bam index file) in the same directory as the input file, with the 
# same name as the input file but ending in the outputSuffix (default _sorted).

## DEPENDENCIES: SAMtools 
## Run speed: Usually under 5 minutes, but depends on inputFile size.

## DIRECTIONS:
# Change the appropriate parameters (inputFile, runName)
# Run the program.

## Example Input:
# inputFile <- "/home/tylerk/massInputData/GATA1_mappedReads_sorted.bam"
# outputDirectory <- "/home/user/outputData/indexBAM"
# outputSuffix <- "_Test1"

################# Set Parameters Here ##################

## The input file and output directory must be complete file paths, with no trailing /
inputFile <- "/home/tylerk/massInputData/GATA1_mappedReads_sorted.bam" # Must be a .bam file
runName <- "_test2" # The suffix of the output file name

################# End Parameter Block ##################

#################     INDEX BAM       ##################
outputFileName <- gsub(".bam", runName, inputFile)
system(paste("samtools sort", inputFile, outputFileName, sep=" "))
sortedBAMfile <- paste(outputFileName, ".bam", sep="")
system(paste("samtools index", sortedBAMfile, sep=" "))
system(paste("mv ", sortedBAMfile, ".bai ", gsub(".bam", ".bai", sortedBAMfile), sep =""))