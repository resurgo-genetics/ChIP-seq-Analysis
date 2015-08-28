## bamtobedTKv06.R
## Version 6.0 Friday, December 12th, 2014
## Created By: Tyler Kolisnik, Nathan Cormier and Mark Bieda.

## Description:
# This is a program to convert a BAMfile to a BED file.

## Run Speed: Usually less than 5 minutes, but dependent on inputFile size. 

## NOTE: This program will not work in RStudio, you must run it from the terminal using Rscript programName.R

## DEPENDENCIES: BEDtools

## DIRECTIONS:
# Change the appropriate parameters (outputDirectory, inputFile).
# Run the program.
# Outputs a .bam file in the outputDirectory.

## EXAMPLE INPUT:
# inputFile <- "/home/user/inputData/Gata1_mappedReads.bam" 
# outputDirectory <- "/home/user/outputData"

## EXAMPLE RUN:
## Rscript bamtobedTKv06.R



#################################   SET PARAMETERS HERE   #####################################

## The input file and output directory must be complete file paths, must be in quotes and do not use a trailing /.
inputFile <- "/home/tylerk/bamtest/Tylerv1AcerH3KFullTest1_mappedReads.bam" #T his should be the .bam file to be converted 
outputDirectory <- "/home/tylerk/bamtest" # This is the directory where output files will be located

###############################   END PARAMETER BLOCK   #######################################

########################################## BAMtoBED ########################################## 

BAMfile <- inputFile
fileName <- tail(strsplit(BAMfile, "/")[[1]], n=1)
BEDfile <- paste(outputDirectory, gsub(".bam", ".bed", fileName), sep="/")
system(paste("bedtools bamtobed -i",BAMfile, ">",BEDfile, sep=" "))

