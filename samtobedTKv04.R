## samtobedTKv04.R
## Version 4.0 Thursday, December 11th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program to convert a SAM file to a BAMfiles then convert the BAMfile to a BED file.

## NOTE: This program will not work in RStudio, you must run it from the terminal using Rscript programName.R

## DEPENDENCIES: SAMtools, BEDtools

## Run Speed: Usually less than 5 minutes, but dependent on inputFile size. 

## DIRECTIONS:
# Change the appropriate parameters (outputDirectory, inputFile).
# Run the program.
# Output is a .bam file (optional) and a .bed file in the outputDirectory.

## EXAMPLE INPUT:
# inputFile <- "/home/user/inputData/Gata1_mappedReads.sam" 
# outputDirectory <- "/home/user/outputData" 
# keepBam <- "yes" 

## EXAMPLE RUN:
# Rscript samtobedv03TK.R

#############################   SET PARAMETERS HERE   #####################################

## The input file and output directory must be complete file paths, with no trailing /

inputFile <- "/home/tylerk/samtest/Tylerv1AcerH3KFullTest1_mappedReads.sam" #This should be a .sam file
outputDirectory <- "/home/tylerk/samtest" #The directory where the .bed file and optional .bam file will be outputted to
keepBam <- "yes" #Specify yes or no as to keep the intermediate Bam file


##############################   END PARAMETER BLOCK   #######################################

########################################## SAMtoBED ########################################## 

SAMfile <- inputFile
fileName <- tail(strsplit(SAMfile, "/")[[1]], n=1)
BAMfile <- paste(outputDirectory, gsub(".sam", ".bam", fileName), sep="/")
BEDfile <- paste(outputDirectory, gsub(".sam", ".bed", fileName), sep="/")
system(paste("samtools view -bS",SAMfile, ">",BAMfile, sep=" "))
system(paste("bedtools bamtobed -i",BAMfile, ">",BEDfile, sep=" "))

if (!keepBam=="yes" ){
  system(paste("rm", BAMfile))
}

