## TSSdensityPlot03NCTKv02.R
## Version 3.2 Wednesday, December 10, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda. 

## DESCRIPTION:
# This is a program that takes in the TSStagDensities.txt output file from Homer and creates a
# plot showing average sequence read coverage relative to the TSS

## Run Speed: Usually less than 30 seconds.

## DIRECTIONS:
# Change the appropriate parameters (inputFile, outputDirectory, runName).
# Run the program.
# Output is one .pdf file into the outputDirectory with the name format: runName_densityPlot.pdf

## EXAMPLE INPUT:
# readFile <- "/home/user/inputData/GATA1_TSStagDensities.txt"
# outputDirectory <- "/home/user/outputDirectory"
# runName <- "test01"

#####################################   SET PARAMETERS HERE   #####################################

## The input file and output directory must be complete file paths with no trailing /
inputFile <- "/home/tylerk/massInputData/GATA1_TSStagDensities.txt"
outputDirectory <- "/home/tylerk/massTestResults/TSSdensity"
runName <- "massTest1" # runName will be used as a prefix in the output file names

####################################  END OF PARAMETER BLOCK  ######################################

########################################## TSSdensityPlot ########################################## 
## errors occur when running in R studio, but works fine when run as an Rscript. This is an Rstudio problem, not a script problem


## creates a graph of the input data
dataIn <- read.table(inputFile, skip = 1)
pdf(paste(outputDirectory, "/", runName, "_densityPlot.pdf", sep=""))
plot(x=dataIn$V1, y=dataIn$V2, type="l", xlab="Distance from TSS", ylab="Tag Density (per bp per TSS)")
dev.off()
