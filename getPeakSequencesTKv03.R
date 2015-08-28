## getPeakSequences03TK.R
## Version 3.0 Wednesday, December 10th, 2014
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program to that takes in a .bed peak file and gives the sequences of the peaks.

## DIRECTIONS:
# Change the appropriate parameters (inputFile, outputDirectory, keplerBinLocation, genomeDirectory, runName)
# Run the program.
# It outputs 1 file into the working directory in the format runName_peakSequences.fa

## Example Input:
# inputFile <- "/home/user/inputData/GATA1_peaks.bed" 
# outputDirectory <- "/home/user/outputData"
# keplerBinLocation <- "/home/user/CWF1/keplerBin"
# genomeDirectory <- "/home/user/CWF1/keplerBin/genomes/hg19"

## Example if using a species or assembly other than hg19:
# The genomeDirectory parameter for genomic data must also be changed, hg19 and mm18 are included by default but for any others they must be downloaded ftp://ftp.ncbi.nih.gov/genomes/ and added to the keplerBin.
# Change the genome directory to the new species, example for mm18:
# genomeDirectory <- "/home/user/CWF1/keplerBin/genomes/mm18"

############################ Parameters ###########################

## The input file, output directory, and keplerBinLocation must be complete file paths with no trailing /
inputFile <- "/home/tylerk/massInputData/GATA1_peaks.bed" # This is the Peak file outputted from MACS (usually ends in "_peaks.bed")
outputDirectory <- "/home/tylerk/massTestResults/getPeakSequences" # This is the location the output files will go to
keplerBinLocation <- "/home/tylerk/CWF1/keplerBin" # The location of the keplerBin where the python scripts are located
genomeDirectory <- "/home/tylerk/CWF1/keplerBin/genomes/hg19" # The directory of your genome file, within the kepler bin. hg19=human genome 19, Change if using different species or assembly!
runName <- "Test1" # runName will be used as a prefix in the output file names

######################End of Parameter Block#######################


pythonCommand <- paste("python ", keplerBinLocation, "/getPeakSequences.py", sep="")
peakSequenceFile=paste(outputDirectory, "/", runName, "_peakSequences.fa", sep="")
system(paste(pythonCommand, inputFile, genomeDirectory, peakSequenceFile, sep=" "))

