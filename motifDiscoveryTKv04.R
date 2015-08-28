## motifDiscoveryTKv04.R
## Version 4.0 Wednesday, December 11th, 2014.
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program that uses DREME to perform de novo motif discovery on a peak sequence file.

## Run Speed: Long, Dependent on maxDremeRunTime variable, at least 1 hour is recommended. 

## DIRECTIONS:
# Change the appropriate parameters (inputFile, outputDirectory, runName, DREMEeValueCutoff, maxMotifsToGenerate, maxDremeRunTime)
# Run the program. Output is a folder with the file name format: runName_DREME_output. It contains 3 files labeled dreme.html, dreme.txt, dreme.xml. 

## EXAMPLE INPUT:
# inputFile <- "/home/user/inputData/GATA1_peakSequences.fa" 
# outputDirectory <- "/home/user/outputData/motifDiscovery"
# runName <- "Test1" 
# DREMEeValueCutoff <- 0.001
# maxMotifsToGenerate <- 10
# maxDremeRuntime <- 0  

##################################### Set Parameters Here  ###################################

## The input file and output directory must be complete file paths, with no trailing /
inputFile <- "/home/tylerk/massInputData/GATA1_peakSequences.fa" #This must be a peak sequence file in .fasta or .fa format
outputDirectory <- "/home/tylerk/massTestResults/motifDiscovery" #Directory where output files will be located
runName <- "massTest" #runName will be used as a prefix in the output file names

## DREME options
DREMEeValueCutoff <- 0.001 #Numeric E value cut off for motif generation
maxMotifsToGenerate <- 10 #Numeric maximum number of motifs to be generated
maxDremeRuntime <- 0  #Numeric maximum motif search time in seconds. Set to 0 for no limit (it will stop when the maximum number of motifs are found). At least 3600 (1 hour) recommended.

#################################### END OF PARAMETER BLOCK ####################################

########################################## runDREME ##################################

dremeOutputDirectory <- paste(outputDirectory, "/", runName, "_DREME_output", sep="")
if (maxDremeRuntime > 0){
  system(paste("dreme -oc", dremeOutputDirectory, "-p", inputFile, "-t", maxDremeRuntime, "-m", maxMotifsToGenerate, "-e", DREMEeValueCutoff, "-desc -l", sep=" "))
} else {
  system(paste("dreme -oc", dremeOutputDirectory, "-p", inputFile, "-m", maxMotifsToGenerate, "-e", DREMEeValueCutoff, "-desc -l", sep=" "))
}

