## runMacsv12TK.R
## Version 12, Friday, December 12th, 2014.
## Created By: Nathan Cormier, Tyler Kolisnik, and Mark Bieda.

## DESCRIPTION:
# This is a program to run MACS (Model-based Analysis of ChIP-Seq) on a .bed or .bam sequence reads file.

## Run Speed: Can take from 30 seconds to 10 minutes depending on dataset size. 

## DIRECTIONS:
# Change the appropriate parameters (workingdirStr, sequenceReads, controlFileName, runName, useControl)
# Run the program.

## Output varies with input file type and if a control file is used or not.
# If using a .bed file:
# Output is: a _peaks.bed file, a _peaks.xls, a _summits.bed file, and if using a control, a _negative_peaks.xls file.
# If using a .bam file:
# Output is: a model.r file, a _peaks.bed file, a _peaks.xls, a _summits.bed file, and if using a control, a _negative_peaks.xls file.

## Example Input 1:
# workingDirectory <- "/home/tylerk/massTestResults/runMacs" 
# sequenceReads <- "/home/tylerk/massTestResults/runMacs/GATA1_peaks.bed" 
# controlFileName <- "" 
# runName <- "Test1"
# useControl <- "no"

## Example Input 2:
# outputDirectory <- "/home/tylerk/massTestResults/runMacs" 
# sequenceReads <- "/home/tylerk/massTestResults/runMacs/GATA1_mappedReads.bam" 
# controlFileName <- "/home/tylerk/massTestResults/runMacs/GATA1_controlreads.bam" 
# runName <- "Test2"
# useControl <- "yes"

################################ Set Parameters Here ####################################
sequenceReads <- "/home/tylerk/H3Kinput/H3K_mappedReads.bed" #Full path of your sequence reads file, may be .bed or .bam
outputDirectory <- "/home/tylerk/H3Koutput/runMacs" #Full path of working directory, do not include a trailing / 
useControl <- "no" #Specify yes if a Control file will be used, or no if a control file won't be used
controlFileName <- "" #Full path of your control file, if used, must be same filetype as sequenceReads, if not used set to ""
runName <- "H3KTest" #The output suffix that will be attached to the names of your output files.
############################# End of Parameter Block ####################################

runDir <- paste(outputDirectory, "/", runName, sep="") #Parameters for part of the string used in invoking the system command
setwd(outputDirectory) #Sets the working directory

fileSuffix <- tail(strsplit(sequenceReads, ".b")[[1]], n=1) #Checks to see if a .bam file or a .bed file was inputted, split returns what comes after the last .b at the end of the file name
if (fileSuffix =="am"){ 
  filetype <- ".bam"
} else if (fileSuffix=="ed"){
  filetype <- ".bed"
} else {
  print("Error: Invalid file type, must be .bed or .bam")
}


if (useControl=="yes" && filetype==".bed"){
   outputString <- paste("macs14 -t", sequenceReads, "-c", controlFileName, "-f BED -g hs -n", runDir)
   system(outputString)
} else if (useControl=="no" && filetype==".bed"){
   outputString2 <- paste("macs14 -t", sequenceReads, "-f BED -g hs -n", runDir) 
   system(outputString2)
} else if (useControl=="yes" && filetype==".bam"){
   outputString3 <- paste("macs14 -t", sequenceReads, "-c", controlFileName, "-f BAM -g hs -n", runDir)
   system(outputString3)
} else if (useControl=="no" && filetype==".bam"){
   outputString4 <- paste("macs14 -t", sequenceReads, "-f BAM -g hs -n", runDir) 
   system(outputString4)
} else {
  print("Error: Incorrect input. Please specify the correct filetype and if a control will be used. ")
}

