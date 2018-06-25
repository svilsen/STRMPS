library("Biostrings")
library("ShortRead")
library("tidyverse")

read_path <- system.file('sampleFiles', package = 'STRMPS')

flankingRegionsPath <- system.file("flankingregions", "", package = "STRMPS")
flankingRegions <-
    STRMPS:::.loadRData(paste(flankingRegionsPath,
                              "flankingRegionsForenSeqSTRsShifted.RData", sep = "/")) %>%
    filter(Type == "AUTOSOMAL")

control = workflow.control()

# Read the file into memory
readfile <- readFastq(read_path)
sread(readfile)
quality(readfile)

# Identify the STR's of the file
identifiedSTRs <- identifySTRRegions(reads = readfile, flankingRegions = flankingRegions,
                                     numberOfMutation = control$numberOfMutations,
                                     control = identifySTRRegions.control(
                                         numberOfThreads = control$numberOfThreads,
                                         includeReverseComplement = FALSE))

