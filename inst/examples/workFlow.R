# ## Basic work-flow of the STRMPSWorkflow function
# library("Biostrings")
# library("ShortRead")
# library("tidyverse")
#
# read_path <- "..."
#
# flankingRegionsPath <- system.file("flankingregions", "", package = "STRMPS")
# flankingRegions <-
#     STRMPS:::.loadRData(paste(flankingRegionsPath,
#                               "flankingRegionsForenSeqSTRsShifted.RData", sep = "/"))
#
# control = workflow.control()
#
# # Read the file into memory
# readfile <- readFastq(read_path)
#
# # Identify the STR's of the file
# identifiedSTRs <- identifySTRRegions(reads = readfile, flankingRegions = flankingRegions,
#                                      numberOfMutation = control$numberOfMutations,
#                                      control = identifySTRRegions.control(
#                                          numberOfThreads = control$numberOfThreads))
#
# # Restructuring
# sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned),
#                                 function(m) which(m == flankingRegions$Marker))
#
# # Aggregate the strings
# stringCoverageList <-
#     stringCoverage(extractedReadsListObject = identifiedSTRs,
#                    control = stringCoverage.control(
#                        motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
#                        Type = flankingRegions$Type[sortedIncludedMarkers],
#                        numberOfThreads = control$numberOfThreads,
#                        trace = control$internalTrace,
#                        simpleReturn = control$simpleReturn))
#
# # In addition string can be trimmed and the number of strings can be reduced.
