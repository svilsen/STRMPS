# Regions identified using 'identifySTRs()'
data("identifiedSTRs")

# Limiting and restructuring
sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned),
                                function(m) which(m == flankingRegions$Marker))

# Aggregate the strings
stringCoverage(extractedReadsListObject = identifiedSTRs,
               motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
               Type = flankingRegions$Type[sortedIncludedMarkers],
               control = stringCoverage.control(numberOfThreads = 1,
                                                trace = FALSE,
                                                simpleReturn = TRUE))
