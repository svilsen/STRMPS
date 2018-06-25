# Regions identified using 'identifySTRs()'
data("identifiedSTRs")

# Limiting and restructuring
sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned),
                                function(m) which(m == flankingRegions$Marker))

# Aggregate the strings
stringCoverage(extractedReadsListObject = identifiedSTRs,
               control = stringCoverage.control(
                   motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
                   Type = flankingRegions$Type[sortedIncludedMarkers],
                   numberOfThreads = control$numberOfThreads,
                   trace = control$internalTrace,
                   simpleReturn = control$simpleReturn))
