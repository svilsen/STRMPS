library("dplyr")

# The object returned by merging a stringCoverageList-Object
# and a genotypeList-Object.
data("stringCoverageGenotypeList")

stutterList <- findStutter(stringCoverageGenotypeList)
stutterTibble <- do.call("rbind", stutterList) %>%
    filter(!is.na(Genotype))

stutterTibble$BlockLengthMissingMotif
stutterTibble$NeighbourRatio
