data("stringCoverageList")
data("genotypeList")
stringCoverageGenotypeList <- mergeGenotypeStringCoverage(stringCoverageList, genotypeList)

data("noiseList")
stringCoverageNoiseList <- mergeNoiseStringCoverage(stringCoverageList, noiseList)
