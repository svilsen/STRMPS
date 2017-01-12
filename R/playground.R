###### Playground ######
if (FALSE) {
    library("microbenchmark")
    library("tidyverse")
    dilution_files <- list.files("~/MPS/Data/Illumina/ForenSeq_Dilution/ForenSeq 2016-01-22 Sensitivity panel A/FastQ/", full.names = T)
    file <- readFastq(dilution_files[1])

    load("~/AAU/PhD/Code/IlluminaForenSeq/PanelAnBSNPs/flankingRegionsForenSeqPanelA.RData")
    ## STR
    flankingRegionsForenSeqPanelASTR <- flankingRegionsForenSeqPanelA %>% filter(Marker != "DYS389II", Type != "SNP")

    # One mismatch
    caseSTR <- identifySTRRegions(reads = file, flankingRegions = flankingRegionsForenSeqPanelASTR, numberOfMutation = 1,
                                    control = identifySTRRegions.control(includeReverseComplement = T, matchPatternMethod = "mclapply"))

    includedMarkersSTR <- which((flankingRegionsForenSeqPanelASTR$Marker %in% names(caseSTR$identifiedMarkersSequencesUniquelyAssigned)))

    sortedSTR <- sapply(names(caseSTR$identifiedMarkersSequencesUniquelyAssigned), function(m) which(m == flankingRegionsForenSeqPanelASTR$Marker))
    motifLengthsSTR <- flankingRegionsForenSeqPanelASTR$MotifLength[sortedSTR]
    typeSTR <- flankingRegionsForenSeqPanelASTR$Type[sortedSTR]

    coverageSTR <- stringCoverage(extractedReadsListObject = caseSTR, control = stringCoverage.control(trace = T, motifLength = motifLengthsSTR, Type = typeSTR))

    noiseSTR <- identifyNoise(coverageSTR, colBelief = "Coverage", thresholdSignal = 0.05, thresholdAbsoluteLowerLimit = 1)
    nonNoisyMarkers <- unique(do.call(rbind, lapply(seq_along(noiseSTR), function(nA) noiseSTR[[nA]] %>% mutate(Marker = names(noiseSTR[nA]))))$Marker)

    # Two mismatch
    caseSTR2 <- identifySTRRegions(reads = file, flankingRegions = flankingRegionsForenSeqPanelASTR, numberOfMutation = 2,
                                        control = identifySTRRegions.control(includeReverseComplement = T, matchPatternMethod = "mclapply"))

    includedMarkersSTR2 <- which((flankingRegionsForenSeqPanelASTR$Marker %in% names(caseSTR2$identifiedMarkersSequencesUniquelyAssigned)))

    sortedSTR2 <- sapply(names(caseSTR2$identifiedMarkersSequencesUniquelyAssigned), function(m) which(m == flankingRegionsForenSeqPanelASTR$Marker))
    motifLengthsSTR2 <- flankingRegionsForenSeqPanelASTR$MotifLength[sortedSTR2]
    typeSTR2 <- flankingRegionsForenSeqPanelASTR$Type[sortedSTR2]

    coverageSTR2 <- stringCoverage(extractedReadsListObject = caseSTR2,
                                        control = stringCoverage.control(trace = T, motifLength = motifLengthsSTR2, Type = typeSTR2))

    noiseSTR2 <- identifyNoise(coverageSTR2, colBelief = "Coverage", thresholdSignal = 0.05, thresholdAbsoluteLowerLimit = 1)
    nonNoisyMarkers2 <- unique(do.call(rbind, lapply(seq_along(noiseSTR2), function(nA) noiseSTR2[[nA]] %>% mutate(Marker = names(noiseSTR2[nA]))))$Marker)

    temp <- flankingRegionsForenSeqPanelASTR
    temp$Found2 <- temp$Found1 <- FALSE
    temp$Found1[temp$Marker %in% nonNoisyMarkers] <- TRUE
    temp$Found2[temp$Marker %in% nonNoisyMarkers2] <- TRUE
    temp[which(temp$Found2 != temp$Found1), ]

    sum(sapply(coverageSTR2, function(cA) sum(cA$Coverage))) / length(sread(file))
    sum(sapply(coverageSTR, function(cA) sum(cA$Coverage))) / length(sread(file))

    ## SNP's
    flankingRegionsForenSeqPanelASNP <- flankingRegionsForenSeqPanelA %>% filter(Marker != "DYF387S1", Type == "SNP")

    caseSNP <- identifySTRRegions(reads = file, flankingRegions = flankingRegionsForenSeqPanelASNP, numberOfMutation = 1,
                                        control = identifySTRRegions.control(includeReverseComplement = T, matchPatternMethod = "mclapply"))

    includedMarkersSNP <- which((flankingRegionsForenSeqPanelASNP$Marker %in% names(caseSNP$identifiedMarkersSequencesUniquelyAssigned)))
    motifLengthsSNP <- flankingRegionsForenSeqPanelASNP$MotifLength[includedMarkersSNP]
    typeSNP <- flankingRegionsForenSeqPanelASNP$Type[includedMarkersSNP]

    coverageSNP <- stringCoverage(extractedReadsListObject = caseSNP,
                                    control = stringCoverage.control(trace = T, motifLength = motifLengthsSNP, Type = typeSNP))

    noiseSNP <- identifyNoise(coverageSNP, colBelief = "Coverage", thresholdSignal = 0.05, thresholdAbsoluteLowerLimit = 1)
    sum(sapply(coverageSNP, function(cA) sum(cA$Coverage))) / length(sread(file))

    ## Collected
    flankingRegionsForenSeqPanelACollected <- flankingRegionsForenSeqPanelA %>% filter(Marker != "DYF387S1")

    caseCollected <- identifySTRRegions(reads = file, flankingRegions = flankingRegionsForenSeqPanelACollected, numberOfMutation = 1,
                                        control = identifySTRRegions.control(includeReverseComplement = T, matchPatternMethod = "mclapply"))

    includedMarkersCollected <- which((flankingRegionsForenSeqPanelACollected$Marker %in% names(caseCollected$identifiedMarkersSequencesUniquelyAssigned)))

    sortedCollected <- sapply(names(caseCollected$identifiedMarkersSequencesUniquelyAssigned), function(m) which(m == flankingRegionsForenSeqPanelACollected$Marker))
    motifLengthsCollected <- flankingRegionsForenSeqPanelACollected$MotifLength[sortedCollected]
    typeCollected <- flankingRegionsForenSeqPanelACollected$Type[sortedCollected]

    coverageCollected <- stringCoverage(extractedReadsListObject = caseCollected,
                                        control = stringCoverage.control(trace = T, motifLength = motifLengthsCollected, Type = typeCollected))

    noiseCollected <- identifyNoise(coverageCollected, colBelief = "Coverage", thresholdSignal = 0.05, thresholdAbsoluteLowerLimit = 1)
    nonNoisyMarkers <- unique(do.call(rbind, lapply(seq_along(noiseCollected), function(nA) noiseCollected[[nA]] %>% mutate(Marker = names(noiseCollected[nA]))))$Marker)
    sum(sapply(coverageSNP, function(cA) sum(cA$Coverage))) / length(sread(file))

    ## Stutter
    genotypeSTR <- getGenotype(coverageSTR, colBelief = "Coverage", thresholdSignal = 0.05, thresholdHeterozygosity = 0.5, thresholdAbsoluteLowerLimit = 1)

    genotypeCoverageSTR <- mergeGenotypeStringCoverage(coverageSTR, genotypeSTR)
    stutterSTR <- findStutter(genotypeCoverageSTR, motifLengthsSTR)
    do.call(rbind, stutterSTR) %>% filter(!is.na(NeighbourAllele))

    genotypeSTRAverageCoverage <- sapply(genotypeSTR, function(s) mean(s$Coverage, na.rm = T))
    genotypeSTR2AverageCoverage <- sapply(getGenotype(coverageSTR2, colBelief = "Coverage", thresholdSignal = 0.05, thresholdHeterozygosity = 0.5, thresholdAbsoluteLowerLimit = 1), function(s) mean(s$Coverage, na.rm = T))
    comparingAverageAlleleCoverage <- tibble(Marker = flankingRegionsForenSeqPanelASTR$Marker) %>%
        left_join(tibble(Marker = names(genotypeSTRAverageCoverage), STRAverageCoverage = genotypeSTRAverageCoverage), by = "Marker") %>%
        left_join(tibble(Marker = names(genotypeSTRAverageCoverage), STR2AverageCoverage = genotypeSTR2AverageCoverage), by = "Marker")


}
