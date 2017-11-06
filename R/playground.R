###### Playground ######
if (FALSE) {
    library("microbenchmark")
    library("tidyverse")
    dilution_files <- list.files("~/MPS/Data/Illumina/ForenSeq_Dilution/ForenSeq 2016-01-22 Sensitivity panel A/FastQ/", full.names = T)
    mixture_file <- readFastq("~/MPS/Data/Illumina/ForenSeq_Mixtures/R701-A501_S1_L001_R1_001.fastq")
    reference_files <- list.files("~/MPS/Data/Illumina/ForenSeq_Danes/FastQ 20150905_CHU_FTAdup1/", full.names = T)

    file <- readFastq(reference_files[19])

    load("~/AAU/PhD/Code/IlluminaForenSeq/PanelAnBSNPs/flankingRegionsForenSeqAll.RData")
    flankingRegionsForenSeqSTRs <- flankingRegionsForenSeqAll %>% filter(Type == "AUTOSOMAL") %>% rename(ForwardFlank = Forward, ReverseFlank = Reverse)
    flankingRegionsForenSeqSTRs <- flankingRegionsForenSeqSTRs %>% mutate(ForwardShift = 0, ReverseShift = 0)

    reads = file; flankingRegions = flankingRegionsForenSeqSTRs; numberOfMutation = 1; control = identifySTRRegions.control(matchPatternMethod = "mclapply")

    seqs <- sread(reads)
    qual <- quality(reads)

    colID <- if(is.null(control$colList)) .getCols(names(flankingRegions)) else colList
    flankSizes <- apply(flankingRegions[, c(colID$forwardCol, colID$reverseCol)], 2, nchar)

    colList = colID; numberOfThreads = control$numberOfThreads; removeEmptyMarkers = control$removeEmptyMarkers

    forwardFlank <- unname(unlist(flankingRegions[, colList$forwardCol]))
    reverseFlank <- unname(unlist(flankingRegions[, colList$reverseCol]))

    seqsCharacter <- as.character(seqs)
    maxNumberOfMismatch <- Biostrings:::normargMaxMismatch(numberOfMutation)

    ss = paste(forwardFlank[1], "AATGAATGAATGAATG", reverseFlank[1], sep = "", collapse = "")
    vmatchMultiPatternSeqAn(forwardFlank, reverseFlank, ss, c(1), 1)

    ################


    microbenchmark(caseSTR <- identifySTRRegions(reads = file, flankingRegions = flankingRegionsForenSeqSTRs, numberOfMutation = 1,
                                                 control = identifySTRRegions.control(matchPatternMethod = "mclapply")),
                   times = 1)

    microbenchmark(caseSTR_SEQAN <- identifySTRRegions(reads = file, flankingRegions = flankingRegionsForenSeqSTRs, numberOfMutation = 1,
                                                       control = identifySTRRegions.control(matchPatternMethod = "seqan")),
                   times = 1)


    caseSTR$n_reads
    caseSTR_SEQAN$n_reads

    includedMarkersSTR <- which((flankingRegionsForenSeqSTRs$Marker %in% names(caseSTR_SEQAN$identifiedMarkersSequencesUniquelyAssigned)))
    sortedSTR <- sapply(names(caseSTR_SEQAN$identifiedMarkersSequencesUniquelyAssigned), function(m) which(m == flankingRegionsForenSeqSTRs$Marker))
    motifLengthsSTR <- flankingRegionsForenSeqSTRs$MotifLength[sortedSTR]
    typeSTR <- flankingRegionsForenSeqSTRs$Type[sortedSTR]

    microbenchmark(coverageSTR <- STRMPS:::.extractedReadsList.stringCoverage(extractedReadsListObject = caseSTR_SEQAN, control = stringCoverage.control(trace = T, motifLength = motifLengthsSTR, Type = typeSTR, includeLUS = FALSE)), times = 1)


    coverageSTR$AMELOGENIN
    names(caseSTR_SEQAN$identifiedMarkersSequencesUniquelyAssigned)

    genotypeSTR <- getGenotype(coverageSTR, thresholdHeterozygosity = 0.3)
    coverageGenotypeSTR <- mergeGenotypeStringCoverage(coverageSTR, genotypeSTR)

    noiseSTR <- identifyNoise(coverageSTR, colBelief = "Coverage", thresholdSignal = 0.05, thresholdAbsoluteLowerLimit = 1)
    nonNoisyMarkers <- unique(do.call(rbind, lapply(seq_along(noiseSTR), function(nA) noiseSTR[[nA]] %>% mutate(Marker = names(noiseSTR[nA]))))$Marker)


    ##

}
