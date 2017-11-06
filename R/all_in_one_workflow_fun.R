workflow.control <- function(numberOfMutations = 1, numberOfThreads = 4, createdThresholdSignal = 0.05, thresholdHomozygote = 0.4,
                             internalTrace = FALSE, simpleReturn = TRUE, identifyNoise = TRUE, identifyStutter = TRUE,
                             identifySTRRegions = identifySTRRegions.control, stringCoverageControl = stringCoverage.control) {
    res <- list(numberOfMutations = numberOfMutations, numberOfThreads = numberOfThreads, createdThresholdSignal = createdThresholdSignal,
                thresholdHomozygote = thresholdHomozygote, internalTrace = internalTrace, simpleReturn = simpleReturn,
                identifyNoise = identifyNoise, identifyStutter = identifyStutter,
                identifySTRRegionsControl = identifySTRRegionsControl, stringCoverageControl = stringCoverageControl)
    return(res)
}

#' @title Workflow function
#'
#' @description The function takes an input file and performs the entire analysis workflow described in (ADD REF).
#' The function creates a series of objects needed for further analyses.
#' An output folder can be provided to store the objects as \code{.RData}-files.
#'
#' @param input A path to a \code{.fastq}-file.
#' @param output A directory where output-files are stored.
#' @param continueCheckpoint Choose a checkpoint to continue from in the workflow. If NULL the function will run the entire workflow.
#' @param flankingRegions The flanking regions used to identify markers and trim alleles.
#' @param control Function controlling non-crucial parameters and other control functions.
#' @return If 'output' not provided the function simply returns the stringCoverageList-object.
#' If an output is provided the function will store ALL created objects at the output-path, but nothing is returned.
workflow <- function(input, output = NULL, continueCheckpoint = NULL, flankingRegions = NULL, control = workflow.control()) {
    if (!is.null(output)) {
        dirExists <- dir.exists(output)
        if (!dirExists) {
            dir.create(output)
        }

        saveCheckpoint = TRUE
    }

    if (is.null(continueCheckpoint))
        continueCheckpoint <- FALSE

    # Reads file
    fileExists <- file.exists(paste(output, "_", "fastqfileloaded", ".RData", sep = ""))
    if (continueCheckpoint & fileExists) {
        load(paste(output, "_", "fastqfileloaded", ".RData", sep = ""))
    }
    else {
        read_path <- input
        readfile <- readFastq(read_path)
        if (saveCheckpoint)
            save(readfile, file = paste(output, "_", "fastqfileloaded", ".RData", sep = ""))
    }

    if (is.null(flankingRegions))
        load("flankingRegions")

    # Identified sequences list
    fileExists <- file.exists(paste(output, "_", "identifiedSTRRegions", ".RData", sep = ""))
    if (continueCheckpoint & fileExists) {
        load(paste(output, "_", "identifiedSTRRegions", ".RData", sep = ""))
    }
    else {
        identifiedSTRs <- identifySTRRegions(reads = read_path, flankingRegions = flankingRegions, nrOfMutations = control$numberOfMutations,
                                             control = control$identifySTRRegionsControl(numberOfThreads = control$numberOfThreads, trace = control$internalTrace))
        if (saveCheckpoint)
            save(identifiedSTRs, file = paste(output, "_", "identifiedSTRRegions", ".RData", sep = ""))
    }

    # String coverage list
    fileExists <- file.exists(paste(output, "/", "stringCoverageList", ".RData", sep = ""))
    if (continueCheckpoint & fileExists) {
            load(paste(output, "_", "stringCoverageList", ".RData", sep = ""))
    }
    else {
        sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned), function(m) which(m == flankingRegions$Marker))
        stringCoverageList <- stringCoverage(extractedReadsListObject = identifiedSTRs,
                                             control = control$stringCoverageControl(motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
                                                                                     Type = flankingRegions$Type[sortedIncludedMarkers],
                                                                                     numberOfThreads = control$numberOfThreads, trace = control$internalTrace,
                                                                                     simpleReturn = control$simpleReturn))
        if (saveCheckpoint)
            save(stringCoverageList, file = paste(output, "_", "stringCoverageList", ".RData", sep = ""))

    }

    # Noise list
    if (control$identifyNoise) {
        fileExists <- file.exists(paste(output, "_", "noiseStringCoverageList", ".RData", sep = ""))
        if (continueCheckpoint & fileExists) {
            load(paste(output, "_", "noiseStringCoverageList", ".RData", sep = ""))
        }
        else {
            identifiedNoise <- identifyNoise(stringCoverageList, "Coverage", thresholdSignal = control$createdThresholdSignal)
            noiseStringCoverageList <- mergeNoiseStringCoverage(stringCoverageList, identifiedNoise)

            if (saveCheckpoint)
                save(noiseStringCoverageList, file = paste(output, "_", "noiseStringCoverageList", ".RData", sep = ""))
        }
    }

    # Stutters
    if (control$identifyStutter) {
        fileExists <- all(file.exists(paste(output, "_", "genotypeList", ".RData", sep = "")),
                          file.exists(paste(output, "_", "stutterTibble", ".RData", sep = "")))
        if (continueCheckpoint & fileExists) {
            load(paste(output, "_", "genotypeList", ".RData", sep = ""))
            load(paste(output, "_", "stutterTibble", ".RData", sep = ""))
        }
        else {
            sortedIncludedMarkers <- sapply(names(stringCoverageList), function(m) which(m == flankingRegions$Marker))
            t_H <- control$thresholdHomozygote + (1 - control$thresholdHomozygote - 0.01) * as.numeric(flankingRegions$Type[sortedIncludedMarkers] == "Y")

            genotypeList = getGenotype(stringCoverageList, thresholdHeterozygosity = t_H, thresholdSignal = 0.01, thresholdAbsoluteLowerLimit = 1)
            stringCoverageGenotypeList <- mergeGenotypeStringCoverage(stringCoverageList, genotypeList)

            stutterTibble <- do.call(rbind, findStutter(stringCoverageGenotypeList, searchDirection = -1)) %>% filter(!is.na(NeighbourAllele))

            if (saveCheckpoint) {
                save(genotypeList, file = paste(output, "_", "genotypeList", ".RData", sep = ""))
                save(stutterTibble, file = paste(output, "_", "stutterTibble", ".RData", sep = ""))
            }
        }
    }

    if (is.null(output) & !saveCheckpoint)
        return(stringCoverageList)
}

#' @title Workflow function for multiple files.
#'
#' @description The function wrap-function for the \link{workflow}-function allowing for more than one input file.
#'
#' @param input A vector of paths to \code{.fastq}-files.
#' @param output A directory where output-files are stored (forced).
#' @param flankingRegions The flanking regions used to identify markers and trim alleles.
#' @param control Function controlling non-crucial parameters and other control functions.
#' @return The function will store ALL created objects at the output-path; nothing is returned.
workflowList <- function(input, output = NULL, flankingRegions = NULL, control = workflow.control()) {
    if (is.null(output)) {
        output = "File"
    }
    if(length(output) == 1) {
        output = paste(output, seq(1, length(input)), sep = "_")
    }

    for (ff in input) {
        workflow(input = input[ff], output = output[ff], flankingRegions = flankingRegions, control = control)
    }
}
