.isInstalled <- function(pkg) is.element(pkg, installed.packages()[,1])

#' Workflow default options
#'
#' Control object for workflow function returning a list of default parameter options.
#'
#' @param numberOfMutations The maximum number of mutations (base-calling errors) allowed during flanking region identification.
#' @param numberOfThreads The number of threads used by mclapply (stuck at '2' on windows).
#' @param createdThresholdSignal Noise threshold.
#' @param thresholdHomozygote Homozygote threshold for genotype identiication.
#' @param internalTrace Show trace.
#' @param simpleReturn TRUE/FALSE: Should the regions be aggregated without including flanking regions?
#' @param identifyNoise TRUE/FALSE: Should noise be identified.
#' @param identifyStutter TRUE/FALSE: Should stutters be identified.
#' @param flankingRegions The flanking regions used to identify the STR regions. If 'NULL' a default set is loaded and used.
#' @param useSTRaitRazor TRUE/FALSE: Should the STRaitRazor command line tool (only linux is implemented) be used for flanking region identification.
#' @param trimRegions TRUE/FALSE: Should the identified regions be further trimmed.
#'
#' @return List of default of options.
workflow.control <- function(numberOfMutations = 1, numberOfThreads = 4, createdThresholdSignal = 0.05, thresholdHomozygote = 0.4,
                             internalTrace = FALSE, simpleReturn = TRUE, identifyNoise = FALSE, identifyStutter = FALSE,
                             flankingRegions = NULL, useSTRaitRazor = TRUE, trimRegions = TRUE) {
    if (useSTRaitRazor) {
        if ((!(.isInstalled("STRaitRazoR"))) | (tolower(Sys.info()['sysname']) != "linux")) {
            useSTRaitRazor = FALSE
        }
    }

    res <- list(numberOfMutations = numberOfMutations, numberOfThreads = numberOfThreads, createdThresholdSignal = createdThresholdSignal,
                thresholdHomozygote = thresholdHomozygote, internalTrace = internalTrace, simpleReturn = simpleReturn,
                identifyNoise = identifyNoise, identifyStutter = identifyStutter, flankingRegions = flankingRegions,
                useSTRaitRazor = useSTRaitRazor, trimRegions = trimRegions)
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
#' @param control Function controlling non-crucial parameters and other control functions.
#' @return If 'output' not provided the function simply returns the stringCoverageList-object.
#' If an output is provided the function will store ALL created objects at the output-path, i.e. nothing is returned.
STRMPSWorkflow <- function(input, output = NULL, continueCheckpoint = NULL, control = workflow.control()) {
    if (!is.null(output)) {
        dirExists <- dir.exists(output)
        if (!dirExists) {
            dir.create(output)
        }

        saveCheckpoint = TRUE
    }
    else {
        saveCheckpoint = FALSE
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

    if (is.null(control$flankingRegions)) {
        flankingRegionsPath <- system.file("flankingregions", "", package = "STRMPS")
        flankingRegions <- .loadRData("flankingRegionsPath/flankingRegionsForenSeqSTRsShifted.RData")
    }
    else {
        flankingRegions <- control$flankingRegions
    }

    if (control$useSTRaitRazor) {
        # STRaitRazor v3
        fileExists <- file.exists(paste(output, "/", "stringCoverageList", ".RData", sep = ""))
        if (continueCheckpoint & fileExists) {
            stringCoverageList = .loadRData(paste(output, "_", "stringCoverageList", ".RData", sep = ""))
        }
        else {
            stringCoverageList <- STRaitRazoR::STRaitRazorSTRMPS(read_path, control = STRaitRazoR::STRaitRazorSTRMPS.control(numberOfThreads = control$numberOfThreads))
            save(stringCoverageList, file = paste(output, "_", "stringCoverageList", ".RData", sep = ""))
        }
    }
    else {
        # Identified sequences list
        fileExists <- file.exists(paste(output, "_", "identifiedSTRRegions", ".RData", sep = ""))
        if (continueCheckpoint & fileExists) {
            load(paste(output, "_", "identifiedSTRRegions", ".RData", sep = ""))
        }
        else {
            identifiedSTRs <- identifySTRRegions(reads = read_path, flankingRegions = flankingRegions, numberOfMutation = control$numberOfMutations,
                                                 control = identifySTRRegions.control(numberOfThreads = control$numberOfThreads))
            if (saveCheckpoint)
                save(identifiedSTRs, file = paste(output, "_", "identifiedSTRRegions", ".RData", sep = ""))
        }

        # String coverage list
        fileExists <- file.exists(paste(output, "/", "stringCoverageList", ".RData", sep = ""))
        if (continueCheckpoint & fileExists) {
            stringCoverageList = .loadRData(paste(output, "_", "stringCoverageList", ".RData", sep = ""))
        }
        else {
            sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned), function(m) which(m == flankingRegions$Marker))
            stringCoverageList <- stringCoverage(extractedReadsListObject = identifiedSTRs,
                                                 control = stringCoverage.control(motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
                                                                                  Type = flankingRegions$Type[sortedIncludedMarkers],
                                                                                  numberOfThreads = control$numberOfThreads,
                                                                                  trace = control$internalTrace,
                                                                                  simpleReturn = control$simpleReturn))
            if (saveCheckpoint)
                save(stringCoverageList, file = paste(output, "_", "stringCoverageList", ".RData", sep = ""))

        }
    }

    if (control$trimRegions) {
        fileExists <- file.exists(paste(output, "/", "stringCoverageList", ".RData", sep = ""))
        if (continueCheckpoint & fileExists) {
            stringCoverageList = .loadRData(paste(output, "_", "stringCoverageListTrimmed", ".RData", sep = ""))
        }
        else {
            stringCoverageListTrimmed <- lapply(stringCoverageList, function(ss) {
                flankingRegions_m <- flankingRegions %>% filter(Marker == unique(ss$Marker))

                if ((flankingRegions_m$ForwardShift == 0) & (flankingRegions_m$ReverseShift == 0)) {
                    res <- ss
                }
                else {
                    res <- ss %>% rename(ExpandedRegion = Region) %>%
                        mutate(BasePairs = nchar(ExpandedRegion), AdjustedBasePairs = BasePairs - flankingRegions_m$Offset,
                               Region = str_sub(ExpandedRegion, start = flankingRegions_m$ForwardShift + 1, end = - flankingRegions_m$ReverseShift - 1)) %>%
                        select(-ExpandedRegion, BasePairs, AdjustedBasePairs) %>% group_by(Marker, Type, Region, MotifLength) %>%
                        summarise(Allele = max(Allele), Coverage = sum(Coverage)) %>% ungroup() %>%
                        select(Marker, Type, Allele, MotifLength, Region, Coverage) %>%
                        arrange(Allele, Region)
                }

                return(res)
            })

            class(stringCoverageListTrimmed) <- "stringCoverageList"
            save(stringCoverageListTrimmed, file = paste(output, "_", "stringCoverageListTrimmed", ".RData", sep = ""))
            stringCoverageList <- stringCoverageListTrimmed
        }
    }

    # Noise list
    if (control$identifyNoise) {
        fileExists <- file.exists(paste(output, "_", "noiseStringCoverageList", ".RData", sep = ""))
        if (!(continueCheckpoint & fileExists)) {
            identifiedNoise <- identifyNoise(stringCoverageList, "Coverage", thresholdSignal = control$createdThresholdSignal)
            noiseStringCoverageList <- mergeNoiseStringCoverage(stringCoverageList, identifiedNoise)

            if (saveCheckpoint)
                save(noiseStringCoverageList, file = paste(output, "_", "noiseStringCoverageList", ".RData", sep = ""))
        }
    }

    # Stutters
    if (control$identifyStutter) {
        fileExists <- file.exists(paste(output, "_", "genotypeList", ".RData", sep = ""))
        if (continueCheckpoint & fileExists) {
            genotypeList = .loadRData(paste(output, "_", "stutterTibble", ".RData", sep = ""))
        }
        else {
            sortedIncludedMarkers <- sapply(names(stringCoverageList), function(m) which(m == flankingRegions$Marker))
            t_H <- control$thresholdHomozygote + (1 - control$thresholdHomozygote - 0.01) * as.numeric(flankingRegions$Type[sortedIncludedMarkers] == "Y")

            genotypeList = getGenotype(stringCoverageList, thresholdHeterozygosity = t_H, thresholdSignal = 0.01, thresholdAbsoluteLowerLimit = 1)

            if (saveCheckpoint) {
                save(genotypeList, file = paste(output, "_", "genotypeList", ".RData", sep = ""))
            }
        }

        fileExists <- file.exists(paste(output, "_", "stutterTibble", ".RData", sep = ""))
        if (!(continueCheckpoint & fileExists)) {
            stringCoverageGenotypeList <- mergeGenotypeStringCoverage(stringCoverageList, genotypeList)

            stutterTibble <- do.call(rbind, findStutter(stringCoverageGenotypeList)) %>% filter(!is.na(NeighbourAllele))

            if (saveCheckpoint) {
                save(stutterTibble, file = paste(output, "_", "stutterTibble", ".RData", sep = ""))
            }
        }
    }

    if (is.null(output) & !saveCheckpoint)
        return(stringCoverageList)
}

#' @title Batch wrapper for the workflow function
#'
#' @description The function takes an input directory and performs the entire analysis workflow described in (ADD REF).
#' The function creates a series of objects needed for further analyses and stores them at the output location.
#'
#' @param input A directory where fastq input-files are stored.
#' @param output A directory where output-files are stored.
#' @param continueCheckpoint Choose a checkpoint to continue from in the workflow. If NULL the function will run the entire workflow.
#' @param control Function controlling non-crucial parameters and other control functions.
#' @return If 'output' not provided the function simply returns the stringCoverageList-object.
#' If an output is provided the function will store ALL created objects at the output-path, i.e. nothing is returned.
STRMPSWorkflowBatch <- function(input, output, continueCheckpoint = NULL, control = workflow.control()) {
    files <- list.files(input, full.names = T, recursive = TRUE)
    files <- files[nchar(gsub("fastq", "", files)) < nchar(files)]

    run_names <- stringr::str_split(list.files(input, recursive = TRUE), "/")
    dir_names <- lapply(run_names, function(ss) gsub("[- ]", "_", ss[-length(ss)]))
    file_names <- gsub("[- ]", "_", sapply(stringr::str_split(sapply(run_names, function(ss) ss[length(ss)]), ".fastq"), function(ss) ss[1]))

    dir_exists <- dir.exists(output)
    if (!dir_exists) {
        dir.create(output)
    }

    for (i in seq_along(files)) {
        dir_names_i <- dir_names[[i]]

        for (j in seq_along(dir_names_i)) {
            dir_exists <- dir.exists(paste(output, dir_names_i[1:j], collapse = "/", sep = "/"))
            if (!dir_exists) {
                dir.create(paste(output, dir_names_i[1:j], collapse = "/", sep = "/"))
            }
        }

        output_i <- paste(output, dir_names_i, file_names[i], collapse = "/", sep = "/")
        STRMPSWorkflow(files[i], output_i, continueCheckpoint, control)
    }
}

#' @title Collect stutters files
#'
#' @description Collects all stutter files.
#'
#' @param stutterDirectory The out most directory containing all stutter files to be collected.
#' @param storeCollection TRUE/FALSE: Should the collected tibble be stored? If 'FALSE' the tibble is returned.
#'
#' @return If 'storeCollection' is TRUE nothing is returned, else the stutter collection is returned.
STRMPSWorkflowCollectStutters <- function(stutterDirectory, storeCollection = TRUE) {
    files <- list.files(stutterDirectory, pattern = "*_stutterTibble.RData", full.names = TRUE, recursive = TRUE)

    collectedStutter <- vector("list", length(files))
    for (i in seq_along(files)) {
        files_i <- files[i]
        collectedStutter[[i]] <- .loadRData(files_i) %>% mutate(Sample = i)
    }

    collectedStutter <- bind_rows(collectedStutter)
    if (storeCollection) {
        save(collectedStutter, file = paste(stutterDirectory, "collectedStutterTibble.RData", sep = "/"))
    }
    else {
        return(collectedStutter)
    }
}
