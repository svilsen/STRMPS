all_in_one.control <- function(identifySTRRegions = identifySTRRegions.control(), stringCoverage = stringCoverage.control(), identifyNoise = TRUE) {
    res <- list(nrOfMutations = nrOfMutations, identifySTRRegionsControl = identifySTRRegionsControl, stringCoverage = stringCoverage,
                identifyNoise = identifyNoise)
    return(res)
}

#' @title All-in-one function
#'
#' @description The function takes an input file and performs the entire analysis workflow described in (ADD REF).
#' The function creates a series of objects needed for further analyses.
#' An output folder can be provided to store the objects as \code{.RData}-files.
#'
#' @param input a path to a \code{.fastq}-file, or directory containing multiple files (NOT IMPLEMENTED) all of which should be analysed.
#' @param output a directory where output-files are stored.
#'
#' @return The function returns all relevant objects created which might be useful for further analysis.
#'
#' @export
all_in_one <- function(input, output = NULL, saveCheckpoint = TRUE, continueCheckpoint = NULL, flankingRegions = NULL, nrThreads = 2, control = all_in_one.control()) {
    if (is.null(output) & saveCheckpoint)
        stop("'saveCheckpoint' is TRUE with no output folder specified.")

    read_path <- input
    readfile <- readFastq(read_path)
    if (saveCheckpoint)
        save(readfile, file = paste(output, "/", "fastqfileloaded", ".RData", sep = ""))

    if (is.null(flankingRegions))
        load("flankingRegions")

    identifiedSTRs <- identifySTRRegions(reads = read_path, flankingRegions = flankingRegions, nrOfMutations = control$nrOfMutations, control = control$identifySTRRegionsControl)
    if (saveCheckpoint)
        save(identifiedSTRs, file = paste(output, "/", "identifiedSTRRegions", ".RData", sep = ""))

    stringCoverageList <- stringCoverage(extractedReadsListObject = identifiedSTRs, control = control$stringCoverage)
    if (saveCheckpoint)
        save(stringCoverageList, file = paste(output, "/", "stringCoverageList", ".RData", sep = ""))

    if (control$identifyNoise) {
        ### create noise threshold function needs to be added...
        createdThresholdSignal <- 0.01
        identifiedNoise <- identifyNoise(stringCoverageList, "Coverage", thresholdSignal = createdThresholdSignal)
        noiseStringCoverageList <- mergeNoiseStringCoverage(stringCoverageList, identifiedNoise)
        if (saveCheckpoint)
            save(noiseStringCoverageList, file = paste(output, "/", "noiseStringCoverageList", ".RData", sep = ""))
    }

    return(0)
}
