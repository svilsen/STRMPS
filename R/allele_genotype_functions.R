stringCoverage.control <- function(motifLength = 4, Type = "AUTOSOMAL", simpleReturn = TRUE, includeLUS = FALSE, numberOfThreads = 4L, meanFunction = mean,
                                   includeAverageBaseQuality = FALSE, trace = FALSE, uniquelyAssigned = TRUE) {
    list(motifLength = motifLength, Type = Type, simpleReturn = simpleReturn, includeLUS = includeLUS, numberOfThreads = numberOfThreads, meanFunction = meanFunction,
         includeAverageBaseQuality = includeAverageBaseQuality, trace = trace, uniquelyAssigned = uniquelyAssigned)
}

.extractedReadsList.stringCoverage <- function(extractedReadsListObject, control = stringCoverage.control()) {
    if (control$uniquelyAssigned) {
        extractedReads <- extractedReadsListObject$identifiedMarkersSequencesUniquelyAssigned
    }
    else {
        warning("Only uniquely assigned sequences extract list should be used")
        extractedReads <- extractedReadsListObject$identifiedMarkers
    }

    if (length(control$motifLength) != length(extractedReads)) {
        if (length(control$motifLength) == 1L) {
            motifLengths <- rep(control$motifLength, length(extractedReads))
        }
        else {
            stop("'motifLength' must have length 1 or the same as 'extractedReads'")
        }

    }
    else {
        motifLengths = control$motifLength
    }

    if (length(control$Type) != length(extractedReads)) {
        if (length(control$Type) == 1L) {
            Types <- rep(control$Type, length(extractedReads))
        }
        else {
            stop("'Type' must have length 1 or the same as 'extractedReads'")
        }

    }
    else {
        Types = control$Type
    }

    alleles <- list()
    nullAlleles <- c()
    j = 1
    for (i in seq_along(extractedReads)) {
        if (dim(extractedReads[[i]]$trimmed)[1] == 0) {
            nullAlleles <- c(nullAlleles, i)
            next
        }

        if (control$trace)
            cat(i, "/", length(extractedReads), ":: Marker:", as.character(matchedFlanks$name), "\n")

        matchedFlanks <- extractedReads[[i]]
        marker = matchedFlanks$name

        stringCoverageQuality <- matchedFlanks$trimmed %>% group_by(ForwardFlank, Region, ReverseFlank) %>%
            summarise(Coverage = n()) %>% ungroup() %>%
            mutate(Marker = marker, MotifLength = motifLengths[i], Type = Types[i], Allele = nchar(Region) / MotifLength) %>%
            select(Marker, Allele, Type, MotifLength, ForwardFlank, Region, ReverseFlank, Coverage)

        if (control$simpleReturn) {
            stringCoverageQuality <- stringCoverageQuality %>% group_by(Marker, Allele, Type, MotifLength, Region) %>%
                summarise(Coverage = sum(Coverage)) %>% ungroup()
        }

        if (control$includeLUS) {
            validLUS = stringCoverageQuality$Allele >= 1
            stringCoverageQuality$LUS <- sapply(seq_along(stringCoverageQuality$Region), function(ss) ifelse(validLUS[ss], LUS(stringCoverageQuality$Region[ss], motifLength = stringCoverageQuality$MotifLength[ss], returnType = "string"), NA))
        }

        if (control$includeAverageBaseQuality) {
            warning("No longer a supported feature. The function will run with 'control$includeAverageBaseQuality' set as 'FALSE'.")
        }

        alleles[[j]] <- stringCoverageQuality
        j = j + 1
    }

    names(alleles) <- names(extractedReads)[!(1:length(extractedReads) %in% nullAlleles)]
    class(alleles) <- "stringCoverageList"
    return(alleles)
}

#' Get string coverage STR identified objects.
#'
#' \code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
#'
#' @param extractedReadsListObject an extractedReadsList-object, created using the \link{identifySTRRegions}-function.
#' @param control an \link{stringCoverage.control}-object.
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
setGeneric("stringCoverage", signature = "extractedReadsListObject",
           function(extractedReadsListObject, control = stringCoverage.control())
               standardGeneric("stringCoverage")
)

setMethod("stringCoverage", "extractedReadsList",
           function(extractedReadsListObject, control = stringCoverage.control())
               .extractedReadsList.stringCoverage(extractedReadsListObject, control)
)

setMethod("stringCoverage", "extractedReadsListReverseComplement",
          function(extractedReadsListObject, control = stringCoverage.control())
              .extractedReadsList.stringCoverage(extractedReadsListObject, control)
)

setMethod("stringCoverage", "extractedReadsListCombined",
          function(extractedReadsListObject, control = stringCoverage.control())
              .extractedReadsList.stringCoverage(extractedReadsListObject, control)
)

setMethod("stringCoverage", "extractedReadsListNonCombined",
          function(extractedReadsListObject, control = stringCoverage.control())
              stop("'stringCoverage' not implemented for 'extractedReadsListNReveseComplementList'. Use lapply on the two elements on the list.")
)

setClass("stringCoverageList")

.stringCoverageList.NoiseGenotype <- function(stringCoverageListObject, colBelief = "Coverage",
                                              thresholdSignal = 0, thresholdHeterozygosity = 0, thresholdAbsoluteLowerLimit = 1,
                                              trueGenotype = NULL, identified = "genotype") {
    if (length(thresholdSignal) == 1L) {
        if(thresholdSignal < 1 & thresholdSignal > 0) {
            thresholdSignal <- unlist(lapply(stringCoverageListObject, function(s) thresholdSignal*max(s[, colBelief])))
            thresholdSignal <- sapply(seq_along(thresholdSignal), function(s) max(thresholdSignal[s], thresholdAbsoluteLowerLimit))
        }
        else {
            thresholdSignal <- rep(max(thresholdSignal, thresholdAbsoluteLowerLimit), length(stringCoverageListObject))
        }
    }

    if (length(thresholdHeterozygosity) == 1L) {
        thresholdHeterozygosity <- rep(thresholdHeterozygosity, length(stringCoverageListObject))
    }

    if (length(thresholdSignal) != length(stringCoverageListObject)) {
        stop("'stringCoverageListObject' and 'thresholdSignal' must have the same length.")
    }

    if (length(thresholdHeterozygosity) != length(stringCoverageListObject)) {
        stop("'stringCoverageListObject' and 'thresholdHeterozygosity' must have the same length.")
    }

    res <- vector("list", length(stringCoverageListObject))
    for (i in seq_along(stringCoverageListObject)) {
        stringCoverage_i <- stringCoverageListObject[[i]]
        if (is.null(trueGenotype)) {
            belief <- unname(stringCoverage_i[, colBelief] %>% as_vector())
            beliefMax <- max(belief)
            beliefKeepers <- which(belief > thresholdSignal[i] & belief > thresholdHeterozygosity[i]*beliefMax)
        }
        else {
            beliefKeepers <- which(stringCoverage_i$Region %in% trueGenotype[[i]])
        }
        res[[i]] <- stringCoverage_i[beliefKeepers, ] %>% mutate(Indices = beliefKeepers)
    }

    names(res) <- names(stringCoverageListObject)
    class(res) <- if(tolower(identified) == "genotype") "genotypeIdentifiedList" else if(tolower(identified) == "noise") "noiseIdentifiedList"
    return(res)
}

#' Assigns genotype.
#'
#' \code{getGenotype} takes an stringCoverageList-object, assumes the sample is a reference file and assings a genotype, based on a heterozygote threshold, for every marker in the provided list.
#'
#' @param stringCoverageListObject an stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param colBelief the name of the coloumn used for identification.
#' @param thresholdSignal threshold applied to the signal (generally the coverage) of every string.
#' @param thresholdHeterozygosity threshold used to determine whether a marker is hetero- or homozygous.
#'
#' @return Returns a list, with an element for every marker in stringCoverageList-object, each element contains the genotype for a given marker.
setGeneric("getGenotype", signature = "stringCoverageListObject",
           function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0, thresholdHeterozygosity = 0.35, thresholdAbsoluteLowerLimit = 1)
               standardGeneric("getGenotype")
)

setMethod("getGenotype", "stringCoverageList",
          function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0, thresholdHeterozygosity = 0.35, thresholdAbsoluteLowerLimit = 1)
              .stringCoverageList.NoiseGenotype(stringCoverageListObject, colBelief, thresholdSignal, thresholdHeterozygosity,
                                                thresholdAbsoluteLowerLimit, NULL, "genotype")
)

#' Idenfities the noise.
#'
#' \code{identifyNoise} takes an stringCoverageList-object and identifies the noise based on a signal threshold for every marker in the provided list.
#'
#' @param stringCoverageListObject an stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param colBelief the name of the coloumn used for identification.
#' @param thresholdSignal threshold applied to the signal (generally the coverage) of every string.
#'
#' @return Returns a list, with an element for every marker in stringCoverageList-object, each element contains the genotype for a given marker.
setGeneric("identifyNoise", signature = "stringCoverageListObject",
           function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0.01, thresholdAbsoluteLowerLimit = 1)
               standardGeneric("identifyNoise")
)

setMethod("identifyNoise", "stringCoverageList",
          function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0.01, thresholdAbsoluteLowerLimit = 1)
              .stringCoverageList.NoiseGenotype(stringCoverageListObject, colBelief, thresholdSignal, 0, thresholdAbsoluteLowerLimit, NULL, "noise")
)

setClass("genotypeIdentifiedList")
setClass("noiseIdentifiedList")

.noiseGenotypeIdentified.stringCoverageList.merge <- function(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "genotype") {
    stringCoverageListObjectMerged <- vector("list", length(stringCoverageListObject))
    indValue <- if(tolower(identified) == "genotype") TRUE else if(tolower(identified) == "noise") FALSE
    indCol <- if(tolower(identified) == "genotype") "AlleleCalled" else if(tolower(identified) == "noise") "Noise"

    for(i in seq_along(stringCoverageListObject)) {
        stringCoverageListObjectMerged[[i]] <- stringCoverageListObject[[i]] %>% mutate(tempName = !indValue, FLAGMoreThanTwoAlleles = FALSE)

        if (!is.null(noiseGenotypeIdentifiedListObject[[i]]) && nrow(noiseGenotypeIdentifiedListObject[[i]]) > 0L) {
            stringCoverageListObjectMerged[[i]]$tempName[noiseGenotypeIdentifiedListObject[[i]]$Indices] <- indValue
        }

        if (nrow(noiseGenotypeIdentifiedListObject[[i]]) > 2L) {
            stringCoverageListObjectMerged[[i]]$FLAGMoreThanTwoAlleles <- TRUE
        }

        names(stringCoverageListObjectMerged[[i]]) <- gsub("tempName", indCol, names(stringCoverageListObjectMerged[[i]]))
    }

    names(stringCoverageListObjectMerged) <- names(stringCoverageListObject)
    class(stringCoverageListObjectMerged) <- if(tolower(identified) == "genotype") "stringCoverageGenotypeList" else if(tolower(identified) == "noise") "stringCoverageNoiseList"
    return(stringCoverageListObjectMerged)
}

#' Merge genotypeIdentifiedList and stringCoverageList.
#'
#' \code{mergeGenotypeStringCoverage} merges genotypeIdentifiedList-objects and stringCoverageList-objects.
#'
#' @param stringCoverageListObject a stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param noiseGenotypeIdentifiedListObject a noiseGenotypeIdentifiedList-object, created using the \link{getGenotype}-function.
#'
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
setGeneric("mergeGenotypeStringCoverage", signature = "noiseGenotypeIdentifiedListObject",
           function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
               standardGeneric("mergeGenotypeStringCoverage")
)

setMethod("mergeGenotypeStringCoverage", "genotypeIdentifiedList",
          function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
              .noiseGenotypeIdentified.stringCoverageList.merge(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "genotype")
)

#' Merge noiseIdentifiedList and stringCoverageList.
#'
#' \code{mergeNoiseStringCoverage} merges noiseIdentifiedList-objects and stringCoverageList-objects.
#'
#' @param stringCoverageListObject a stringCoverageList-object, created using the \link{stringCoverage}-function.
#' @param noiseGenotypeIdentifiedListObject a noiseGenotypeIdentifiedList-object, created using the \link{identifyNoise}-function.
#'
#' @return Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
setGeneric("mergeNoiseStringCoverage", signature = "noiseGenotypeIdentifiedListObject",
           function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
               standardGeneric("mergeNoiseStringCoverage")
)

setMethod("mergeNoiseStringCoverage", "noiseIdentifiedList",
          function(stringCoverageListObject, noiseGenotypeIdentifiedListObject)
              .noiseGenotypeIdentified.stringCoverageList.merge(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "noise")
)

setClass("stringCoverageGenotypeList")
setClass("stringCoverageNoiseList")
