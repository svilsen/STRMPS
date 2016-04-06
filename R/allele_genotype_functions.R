.findStringCoverage <- function(matchedFlanks, matchedFlanksSplit, matchedFlanksReverseComplement = NULL, matchedFlanksQuality = NULL,
                                motifLength, meanFunction, includeLUS, numberOfThreads) {
    res <- mclapply(seq_along(matchedFlanksSplit), function(j) {
        strings_j <- split(seq_along(as.character(matchedFlanks$trimmedIncludingFlanks[matchedFlanksSplit[[j]]])),
                       as.character(matchedFlanks$trimmedIncludingFlanks[matchedFlanksSplit[[j]]]))

        strRegions_j <- split(seq_along(as.character(matchedFlanks$trimmed[matchedFlanksSplit[[j]]])),
                           as.character(matchedFlanks$trimmed[matchedFlanksSplit[[j]]]))

        stringOfSTRRegion <- unname(unlist(lapply(strings_j, function(s) {
            whichString <- unlist(lapply(strRegions_j, function(ss) sum(ss %in% s) > 0))
            names(whichString)[whichString]
        })))

        stringCoverage_j <- structure(data.frame(as.numeric(names(matchedFlanksSplit[j])), names(strings_j),
                                                 stringOfSTRRegion, unname(unlist(lapply(strings_j, length))),
                                                 stringsAsFactors = FALSE),
                           .Names = c("Allele", "String", "STRRegion", "Coverage"))

        if (!is.null(matchedFlanksReverseComplement)) {
            stringCoverage_j$RCPercentage <- unname(unlist(lapply(strings_j, function(x) sum(matchedFlanksReverseComplement[[j]][x, 2]) / (length(matchedFlanksReverseComplement[[j]][x, 2])))))
        }

        if (!is.null(matchedFlanksQuality)) {
            quality_j <- as.matrix(PhredQuality(matchedFlanksQuality[matchedFlanksSplit[[j]]]))
            PhredQuality_j <- lapply(strings_j, function(s) {
                quality_s <- if(is.null(dim(quality_j))) quality_j else quality_j[s, ]
                qualityAvg_s <- if(is.null(dim(quality_s))) quality_s else apply(quality_s, 2, meanFunction)
                probAvg_s <- 10^(qualityAvg_s/(-10))
                phredQuality <- PhredQuality(probAvg_s[!is.na(probAvg_s)])
                return(phredQuality)
            })
            stringCoverage_j$AverageStringPhredQuality <- unname(unlist(lapply(PhredQuality_j, as.character)))
        }

        if (includeLUS) {
            stringCoverage_j$LUS <- unlist(lapply(stringCoverage_j$STRRegion, function(s) LUS(s, motifLength = 4, returnType = "string")))
        }

        return(stringCoverage_j)
    }, mc.cores = numberOfThreads)

    return(res)
}

stringCoverage.control <- function(motifLength = 4, includeLUS = TRUE, numberOfThreads = 4L, meanFunction = mean,
                                   includeAverageBaseQuality = TRUE, trace = FALSE, uniquelyAssigned = TRUE) {
    list(motifLength = motifLength, includeLUS = includeLUS, numberOfThreads = numberOfThreads, meanFunction = meanFunction,
         includeAverageBaseQuality = includeAverageBaseQuality, trace = trace, uniquelyAssigned = uniquelyAssigned)
}

.extractedReadsList.stringCoverage <- function(extractedReadsListObject, control = stringCoverage.control()) {
    if (control$uniquelyAssigned) {
        extractedReads <- extractedReadsListObject$identifiedMarkersSequencesUniquelyAssigned
    }
    else {
        extractedReads <- extractedReadsListObject$identifiedMarkers
    }

    if (length(control$motifLength) == 1L) {
        motifLengths <- rep(control$motifLength, length(extractedReads))
    }
    else if (length(motifLengths) == length(extractedReads)) {
        motifLengths = control$motifLength
    }
    else {
        stop("'motifLenght' must have length 1 or the same as 'extractedReads'")
    }

    alleles <- list()
    nullAlleles <- c()
    j = 1
    for (i in seq_along(extractedReads)) {
        if (is.null(extractedReadsListObject$identifiedMarkers[[i]])) {
            nullAlleles <- c(nullAlleles, i)
            next
        }

        matchedFlanks <- extractedReads[[i]]

        if (control$trace)
            cat(i, "/", length(extractedReads), ":: Marker:", as.character(matchedFlanks$name), "\n")

        widths <- nchar(matchedFlanks$trimmed) / motifLengths[i]
        matchedFlanksSplit <- split(seq_along(widths), widths)

        if (class(extractedReadsListObject) == "extractedReadsListCombined") {
            matchedFlanksReverseComplement <- matchedFlanks$ReverseComplement
            matchedFlanksReverseComplementSplit <- lapply(matchedFlanksSplit, function(x) matchedFlanksReverseComplement[x, ])
        } else {
            matchedFlanksReverseComplementSplit <- NULL
        }

        if (control$includeAverageBaseQuality) {
            matchedFlanksQuality <- matchedFlanks$trimmedQuality
        } else {
            matchedFlanksQuality <- NULL
        }

        stringCoverageQuality <- .findStringCoverage(matchedFlanks, matchedFlanksSplit, matchedFlanksReverseComplement = matchedFlanksReverseComplementSplit,
                                                matchedFlanksQuality = matchedFlanksQuality, motifLength = motifLengths[i], meanFunction = control$meanFunction,
                                                includeLUS = control$includeLUS, numberOfThreads = control$numberOfThreads)

        alleles[[j]] <- do.call(rbind, stringCoverageQuality)
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
#'
#' @export
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
                                              thresholdSignal = 0, thresholdHeterozygosity = 0,
                                              trueGenotype = NULL, identified = "genotype") {
    if (length(thresholdSignal) == 1L) {
        if(thresholdSignal < 1 & thresholdSignal > 0) {
            thresholdSignal <- unlist(lapply(stringCoverageListObject, function(s) thresholdSignal*max(s[, colBelief])))
        }

        thresholdSignal <- rep(thresholdSignal, length(stringCoverageListObject))
    }

    if (length(thresholdSignal) != length(stringCoverageListObject)) {
        stop("alleles and thresholdSignal must have the same length.")
    }

    res <- vector("list", length(stringCoverageListObject))
    for (i in seq_along(stringCoverageListObject)) {
        stringCoverage_i <- stringCoverageListObject[[i]]
        if (is.null(trueGenotype)) {
            belief <- stringCoverage_i[, colBelief]
            beliefMax <- max(belief)
            beliefKeepers <- which(belief > thresholdSignal[i] & belief > thresholdHeterozygosity*beliefMax)
        }
        else {
            beliefKeepers <- which(stringCoverage_i$String %in% trueGenotype[[i]])
        }
        res[[i]] <- cbind(stringCoverage_i[beliefKeepers, ], Indices = beliefKeepers)
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
#'
#' @export
setGeneric("getGenotype", signature = "stringCoverageListObject",
           function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0, thresholdHeterozygosity = 0.35)
               standardGeneric("getGenotype")
)

setMethod("getGenotype", "stringCoverageList",
          function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0, thresholdHeterozygosity = 0.35)
              .stringCoverageList.NoiseGenotype(stringCoverageListObject, colBelief, thresholdSignal, thresholdHeterozygosity, identified = "genotype")
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
#'
#' @export
setGeneric("identifyNoise", signature = "stringCoverageListObject",
           function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0.01)
               standardGeneric("identifyNoise")
)

setMethod("identifyNoise", "stringCoverageList",
          function(stringCoverageListObject, colBelief = "Coverage", thresholdSignal = 0.01)
              .stringCoverageList.NoiseGenotype(stringCoverageListObject, colBelief, thresholdSignal, thresholdHeterozygosity = 0, identified = "noise")
)

setClass("genotypeIdentifiedList")
setClass("noiseIdentifiedList")

.noiseGenotypeIdentified.stringCoverageList.merge <- function(stringCoverageListObject, noiseGenotypeIdentifiedListObject, identified = "genotype") {
    stringCoverageListObjectMerged <- vector("list", length(stringCoverageListObject))
    indValue <- if(tolower(identified) == "genotype") TRUE else if(tolower(identified) == "noise") FALSE
    indCol <- if(tolower(identified) == "genotype") "AlleleCalled" else if(tolower(identified) == "noise") "Noise"

    for(i in seq_along(stringCoverageListObject)) {
        stringCoverageListObjectMerged[[i]] <- cbind(stringCoverageListObject[[i]], tempName = !indValue)

        if (!is.null(noiseGenotypeIdentifiedListObject[[i]]) && nrow(noiseGenotypeIdentifiedListObject[[i]]) > 0L) {
            stringCoverageListObjectMerged[[i]]$tempName[noiseGenotypeIdentifiedListObject[[i]]$Indices] <- indValue
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
#'
#' @export
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
#'
#' @export
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
