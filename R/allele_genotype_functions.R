.findStringCoverage <- function(matchedFlanks, matchedFlanksSplit, matchedFlanksReverseComplement = NULL, matchedFlanksQuality = NULL,
                                motifLength, Type, flankingRegionLength, meanFunction, includeLUS, numberOfThreads) {
    res <- mclapply(seq_along(matchedFlanksSplit), function(j) {
        strings_j <- split(seq_along(as.character(matchedFlanks$trimmedIncludingFlanks[matchedFlanksSplit[[j]]])),
                       as.character(matchedFlanks$trimmedIncludingFlanks[matchedFlanksSplit[[j]]]))

        strRegions_j <- split(seq_along(as.character(matchedFlanks$trimmed[matchedFlanksSplit[[j]]])),
                           as.character(matchedFlanks$trimmed[matchedFlanksSplit[[j]]]))


        regionPointer_j <- unlist(lapply(strings_j, function(s) unname(which(sapply(strRegions_j, function(r) all(s %in% r))))))

        stringOfSTRRegion <- lapply(strings_j, function(s) {
            whichString <- unlist(lapply(strings_j, function(ss) sum(s %in% ss) > 0))
            string <- names(strings_j)[whichString]
            region <- names(strRegions_j)[regionPointer_j[whichString]]

            regionFlankingTest <- substr(string, flankingRegionLength[1] + 1, nchar(string) - flankingRegionLength[2])
            forwardFlank = NA
            reverseFlank = NA
            if (region == regionFlankingTest) {
                forwardFlank <- substr(string, 1, flankingRegionLength[1])
                reverseFlank <- substr(string, nchar(string) - flankingRegionLength[2] + 1, nchar(string))
            }
            else {
                flanks <- strsplit(string, region)
                if (length(flanks[[1]]) == 2) {
                    forwardFlank <- flanks[[1]][1]
                    reverseFlank <- flanks[[1]][2]
                }
            }

            string_coverage <- sapply(strings_j[whichString], length)
            tibble(Allele = nchar(region) / motifLength, Type = Type, MotifLength = motifLength,
                   ForwardFlank = forwardFlank, Region = region, ReverseFlank = reverseFlank, Coverage = unname(string_coverage))
        })

        stringCoverage_j <- bind_rows(stringOfSTRRegion)

        if (!is.null(matchedFlanksReverseComplement[[j]])) {
            stringCoverage_j$RCPercentage <- unname(unlist(lapply(strings_j, function(x) sum(matchedFlanksReverseComplement[[j]][x, 2]) / (length(matchedFlanksReverseComplement[[j]][x, 2])))))
        }
        else {
            stringCoverage_j$RCPercentage <- NA
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
            stringCoverage_j$RegionAveragePhredQuality <- unname(unlist(lapply(PhredQuality_j, as.character)))
        }

        if (includeLUS) {
            stringCoverage_j$LUS <- NA
            if (stringCoverage_j$Allele >= 1) {
                stringCoverage_j$LUS <- unlist(lapply(stringCoverage_j$Region, function(s) LUS(s, motifLength = motifLength, returnType = "string")))
            }
        }

        return(stringCoverage_j)
    }, mc.cores = numberOfThreads)

    return(res)
}

stringCoverage.control <- function(motifLength = 4, Type = "AUTOSOMAL", flankingRegionLength = 12, includeLUS = TRUE, numberOfThreads = 4L, meanFunction = mean,
                                   includeAverageBaseQuality = FALSE, trace = FALSE, uniquelyAssigned = TRUE) {
    list(motifLength = motifLength, Type = Type, flankingRegionLength = flankingRegionLength, includeLUS = includeLUS, numberOfThreads = numberOfThreads, meanFunction = meanFunction,
         includeAverageBaseQuality = includeAverageBaseQuality, trace = trace, uniquelyAssigned = uniquelyAssigned)
}

.extractedReadsList.stringCoverage <- function(extractedReadsListObject, control = stringCoverage.control()) {
    if (control$uniquelyAssigned) {
        extractedReads <- extractedReadsListObject$identifiedMarkersSequencesUniquelyAssigned
    }
    else {
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

    if (is.matrix(control$flankingRegionLength) & (dim(control$flankingRegionLength)[1] == length(extractedReads))) {
        flankingRegionLengths = control$flankingRegionLength
    }
    else if (is.numeric(control$flankingRegionLength) & (length(control$flankingRegionLength) == 1)) {
        flankingRegionLengths = matrix(control$flankingRegionLength, nrow = length(extractedReads), ncol = 2)
    }
    else {
        stop("'flankingRegionLengths' must be a vector of length 1 or a matrix with 2 columns and number of rows equal to the length of 'extractedReads'")
    }


    alleles <- list()
    nullAlleles <- c()
    j = 1
    for (i in seq_along(extractedReads)) {
        if (length(extractedReads[[i]]$trimmed) == 0) {
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
        }
        else {
            matchedFlanksReverseComplementSplit <- NULL
        }

        if (control$includeAverageBaseQuality) {
            matchedFlanksQuality <- matchedFlanks$trimmedQuality
        }
        else {
            matchedFlanksQuality <- NULL
        }

        stringCoverageQuality <- .findStringCoverage(matchedFlanks, matchedFlanksSplit, matchedFlanksReverseComplement = matchedFlanksReverseComplementSplit,
                                                matchedFlanksQuality = matchedFlanksQuality, motifLength = motifLengths[i], Type = Types[i],
                                                flankingRegionLength = flankingRegionLengths[i, ], meanFunction = control$meanFunction, includeLUS = control$includeLUS, numberOfThreads = control$numberOfThreads)

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

    if (length(thresholdSignal) != length(stringCoverageListObject)) {
        stop("alleles and thresholdSignal must have the same length.")
    }

    res <- vector("list", length(stringCoverageListObject))
    for (i in seq_along(stringCoverageListObject)) {
        stringCoverage_i <- stringCoverageListObject[[i]]
        if (is.null(trueGenotype)) {
            belief <- unname(stringCoverage_i[, colBelief] %>% as_vector())
            beliefMax <- max(belief)
            beliefKeepers <- which(belief > thresholdSignal[i] & belief > thresholdHeterozygosity*beliefMax)
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
#'
#' @export
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
#'
#' @export
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
