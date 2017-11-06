## .mclapplyvmatchPattern and .vmatchMultiplePatternsAlternate are used for diagnostic purposes only
.mclapplyvmatchPattern <- function(flanks, seqs, max.mismatch, numberOfThreads, limitSeqRange = NULL) {
    if (is.null(limitSeqRange)) {
        resList <- mclapply(seq_along(flanks), function(f) vmatchPattern(flanks[[f]], seqs, max.mismatch = max.mismatch), mc.cores = numberOfThreads)
    }
    else {
        resList <- mclapply(seq_along(flanks), function(f) vmatchPattern(flanks[[f]], seqs[limitSeqRange[[f]]], max.mismatch = max.mismatch), mc.cores = numberOfThreads)
    }
    return(resList)
}

.getCols <- function(colNames) {
    structure(lapply(c("marker", "forwardflank", "reverseflank", "forwardshift", "reverseshift"), function(n) grep(n, tolower(colNames))), .Names = c("markerCol", "forwardCol", "reverseCol", "forwardShiftCol", "reverseShiftCol"))
}

.identifyFlankingRegions <- function(seqs, flankingRegions, colList, numberOfMutation = 1, numberOfThreads = 4, removeEmptyMarkers = TRUE) {
    forwardFlank <- DNAStringSet(unname(unlist(flankingRegions[, colList$forwardCol])))
    reverseFlank <- DNAStringSet(unname(unlist(flankingRegions[, colList$reverseCol])))

    # Forward match
    identifiedForward <- STRMPS:::.mclapplyvmatchPattern(forwardFlank, seqs, numberOfMutation, numberOfThreads, limitSeqRange = NULL)
    forwardMatch <- lapply(1:nrow(flankingRegions), function(i) which(elementNROWS(identifiedForward[[i]]) >= 1L))

    # Reverse match
    identifiedReverse <- STRMPS:::.mclapplyvmatchPattern(reverseFlank, seqs, numberOfMutation, numberOfThreads, limitSeqRange = forwardMatch)
    reverseMatch <- lapply(1:nrow(flankingRegions), function(i) which(elementNROWS(identifiedReverse[[i]]) >= 1L))

    matchedSeq <- lapply(seq_along(reverseMatch), function(i) forwardMatch[[i]][reverseMatch[[i]]])

    if(removeEmptyMarkers) {
        nonEmptyEntries <- which(sapply(matchedSeq, length) != 0)
        rList <- list(markers = unname(as_vector(flankingRegions[nonEmptyEntries, colList$markerCol])),
                      identifiedForward = identifiedForward[nonEmptyEntries], identifiedReverse = identifiedReverse[nonEmptyEntries],
                      matchedSeq = matchedSeq[nonEmptyEntries], reverseMatchedSeq = reverseMatch[nonEmptyEntries])
    }
    else {
        rList <- list(markers = unname(unlist((flankingRegions[, colList$markerCol]))), identifiedForward = identifiedForward,
                      identifiedReverse = identifiedReverse, matchedSeq = matchedSeq, reverseMatchedSeq = reverseMatch)
    }
    return(rList)
}

.extractAndTrimMarkerIdentifiedReads.mclapply <- function(seqs, qual, flankingRegions, colList, numberOfMutation, removeEmptyMarkers,
                                                          flankSizes, numberOfThreads = 4, reversed = FALSE, reverseComplementRun = FALSE,
                                                          trimming = "flankshift") {
    identifiedFlanksObj <- STRMPS:::.identifyFlankingRegions(seqs, flankingRegions, colList = colList, numberOfMutation = numberOfMutation, numberOfThreads = numberOfThreads, removeEmptyMarkers = removeEmptyMarkers)

    identifiedMarkers <- structure(mclapply(seq_along(identifiedFlanksObj$matchedSeq), function(i) {
        marker <- identifiedFlanksObj$markers[i]
        identifiedForward <- identifiedFlanksObj$identifiedForward[[i]]
        identifiedReverse <- identifiedFlanksObj$identifiedReverse[[i]]
        matchedSeq <- identifiedFlanksObj$matchedSeq[[i]]
        reverseMatchedSeq <- identifiedFlanksObj$reverseMatchedSeq[[i]]

        # F: Take the first in the IRanges
        endForward <- unlist(lapply(endIndex(identifiedForward)[matchedSeq], function(v) v[1L]))
        startForward <- endForward - (flankSizes[i, 1] - 1)

        # R: Take the last in the IRanges
        startReverse <- unlist(lapply(startIndex(identifiedReverse)[reverseMatchedSeq], function(v) v[length(v)]))
        endReverse <- startReverse + (flankSizes[i, 2] - 1)

        # Removes reads where forward is NOT observed before reverse
        keepSeq <- which(startReverse > (endForward + 1) & endReverse <= nchar(seqs[matchedSeq]) & startForward > 0)
        matchedSeq <- matchedSeq[keepSeq]

        startForward <- startForward[keepSeq]
        endForward <- endForward[keepSeq]

        startReverse <- startReverse[keepSeq]
        endReverse <- endReverse[keepSeq]

        if (length(keepSeq) == 0L) {
            rList <- NULL
        }
        else {
            flankShift = unlist(flankingRegions %>% filter(Marker == marker) %>% select(contains("Shift")))
            trimmedIncludingFlanks <- subseq(seqs[matchedSeq], start = startForward, end = endReverse)
            trimmedQualityIncludingFlanks <- subseq(quality(qual[matchedSeq]), start = startForward, end = endReverse)

            if (tolower(trimming) == "flankshift") {
                if (!(reversed && reverseComplementRun)) {
                    forwardFlankingRegion <- subseq(seqs[matchedSeq], start = startForward, end = endForward + flankShift[1] + 1)
                    reverseFlankingRegion <- subseq(seqs[matchedSeq], start = startReverse - flankShift[2] - 1, end = endReverse)

                    trimmed <- subseq(seqs[matchedSeq], start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1)

                    forwardFlankingRegionQuality <- subseq(quality(qual[matchedSeq]), start = startForward, end = endForward + flankShift[1] + 1)
                    reverseFlankingRegionQuality <- subseq(quality(qual[matchedSeq]), start = startReverse - flankShift[2] - 1, end = endReverse)

                    trimmedQuality <- subseq(quality(qual[matchedSeq]), start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1)
                }
                else {
                    forwardFlankingRegion <- reverseComplement(subseq(seqs[matchedSeq], start = startReverse - flankShift[2] - 1, end = endReverse))
                    reverseFlankingRegion <- reverseComplement(subseq(seqs[matchedSeq], start = startForward, end = endForward + flankShift[1] + 1))

                    trimmed <- reverseComplement(subseq(seqs[matchedSeq], start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1))

                    forwardFlankingRegionQuality <- reverse(subseq(quality(qual[matchedSeq]), start = startReverse - flankShift[2] - 1, end = endReverse))
                    reverseFlankingRegionQuality <- reverse(subseq(quality(qual[matchedSeq]), start = startForward, end = endForward + flankShift[1] + 1))

                    trimmedQuality <- reverse(subseq(quality(qual[matchedSeq]), start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1))
                }
            } else if (tolower(trimming) == "directflankingregion") {
                stop("Trimming by 'directflankingregion' is not implemented... yet.")
            }

            trimmedTibble <- tibble(ForwardFlank = as.character(forwardFlankingRegion), Region = as.character(trimmed), ReverseFlank = as.character(reverseFlankingRegion))
            trimmedQualityTibble <- tibble(ForwardFlank = as.character(forwardFlankingRegionQuality), Region = as.character(trimmedQuality), ReverseFlank = as.character(reverseFlankingRegionQuality))
            rList <- list(name = marker, matchedSeq = matchedSeq,
                          startForward = startForward, endForward = endForward,
                          startReverse = startReverse, endReverse = endReverse,
                          trimmedIncludingFlanks = trimmedIncludingFlanks, trimmedQualityIncludingFlanks = trimmedQualityIncludingFlanks,
                          trimmed = trimmedTibble, trimmedQuality = trimmedQualityTibble)
        }

        return(rList)
    }, mc.cores = numberOfThreads), .Names = as.character(identifiedFlanksObj$markers))

    # Remove sequences occuring at several loci
    identifiedSeqTable <- sort(table(unlist(lapply(identifiedMarkers, function(r) r$matchedSeq))), decreasing = TRUE)
    identifiedSeqNotUnique <- as.integer(names(identifiedSeqTable[identifiedSeqTable > 1L]))
    identifiedMarkersUniqueSequences <- lapply(identifiedMarkers, function(r) {
        keepSeq <- which(!(r$matchedSeq %in% identifiedSeqNotUnique))
        return(list(name = r$name, matchedSeq = r$matchedSeq[keepSeq], startForward = r$startForward[keepSeq], endForward = r$endForward[keepSeq],
                    startReverse = r$startReverse[keepSeq], endReverse = r$endReverse[keepSeq],
                    trimmedIncludingFlanks = r$trimmedIncludingFlanks[keepSeq], trimmedQualityIncludingFlanks = r$trimmedQualityIncludingFlanks[keepSeq],
                    trimmed = r$trimmed[keepSeq, ], trimmedQuality = r$trimmedQuality[keepSeq, ]))
    })

    extractedReadsList <- list(n_reads = sum(sapply(identifiedMarkersUniqueSequences, function(vv) length(vv[["matchedSeq"]]))),
                               reverseComplement = reverseComplementRun, identifiedMarkers = identifiedMarkers,
                               identifiedMarkersSequencesUniquelyAssigned = identifiedMarkersUniqueSequences)
    class(extractedReadsList) <- "extractedReadsList"
    return(extractedReadsList)
}


.extractAndTrimMarkerIdentifiedReads.seqAn <- function(seqs, qual, flankingRegions, colList, numberOfMutation = 1, removeEmptyMarkers = TRUE,
                                                       flankSizes, numberOfThreads = 4, reversed = FALSE, reverseComplementRun = FALSE,
                                                       trimming = "flankShift") {
    markers <- unname(unlist(flankingRegions[, colList$markerCol]))
    forwardFlank <- unname(unlist(flankingRegions[, colList$forwardCol]))
    reverseFlank <- unname(unlist(flankingRegions[, colList$reverseCol]))

    seqsCharacter <- as.character(seqs)
    maxNumberOfMismatch <- Biostrings:::normargMaxMismatch(numberOfMutation)

    matched <- do.call("rbind", mclapply(seq_along(seqs), function(ss) vmatchMultiPatternSeqAn(forwardFlank, reverseFlank, seqsCharacter[ss], ss, maxNumberOfMismatch), mc.cores = numberOfThreads))
    colnames(matched) <- c("Sequence", "Marker", "ForwardEnd", "ReverseEnd")
    matched <- as_tibble(matched) %>% filter(Marker > 0)

    identifiedMarkers <- structure(mclapply(seq_along(markers), function(m) {
        matched_m <- matched %>% filter(Marker == m, ForwardEnd < (ReverseEnd - (flankSizes[m, 2] - 1)))

        #cat("Marker:", markers[m], ":: ID'ed:", length((matched %>% filter(Marker == m))$Sequence), ":: Kept:", length(matched_m$Sequence), "\n")
        if (length(matched_m$Sequence) == 0L) {
            rList <- NULL
        }
        else {
            flankShift = unlist(flankingRegions %>% filter(Marker == markers[unique(matched_m$Marker)]) %>% select(contains("Shift")))

            endForward <- unname(unlist(matched_m$ForwardEnd))
            startForward <- endForward - (flankSizes[m, 1] - 1)

            endReverse <- unname(unlist(matched_m$ReverseEnd)) - 1
            startReverse <- endReverse - (flankSizes[m, 2] - 1)

            trimmedIncludingFlanks <- subseq(seqs[matched_m$Sequence], start = startForward, end = endReverse)
            trimmedQualityIncludingFlanks <- subseq(quality(qual[matched_m$Sequence]), start = startForward, end = endReverse)

            if (tolower(trimming) == "flankshift") {
                if (!(reversed && reverseComplementRun)) {
                    forwardFlankingRegion <- subseq(seqs[matched_m$Sequence], start = startForward, end = endForward + flankShift[1] + 1)
                    reverseFlankingRegion <- subseq(seqs[matched_m$Sequence], start = startReverse - flankShift[2] - 1, end = endReverse)

                    trimmed <- subseq(seqs[matched_m$Sequence], start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1)

                    forwardFlankingRegionQuality <- subseq(quality(qual[matched_m$Sequence]), start = startForward, end = endForward + flankShift[1] + 1)
                    reverseFlankingRegionQuality <- subseq(quality(qual[matched_m$Sequence]), start = startReverse - flankShift[2] - 1, end = endReverse)

                    trimmedQuality <- subseq(quality(qual[matched_m$Sequence]), start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1)
                } else {
                    forwardFlankingRegion <- reverseComplement(subseq(seqs[matched_m$Sequence], start = startReverse - flankShift[2] - 1, end = endReverse))
                    reverseFlankingRegion <- reverseComplement(subseq(seqs[matched_m$Sequence], start = startForward, end = endForward + flankShift[1] + 1))

                    trimmed <- reverseComplement(subseq(seqs[matched_m$Sequence], start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1))

                    forwardFlankingRegionQuality <- reverse(subseq(quality(qual[matched_m$Sequence]), start = startReverse - flankShift[2] - 1, end = endReverse))
                    reverseFlankingRegionQuality <- reverse(subseq(quality(qual[matched_m$Sequence]), start = startForward, end = endForward + flankShift[1] + 1))

                    trimmedQuality <- reverse(subseq(quality(qual[matched_m$Sequence]), start = endForward + flankShift[1] + 1, end = startReverse - flankShift[2] - 1))
                }

            }
            else if (tolower(trimming) == "directflankingregion") {
                stop("Trimming by 'directflankingregion' is not implemented... yet.")
            }

            trimmedTibble <- tibble(ForwardFlank = as.character(forwardFlankingRegion), Region = as.character(trimmed), ReverseFlank = as.character(reverseFlankingRegion))
            trimmedQualityTibble <- tibble(ForwardFlank = as.character(forwardFlankingRegionQuality), Region = as.character(trimmedQuality), ReverseFlank = as.character(reverseFlankingRegionQuality))
            rList <- list(name = markers[m], matchedSeq = matched_m$Sequence,
                          startForward = startForward, endForward = endForward,
                          startReverse = startReverse, endReverse = endReverse,
                          trimmedIncludingFlanks = trimmedIncludingFlanks, trimmedQualityIncludingFlanks = trimmedQualityIncludingFlanks,
                          trimmed = trimmedTibble, trimmedQuality = trimmedQualityTibble)
        }
    }, mc.cores = numberOfThreads), .Names = markers) #

    numberOfReadsID = sum(sapply(identifiedMarkers, function(vv) length(vv$matchedSeq)))
    extractedReadsList <- list(n_reads = numberOfReadsID, reverseComplement = reverseComplementRun, identifiedMarkers = NULL,
                               identifiedMarkersSequencesUniquelyAssigned = identifiedMarkers)
    class(extractedReadsList) <- "extractedReadsList"
    return(extractedReadsList)
}



.extractAndTrimMarkerIdentifiedReadsReverseComplement <- function(seqs, qual, flankingRegions, colList = NULL,
                                                                  numberOfMutation, removeEmptyMarkers = TRUE,
                                                                  numberOfThreads = numberOfThreads, reversed = TRUE,
                                                                  trimming = "flankshift", matchPatternMethod = "seqan") {
    flankingRegionsReverseComplement <- structure(data.frame(flankingRegions[, colList$markerCol],
                                                             unname(as.character(reverseComplement(DNAStringSet(as_vector(flankingRegions[, colList$reverseCol]))))),
                                                             unname(as.character(reverseComplement(DNAStringSet(as_vector(flankingRegions[, colList$forwardCol]))))),
                                                             flankingRegions[, colList$reverseShiftCol],
                                                             flankingRegions[, colList$forwardShiftCol],
                                                             stringsAsFactors = FALSE),
                                                  .Names = c("Marker", "ForwardRC", "ReverseRC", "ForwardShiftRC", "ReverseShiftRC"))

    extractAndTrimMarkerIdentifiedReads = switch(matchPatternMethod,
                                                 "mclapply" = STRMPS:::.extractAndTrimMarkerIdentifiedReads.mclapply,
                                                 "seqan" = STRMPS:::.extractAndTrimMarkerIdentifiedReads.seqAn)

    extractedTrimmedRC <- extractAndTrimMarkerIdentifiedReads(seqs, qual, flankingRegionsReverseComplement,
                                                              colList = list(markerCol = 1, forwardCol = 2, reverseCol = 3, forwardShiftCol = 4, reverseShiftCol = 5),
                                                              numberOfMutation = numberOfMutation, removeEmptyMarkers = removeEmptyMarkers,
                                                              flankSizes = apply(flankingRegionsReverseComplement[, -1], 2, nchar),
                                                              numberOfThreads = numberOfThreads,
                                                              reversed = reversed, reverseComplementRun = TRUE)

    class(extractedTrimmedRC) <- "extractedReadsListReverseComplement"
    return(extractedTrimmedRC)
}

.combineMarkerIdentifiedReadsLists <- function(extractedReadsList1, extractedReadsList2) {
    uniqueMarkers <- unique(c(names(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned), names(extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned)))
    combinedList <- list(n_reads = extractedReadsList1$n_reads + extractedReadsList2$n_reads)

    combinedList$identifiedMarkers <- structure(lapply(uniqueMarkers, function(j) {
        if (is.null(extractedReadsList2$identifiedMarkers[[j]])) return(extractedReadsList1$identifiedMarkers[[j]])
        else if (is.null(extractedReadsList1$identifiedMarkers[[j]])) return(extractedReadsList2$identifiedMarkers[[j]])
        else return(.appendExtractLists(extractedReadsList1$identifiedMarkers[[j]], extractedReadsList2$identifiedMarkers[[j]]))
    }), .Names = c(uniqueMarkers))

    combinedList$identifiedMarkersSequencesUniquelyAssigned <- structure(lapply(uniqueMarkers, function(j) {
        if (is.null(extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned[[j]])) return(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned[[j]])
        else if (is.null(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned[[j]])) return(extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned[[j]])
        else return(.appendExtractLists(extractedReadsList1$identifiedMarkersSequencesUniquelyAssigned[[j]], extractedReadsList2$identifiedMarkersSequencesUniquelyAssigned[[j]]))
    }), .Names = c(uniqueMarkers))

    class(combinedList) <- "extractedReadsListCombined"
    return(combinedList)
}

#' Control function for identifySTRRegions
#'
#' \code{identifySTRRegions.control}
#'
#' @return A control list setting default behaviour.
identifySTRRegions.control <- function(colList = NULL, numberOfThreads = 4L, reversed = TRUE,
                                        includeReverseComplement = TRUE, combineLists = TRUE, removeEmptyMarkers = TRUE, matchPatternMethod = "mclapply") {
    controlList <- list(colList = NULL, numberOfThreads = numberOfThreads, removeEmptyMarkers = removeEmptyMarkers,
                        reversed = reversed, includeReverseComplement = includeReverseComplement,
                        combineLists = combineLists, matchPatternMethod = matchPatternMethod)
    return(controlList)
}


.ShortReadQ.identifySTRRegions <- function(reads, flankingRegions, numberOfMutation, control = identifySTRRegions.control()) {
    seqs <- sread(reads)
    qual <- quality(reads)

    colID <- if(is.null(control$colList)) STRMPS:::.getCols(names(flankingRegions)) else colList
    flankSizes <- apply(flankingRegions[, c(colID$forwardCol, colID$reverseCol)], 2, nchar)

    extractAndTrimMarkerIdentifiedReads = switch(control$matchPatternMethod,
                                                 "mclapply" = STRMPS:::.extractAndTrimMarkerIdentifiedReads.mclapply,
                                                 "seqan" = STRMPS:::.extractAndTrimMarkerIdentifiedReads.seqAn)

    extractedSTRs <- extractAndTrimMarkerIdentifiedReads(seqs, qual, flankingRegions, colList = colID, numberOfMutation, control$removeEmptyMarkers,
                                                         flankSizes = flankSizes,
                                                         numberOfThreads = control$numberOfThreads, reversed = FALSE,
                                                         reverseComplementRun = FALSE, trimming = "flankshift")

    if (control$includeReverseComplement) {
        extractedSTRs_RC <- STRMPS:::.extractAndTrimMarkerIdentifiedReadsReverseComplement(seqs, qual, flankingRegions, colList = colID,
                                                                                  numberOfMutation = numberOfMutation,
                                                                                  removeEmptyMarkers = control$removeEmptyMarkers,
                                                                                  numberOfThreads = control$numberOfThreads,
                                                                                  reversed = control$reversed, trimming = "flankshift",
                                                                                  matchPatternMethod = control$matchPatternMethod)

        if (control$combineLists) {
            combinedRegions <- STRMPS:::.combineMarkerIdentifiedReadsLists(extractedSTRs, extractedSTRs_RC)
            return(combinedRegions)
        }
        else {
            extractedSTRs_Dual <- list(identifiedReads = extractedSTRs, identifiedReverseComplementReads = extractedSTRs_RC)
            class(extractedSTRs_Dual) <- "extractedReadsListNonCombined"
            return(extractedSTRs_Dual)
        }
    }
    else {
        return(extractedSTRs)
    }
}

.character.identifySTRRegions <- function(reads, flankingRegions, numberOfMutation, control = identifySTRRegions.control()) {
    if(!file.exists(reads))
        stop("File does not exist.")

    .ShortReadQ.identifySTRRegions(reads = readFastq(reads), flankingRegions = flankingRegions, numberOfMutation = numberOfMutation, control = control)
}

#' Identify the STR regions of a fastq-file or ShortReadQ-object.
#'
#' \code{identifySTRRegions} takes a fastq-file location or a ShortReadQ-object and identifies the STR regions
#' based on a directly adjacent flanking regions.
#' The function allows for mutation in the flanking regions through the numberOfMutation argument.
#'
#' @param reads either a fastq-file location or a ShortReadQ-object
#' @param flankingRegions containing marker ID/name, the directly adjacent forward and reverse flanking regions, used for identification.
#' @param control an \link{identifySTRRegions.control}-object.
#' @return The returned object is a list of lists. If the reverse complement strings are not included or if the \code{control$combineLists == TRUE},
#' a list, contains lists of untrimmed and trimmed strings for each row in \code{flankingRegions}. If \code{control$combineLists == FALSE}, the function returns a list of two such lists,
#' one for forward strings and one for the reverse complement strings.
setGeneric("identifySTRRegions", signature = "reads",
          function(reads, flankingRegions, numberOfMutation, control)
              standardGeneric("identifySTRRegions")
)

setMethod("identifySTRRegions", "ShortReadQ",
          function(reads, flankingRegions, numberOfMutation = 1, control = identifySTRRegions.control())
              .ShortReadQ.identifySTRRegions(reads = reads, flankingRegions = flankingRegions, numberOfMutation = numberOfMutation, control = control)
)

setMethod("identifySTRRegions", "character",
          function(reads, flankingRegions, numberOfMutation = 1, control = identifySTRRegions.control())
              .character.identifySTRRegions(reads = reads, flankingRegions = flankingRegions, numberOfMutation = numberOfMutation, control = control)
)

setClass("extractedReadsList")
setClass("extractedReadsListReverseComplement")
setClass("extractedReadsListCombined")
setClass("extractedReadsListNonCombined", representation(identifiedReads = "extractedReadsList", identifiedReverseComplementReads = "extractedReadsListReverseComplement"))
