## .mclapplyvmatchPattern and .vmatchMultiplePatternsAlternate are used for diagnostic purposes only
.mclapplyvmatchPattern <- function(flanks, seqs, max.mismatch, numberOfThreads) {
    mclapply(flanks, function(f) vmatchPattern(f, seqs, max.mismatch = max.mismatch), mc.cores = numberOfThreads)
}

.vmatchMultiplePatternsAlternate <- function(flanks, seqs, max.mismatch = 1, numberOfThreads) {
    numberOfThreads <- if(is.integer(numberOfThreads)) numberOfThreads else as.integer(numberOfThreads)
    max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
    min.mismatch <- Biostrings:::normargMaxMismatch(0)
    with.indels <- Biostrings:::normargWithIndels(FALSE)
    fixed <- Biostrings:::normargFixed(TRUE, seqs)
    alg <- Biostrings:::selectAlgo("auto", flanks, max.mismatch, min.mismatch, with.indels = with.indels, fixed = fixed)

    match_s <- vmatchMultiPatternAlternate(flanks, seqs, max_mismatch = max.mismatch, min_mismatch = min.mismatch,
                                           with_indels = with.indels, fixed = fixed, algorithm = alg,
                                           matches_as = "MATCHES_AS_ENDS", envir = NULL)

    resList <- lapply(seq_along(match_s), function(i) new("ByPos_MIndex", width0=rep.int(width(flanks[i]), length(seqs)), NAMES=names(seqs), ends=match_s[[i]]))
    return(resList)
}

.vmatchMultiplePatterns <- function(flanks, seqs, max.mismatch = 1, numberOfThreads) {
    numberOfThreads <- if(is.integer(numberOfThreads)) numberOfThreads else as.integer(numberOfThreads)
    max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)
    min.mismatch <- Biostrings:::normargMaxMismatch(0)
    with.indels <- Biostrings:::normargWithIndels(FALSE)
    fixed <- Biostrings:::normargFixed(TRUE, seqs)
    alg <- Biostrings:::selectAlgo("auto", flanks, max.mismatch, min.mismatch, with.indels = with.indels, fixed = fixed)

    dataSplit <- if(numberOfThreads == 1) 1:length(flanks) else unname(split(1:length(flanks), cut(1:length(flanks), numberOfThreads, labels = FALSE)))
    match_s <- unlist(mclapply(dataSplit, function(d) vmatchMultiPattern(flanks[d], seqs, max_mismatch = max.mismatch, min_mismatch = min.mismatch,
                                   with_indels = with.indels, fixed = fixed, algorithm = alg,
                                   matches_as = "MATCHES_AS_ENDS", envir = NULL), mc.cores = numberOfThreads), recursive = FALSE)

    match_s <- split(match_s, rep(1:length(flanks), each = length(seqs)))
    resList <- lapply(seq_along(match_s), function(i) new("ByPos_MIndex", width0=rep.int(width(flanks[i]), length(seqs)), NAMES=names(seqs), ends=match_s[[i]]))
    return(resList)
}

.vmatchMultiplePatternsSeqAn <- function(flanks, seqs, max.mismatch = 1, numberOfThreads) {
    flanks_char <- as.character(flanks)
    seqs_char <- as.character(seqs)
    max.mismatch <- Biostrings:::normargMaxMismatch(max.mismatch)

    dataSplit <- if(numberOfThreads == 1) 1:length(flanks) else split(1:length(flanks), cut(1:length(flanks), numberOfThreads, labels = FALSE))
    match_s <- structure(unlist(mclapply(dataSplit, function(i) vmatchMultiPatternSeqAn(flanks_char[i], seqs_char, max.mismatch), mc.cores = numberOfThreads), recursive = FALSE), .Names = c())

    resList <- lapply(seq_along(match_s), function(i) new("ByPos_MIndex", width0=rep.int(width(flanks[i]), length(seqs)), NAMES=names(seqs), ends=match_s[[i]]))
    return(resList)
}

.getCols <- function(colNames) {
    structure(lapply(c("marker", "forward", "reverse"), function(n) grep(n, tolower(colNames))), .Names = c("markerCol", "forwardCol", "reverseCol"))
}

.identifyFlankingRegions <- function(seqs, flankingRegions, matchPatternMethod = c("vmatch", "seqan"),
                                     colList, nrOfMutations = 1L, numberOfThreads = 4, removeEmptyMarkers = TRUE) {
    forwardFlank <- DNAStringSet(flankingRegions[, colList$forwardCol])
    reverseFlank <- DNAStringSet(flankingRegions[, colList$reverseCol])

    matchPatternMethod <- match.arg(matchPatternMethod)
    parallelvmatchPattern <- switch(tolower(matchPatternMethod), seqan = .vmatchMultiplePatternsSeqAn, vmatch = .vmatchMultiplePatterns)

    identifiedForward <- parallelvmatchPattern(forwardFlank, seqs, nrOfMutations, numberOfThreads)
    identifiedReverse <- parallelvmatchPattern(reverseFlank, seqs, nrOfMutations, numberOfThreads)

    matchedSeq <- lapply(1:nrow(flankingRegions), function(i) which(elementLengths(identifiedForward[[i]]) >= 1L & elementLengths(identifiedReverse[[i]]) >= 1L))
    if(removeEmptyMarkers) {
        nonEmptyEntries <- which(sapply(matchedSeq, length) != 0)
        rList <- list(markers = flankingRegions[nonEmptyEntries, colList$markerCol], identifiedForward = identifiedForward[nonEmptyEntries], identifiedReverse = identifiedReverse[nonEmptyEntries], matchedSeq = matchedSeq[nonEmptyEntries])
    }
    else {
        rList <- list(markers = flankingRegions[, colList$markerCol], identifiedForward = identifiedForward, identifiedReverse = identifiedReverse, matchedSeq = matchedSeq)
    }
    return(rList)
}

.extractAndTrimMarkerIdentifiedReads <- function(seqs, qual, identifiedFlanksObj, flankSizes, flankShift = NULL, numberOfThreads = 4, reversed = FALSE, reverseComplement = FALSE) {
    if (is.null(flankShift))
        flankShift <- data.frame(matrix(0, nrow = length(identifiedFlanksObj$matchedSeq), ncol = 2))

    identifiedMarkers <- structure(lapply(seq_along(identifiedFlanksObj$matchedSeq), function(i) {
        marker <- identifiedFlanksObj$markers[i]
        identifiedForward <- identifiedFlanksObj$identifiedForward[[i]]
        identifiedReverse <- identifiedFlanksObj$identifiedReverse[[i]]
        matchedSeq <- identifiedFlanksObj$matchedSeq[[i]]

        # F: Take the first in the IRanges
        endForward <- unlist(lapply(endIndex(identifiedForward)[matchedSeq], function(v) v[1L]))
        startForward <- endForward - (flankSizes[i, 1] + flankShift[i, 1] - 1)

        # R: Take the last in the IRanges
        startReverse <- unlist(lapply(startIndex(identifiedReverse)[matchedSeq], function(v) v[length(v)]))
        endReverse <- startReverse + (flankSizes[i, 2] - flankShift[i, 2] - 1)

        # Removes reads were forward is NOT observed before reverse
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
            trimmed <- subseq(seqs[matchedSeq], start = endForward + 1, end = startReverse - 1)
            trimmedIncludingFlanks <- subseq(seqs[matchedSeq], start = startForward, end = endReverse)
            trimmedQuality <- subseq(quality(qual[matchedSeq]), start = endForward + 1, end = startReverse - 1)

            if ((reversed && reverseComplement)) {
                trimmed <- reverseComplement(trimmed)
                trimmedIncludingFlanks <- reverseComplement(trimmedIncludingFlanks)
                trimmedQuality <- reverse(trimmedQuality)
            }

            rList <- list(name = marker, matchedSeq = matchedSeq, startForward = startForward, endForward = endForward, startReverse = startReverse, endReverse = endReverse,
                          trimmed = trimmed, trimmedIncludingFlanks = trimmedIncludingFlanks, trimmedQuality = trimmedQuality)
        }

        return(rList)
    }), .Names = identifiedFlanksObj$markers)

    # Remove sequences occuring at several loci
    identifiedSeqTable <- sort(table(unlist(lapply(identifiedMarkers, function(r) r$matchedSeq))), decreasing = TRUE)
    identifiedSeqNotUnique <- as.integer(names(identifiedSeqTable[identifiedSeqTable > 1L]))
    identifiedMarkersUniqueSequences <- lapply(identifiedMarkers, function(r) {
        keepSeq <- which(!(r$matchedSeq %in% identifiedSeqNotUnique))
        return(list(name = r$name, matchedSeq = r$matchedSeq[keepSeq], endForward = r$endForward[keepSeq], startReverse = r$startReverse[keepSeq],
                    trimmed = r$trimmed[keepSeq], trimmedIncludingFlanks = r$trimmedIncludingFlanks[keepSeq], trimmedQuality = r$trimmedQuality[keepSeq]))
    })

    extractedReadsList <- list(n_reads = length(seqs), reverseComplement = reverseComplement, identifiedMarkers = identifiedMarkers, identifiedMarkersSequencesUniquelyAssigned = identifiedMarkersUniqueSequences)
    class(extractedReadsList) <- "extractedReadsList"
    return(extractedReadsList)
}

.extractAndTrimMarkerIdentifiedReadsReverseComplement <- function(seqs, qual, flankingRegions, colList = NULL, matchPatternMethod, nrOfMutations, numberOfThreads = numberOfThreads, removeEmptyMarkers = TRUE, reversed = TRUE) {
    flankingRegionsReverseComplement <- structure(data.frame(flankingRegions[, colID$markerCol],
                                                             as.character(reverseComplement(DNAStringSet(flankingRegions[, colList$reverseCol]))),
                                                             as.character(reverseComplement(DNAStringSet(flankingRegions[, colList$forwardCol]))),
                                                             stringsAsFactors = FALSE),
                                                  .Names = c("Marker", "ForwardRC", "ReverseRC"))

    identifiedFlanksReverseComplement <- .identifyFlankingRegions(seqs, flankingRegionsReverseComplement, matchPatternMethod = matchPatternMethod,
                                                                  colList = list(markerCol = 1, forwardCol = 2, reverseCol = 3), nrOfMutations = nrOfMutations,
                                                                  numberOfThreads = numberOfThreads, removeEmptyMarkers = removeEmptyMarkers)

    extractedTrimmedRC <- .extractAndTrimMarkerIdentifiedReads(seqs, qual, identifiedFlanksReverseComplement, flankSizes = apply(flankingRegionsReverseComplement[, -1], 2, nchar), numberOfThreads = numberOfThreads, reversed = reversed, reverseComplement = TRUE)
    class(extractedTrimmedRC) <- "extractedReadsListReverseComplement"
    return(extractedTrimmedRC)
}

.combineMarkerIdentifiedReadsLists <- function(extractedReadsList1, extractedReadsList2) {
    uniqueMarkers <- unique(c(names(extractedReadsList1$identifiedMarkers), names(extractedReadsList2$identifiedMarkers)))
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
#' @return A list containing...
#' @export
identifySTRRegions.control <- function(colList = NULL, numberOfThreads = 4L, reversed = TRUE,
                                        includeReverseComplement = TRUE, combineLists = TRUE, removeEmptyMarkers = TRUE, matchPatternMethod = "vmatch") {
    controlList <- list(colList = NULL, numberOfThreads = numberOfThreads, removeEmptyMarkers = removeEmptyMarkers,
                        reversed = reversed, includeReverseComplement = includeReverseComplement,
                        combineLists = combineLists, matchPatternMethod = matchPatternMethod)
    return(controlList)
}

.ShortReadQ.identifySTRRegions <- function(reads, flankingRegions, nrOfMutations, control = identifySTRRegions.control()) {
    seqs <- sread(reads)
    qual <- quality(reads)

    colID <- if(is.null(colList)) .getCols(names(flankingRegions)) else colList

    identifiedRegions <- .identifyFlankingRegions(seqs, flankingRegions, matchPatternMethod = control$matchPatternMethod,
                                                  colList = colID, nrOfMutations = nrOfMutations,
                                                  numberOfThreads = control$numberOfThreads, removeEmptyMarkers = control$removeEmptyMarkers)
    extractedSTRs <- .extractAndTrimMarkerIdentifiedReads(seqs, qual, identifiedFlanksObj = identifiedRegions, flankSizes = apply(flankingRegions[, c(colID$forwardCol, colID$reverseCol)], 2, nchar),
                                                          numberOfThreads = control$numberOfThreads, reverseComplement = FALSE)

    if (control$includeReverseComplement) {
        extractedSTRs_RC <- .extractAndTrimMarkerIdentifiedReadsReverseComplement(seqs, qual, flankingRegions, colList = colID,
                                                                                  matchPatternMethod = control$matchPatternMethod,
                                                                                  nrOfMutations = control$nrOfMutations,
                                                                                  numberOfThreads = control$numberOfThreads,
                                                                                  removeEmptyMarkers = control$removeEmptyMarkers,
                                                                                  reversed = control$reversed)

        if (control$combineLists) {
            combinedRegions <- .combineMarkerIdentifiedReadsLists(extractedSTRs, extractedSTRs_RC)
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

.character.identifySTRRegions <- function(reads, flankingRegions, nrOfMutations, control = identifySTRRegions.control()) {
    if(!file.exists(reads))
        stop("File does not exist.")

    .ShortReadQ.identifySTRRegions(reads = readFastq(reads), flankingRegions = flankingRegions, nrOfMutations = nrOfMutations, control = control)
}

#' Identify the STR regions of a fastq-file or ShortReadQ-object.
#'
#' \code{identifySTRRegions} takes a fastq-file location or a ShortReadQ-object and identifies the STR regions
#' based on a directly adjacent flanking regions.
#' The function allows for mutation in the flanking regions through the nrOfMutations argument.
#'
#' @param reads either a fastq-file location or a ShortReadQ-object
#' @param flankingRegions containing marker ID/name, the directly adjacent forward and reverse flanking regions, used for identification.
#' @param control an \link{identifySTRRegions.control}-object.
#' @return The returned object is a list of lists. If the reverse complement strings are not included or if the \code{control$combineLists == TRUE},
#' a list, contains lists of untrimmed and trimmed strings for each row in \code{flankingRegions}. If \code{control$combineLists == FALSE}, the function returns a list of two such lists,
#' one for forward strings and one for the reverse complement strings.
#' @export
setGeneric("identifySTRRegions", signature = "reads",
          function(reads, flankingRegions, nrOfMutations, control)
              standardGeneric("identifySTRRegions")
)

setMethod("identifySTRRegions", "ShortReadQ",
          function(reads, flankingRegions, nrOfMutations = 1, control = identifySTRRegions.control())
              .ShortReadQ.identifySTRRegions(reads = reads, flankingRegions = flankingRegions, nrOfMutations = nrOfMutations, control = control)
)

setMethod("identifySTRRegions", "character",
          function(reads, flankingRegions, nrOfMutations = 1, control = identifySTRRegions.control())
              .character.identifySTRRegions(reads = reads, flankingRegions = flankingRegions, nrOfMutations = nrOfMutations, control = control)
)

setClass("extractedReadsList")
setClass("extractedReadsListReverseComplement")
setClass("extractedReadsListCombined")
setClass("extractedReadsListNonCombined", representation(identifiedReads = "extractedReadsList", identifiedReverseComplementReads = "extractedReadsListReverseComplement"))
