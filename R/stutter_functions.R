.findNeighbourStrings <- function(strings, alleles_i, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty) {
    motifDifference <- motifLength*abs(searchDirection)

    df <- vector("list", length(alleles_i))
    for(j in seq_along(alleles_i)) {
        alleles_j <- alleles_i[j]
        entireParentRepeatStructure <- LUS(strings[alleles_j, "Region"], motifLength, returnType = "fullList")

        alleleRepeatLength <- strings$Allele[alleles_j]
        neighbourRepeatLength <- strings$Allele[alleles_j] + searchDirection
        neighbours_j <- which(abs(strings$Allele - neighbourRepeatLength) < 1e-10)

        subMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -nchar(strings$Region[alleles_j]), baseOnly = FALSE)
        stutterAligned <- pairwiseAlignment(DNAStringSet(strings[neighbours_j, "Region"]), strings$Region[alleles_j], substitutionMatrix = subMatrix,
                                            gapOpening = -gapOpeningPenalty, gapExtension = -gapExtensionPenalty)

        trueStutters <- which(stutterAligned@score == (nchar(strings[alleles_j, "Region"]) - motifDifference - (gapOpeningPenalty + motifDifference*gapExtensionPenalty)))

        df_j <- vector("list", length = length(trueStutters))
        if (length(trueStutters) > 0) {
            for (k in seq_along(trueStutters)) {
                if ((j > 1) && any(strings[neighbours_j[trueStutters], "AlleleCalled"])) {
                    next
                }
                if (searchDirection == -1) {
                    missingRepeatUnit <- paste(unlist(strsplit(strings[alleles_j, "Region"], ""))[which(as.matrix(stutterAligned)[trueStutters[k],] == "-")], collapse = "")
                    missingRepeatUnitStartPosition <- which(unlist(strsplit(as.character(aligned(stutterAligned)[trueStutters[k]]), "")) == "-")[1]
                    entireParentRepeatStructure_k <- entireParentRepeatStructure[[missingRepeatUnit]]
                    occurenceInParent <- subset(entireParentRepeatStructure_k, Start <= missingRepeatUnitStartPosition & End >= missingRepeatUnitStartPosition + motifLength - 1L)$Repeats
                }
                else {
                    occurenceInParent = NA
                }

                if (j == 1 && length(alleles_i) > 1L) {
                    neighbour_fraction <- ifelse(strings$Allele[alleles_i[j + 1]] - strings$Allele[alleles_i[j]] > abs(searchDirection),
                                                 strings$Coverage[neighbours_j[trueStutters[k]]] / strings$Coverage[alleles_j],
                                                 strings$Coverage[neighbours_j[trueStutters[k]]] / (0.5*strings$Coverage[alleles_j] + 0.5*strings$Coverage[alleles_i[j + 1]]))
                }
                else {
                    neighbour_fraction <- strings$Coverage[neighbours_j[trueStutters[k]]] / strings$Coverage[alleles_j]
                }

                df_j[[k]] <- data.frame(Genotype = paste(strings$Allele[alleles_i], collapse = ",", sep = ""),
                                        ParentAllele = alleleRepeatLength,
                                        NeighbourAllele = neighbourRepeatLength,
                                        ParentLUS = as.character(strings$LUS[alleles_j]),
                                        ParentLUSLength = as.numeric(unlist(strsplit(strings$LUS[alleles_j], "]"))[length(unlist(strsplit(strings$LUS[alleles_j], "]")))]),
                                        NeighbourLUS = LUS(as.character(strings$Region[trueStutters[k]]), motifLength, returnType = "string"),
                                        NeighbourMissingUnit = paste("[", missingRepeatUnit, "]", sep = ""),
                                        NeighbourMissingUnitAppearanceParentCount = occurenceInParent,
                                        ParentCoverage = strings$Coverage[alleles_j],
                                        NeighbourCoverage = strings$Coverage[neighbours_j[trueStutters[k]]],
                                        NeighbourRatio = neighbour_fraction, stringsAsFactors = FALSE)
            }

            df[[j]] <- do.call(rbind, df_j)
        }
    }

    df_res <- do.call(rbind, df)
    if (is.null(df_res)) {
        df_res <- data.frame(Genotype = paste(strings$Allele[alleles_i], collapse = ",", sep = ""),
                             ParentAllele = alleleRepeatLength,
                             NeighbourAllele = NA,
                             ParentLUS = as.character(strings$LUS[alleles_j]),
                             ParentLUSLength = as.numeric(unlist(strsplit(strings$LUS[alleles_j], ""))[length(unlist(strsplit(strings$LUS[alleles_j], "")))]),
                             NeighbourLUS = NA, NeighbourMissingUnit = NA, NeighbourMissingUnitAppearanceParentCount = NA,
                             ParentCoverage = strings$Coverage[alleles_j],
                             NeighbourCoverage = NA, NeighbourRatio = NA, stringsAsFactors = FALSE)
    }

    return(df_res)
}

# motifLength = 4; searchDirection = -1; gapOpeningPenalty = 6; gapExtensionPenalty = 1; i = 1; j = 1; k = 1
.findNeighbours <- function(stringCoverageGenotypeListObject, motifLengths, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1) {
    if (length(motifLengths) == 1) {
        motifLengths <- rep(motifLengths, length(stringCoverageGenotypeListObject))
    }
    if (length(motifLengths) != length(stringCoverageGenotypeListObject)) {
        stop("'motifLengths' and 'stringCoverageGenotypeListObject' has to be of equal length.")
    }

    res <- vector("list", length(stringCoverageGenotypeListObject))
    for (i in seq_along(stringCoverageGenotypeListObject)) {
        strings <- stringCoverageGenotypeListObject[[i]]
        alleles_i <- which(strings$AlleleCalled)

        if (length(alleles_i) > 0) {
            df <- STRMPS:::.findNeighbourStrings(strings, alleles_i, motifLengths[i], searchDirection, gapOpeningPenalty, gapExtensionPenalty)
            res[[i]] <- data.frame(Marker = names(stringCoverageGenotypeListObject[i]), df, stringsAsFactors = FALSE)
        }
    }

    class(res) <- "neighbourList"
    return(res)
}

setClass("neighbourList")

#' @title Find stutters
#'
#' @export
setGeneric("findStutter", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, motifLength, searchDirection = -1, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
               standardGeneric("findStutter")
)

setMethod("findStutter", "stringCoverageGenotypeList",
          function(stringCoverageGenotypeListObject, motifLength = 4, searchDirection = -1, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
              .findNeighbours(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty)
)

#' @title Find left shoulder
#'
#' @export
setGeneric("findLeftShoulder", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
               standardGeneric("findLeftShoulder")
)

setMethod("findLeftShoulder", "stringCoverageGenotypeList",
          function(stringCoverageGenotypeListObject, motifLength = 4, searchDirection = -0.25, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
              .findNeighbours(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty)
)

#' @title Find right shoulder
#'
#' @export
setGeneric("findRightShoulder", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
               standardGeneric("findRightShoulder")
)

setMethod("findRightShoulder", "stringCoverageGenotypeList",
          function(stringCoverageGenotypeListObject, motifLength = 4, searchDirection = 0.25, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
              .findNeighbours(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty)
)
