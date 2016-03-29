.findNeighbourStrings <- function(strings, alleles_i, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty) {
    df <- vector("list", length(alleles_i))
    for(j in seq_along(alleles_i)) {
        alleles_j <- alleles_i[j]
        entireParentRepeatStructure <- LUS(strings[alleles_j, "String"], motifLength, returnType = "fullList")

        alleleRepeatLength <- strings$Allele[alleles_j]
        neighbourRepeatLength <- strings$Allele[alleles_j] + searchDirection
        neighbours_j <- which(abs(strings$Allele - neighbourRepeatLength) < 1e-10)

        subMatrix <- nucleotideSubstitutionMatrix(match = 1, mismatch = -nchar(strings$String[alleles_j]), baseOnly = TRUE)
        stutterAligned <- pairwiseAlignment(DNAStringSet(strings$String[neighbours_j]), strings$String[alleles_j], substitutionMatrix = subMatrix,
                                            gapOpening = -gapOpeningPenalty, gapExtension = -gapExtensionPenalty)

        trueStutters <- which(stutterAligned@score == (motifLength*alleleRepeatLength - motifDifference - (gapOpeningPenalty + motifDifference*gapExtensionPenalty) +
                                                           (nchar(strings[alleles_j, "String"]) - nchar(strings[alleles_j, "STRRegion"]))))

        df_j <- vector("list", length = length(trueStutters))
        if (length(trueStutters) > 0) {
            for (k in seq_along(trueStutters)) {
                if (searchDirection == -1) {
                    missingRepeatUnit <- paste(unlist(strsplit(strings[alleles_j, "String"], ""))[which(as.matrix(stutterAligned)[trueStutters[k],] == "-")], collapse = "")
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

                df_j[[k]] <- data.frame(
                    Genotype = paste(strings$Allele[alleles_i], collapse = ",", sep = ""),
                    ParentAllele = alleleRepeatLength,
                    NeighbourAllele = neighbourRepeatLength,
                    ParentLUS = as.character(strings$LUS[alleles_j]),
                    NeighbourLUS = LUS(as.character(strings$String[trueStutters[k]]), returnType = "string"),
                    NeighbourMissingUnit = paste("[", missingRepeatUnit, "]", sep = ""),
                    NeighbourMissingUnitAppearanceParentCount = occurenceInParent,
                    ParentCoverage = strings$Coverage[alleles_j],
                    NeighbourCoverage = strings$Coverage[neighbours_j[trueStutters[k]]],
                    NeighbourRatio = neighbour_fraction)
            }

            df[[j]] <- do.call(rbind, df_j)
        }
    }

    df_res <- do.call(rbind, df)
    return(df_res)
}

# motifLength = 4; searchDirection = -1; gapOpeningPenalty = 6; gapExtensionPenalty = 1; i = 1; j = 1; k = 1
.findNeighbours <- function(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1) {
    res <- vector("list", length(stringCoverageGenotypeListObject))
    motifDifference <- motifLength*abs(searchDirection)

    for (i in seq_along(stringCoverageGenotypeListObject)) {
        strings <- stringCoverageGenotypeListObject[[i]]
        alleles_i <- which(strings$AlleleCalled)

        if (length(alleles_i) > 0) {
            df <- .findNeighbourStrings(strings, alleles_i, motifLength, searchDirection, gapOpeningPenalty, gapExtensionPenalty)
            res[[i]] <- cbind(Marker = i, MarkerName = names(stringCoverageGenotypeListObject[i]), df)
        }
    }

    res <- do.call(rbind, res)
    class(res) <- "neighbourDataFrame"
    return(res)
}

setClass("neighbourDataFrame")

#' @title Find stutters
#'
#' @export
setGeneric("findStutter", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, motifLength, searchDirection, gapOpeningPenalty = 6, gapExtensionPenalty = 1)
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
