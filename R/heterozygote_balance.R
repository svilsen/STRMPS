.stringCoverageGenotypeList.heterozygoteBalance <- function(stringCoverageGenotypeListObject, trueAllelesObject = NULL, beliefCol) {
    if (!is.null(trueAllelesObject))
        stop("'trueAllelesObject' is not NULL, but a 'stringCoverageGenotypeListObject' is supplied.
             Set 'trueAllelesObject' as NULL or use a 'stringCoverageGenotypeListObject' as the first argument.")

    ans <- vector("list", length = stringCoverageGenotypeListObject)
    for (i in seq_along(stringCoverageGenotypeListObject)) {
        alleles <- stringCoverageGenotypeListObject[[i]]
        heights <- alleles[alleles$AlleleCalled, beliefCol]
        calledAlleles <- alleles[alleles$AlleleCalled, "Allele"]

        ans[[i]] <- data.frame(Locus = names(stringCoverageGenotypeListObject)[i],
                               Allele1Belief = ifelse(length(heights) >= 1L, heights[1L], NA),
                               Allele2Belief = ifelse(length(heights) >= 2L, heights[2L], NA),
                               HeterozygosityBalance = ifelse(length(heights) == 2L, heights[2L] / heights[1L], NA),
                               RepeatNumberDifference = ifelse(length(heights) == 2L, calledAlleles[2L] - calledAlleles[1L], 0))
    }

    ans <- do.call(rbind, ans)
    names(ans) <- names(stringCoverageGenotypeListObject)
    class(ans) <- "heterozygoteBalanceUnknownList"
    return(ans)
}

.character.trueAlleles <- function(alleles, strings, markers, motifLength) {
    ans <- vector("list", length(markers))
    for (i in seq_along(markers)) {
        alleles_i <- as.numeric(unlist(strsplit(alleles[which(alleles$Marker == markers[i]), ], ",")))
        if (!all(alleles_i == floor(alleles_i))) {
            fixAllelesIndex <- which(alleles_i != floor(alleles_i))
            fixAlleles <- as.numeric(unlist(strsplit(as.character(alleles_i[fixAllelesIndex]), "\\.(?=\\d)", perl = TRUE)))
            alleles_i[fixAllelesIndex] <- fixAlleles[1] + fixAlleles[2]/motifLength
        }

        ans[[i]] <- list(Alleles = alleles_i, Strings = as.character(strings[[markers[i]]]))
    }

    names(ans) <- markers
    class(ans) <- "trueAllelesList"
    return(ans)
}

.list.trueAlleles <- function(alleles, strings, markers) {
    ans <- vector("list", length(markers))
    for (i in seq_along(markers)) {
        if (!is.numeric(alleles[[markers[i]]]))
            stop("The function assumes the elements of the allele list are numeric.")
        ans[[i]] <- list(Alleles = alleles[[markers[i]]], Strings = as.character(strings[[markers[i]]]))
    }

    names(ans) <- markers
    class(ans) <- "trueAllelesList"
    return(ans)
}

setClass("trueAllelesList")

#' Creates trueAllelesList-object.
#'
#' @description The function creates a list, with an element for each marker, containing two lists, with the true alleles and strings, respectively.
#'
#' @param alleles the true alleles provided as a list of numerics or a vector of strings with the alleles seperated by a ','.
#' @param strings the true strings of an individual (has to have same ordering as alleles) provided as 'character'.
#' @param markers the markers needed (use e.g. names of stringCoverage-object).
#' @param motifLength the length of the motif.
#'
#' @return An object of class \code{trueAllelesList}.
setGeneric("trueAlleles", signature = "alleles",
           function(alleles, strings, markers, motifLength)
               standardGeneric("trueAlleles")
)

setMethod("trueAlleles", "character",
          function(alleles, strings, markers, motifLength = 4)
              .character.trueAlleles(alleles, strings, markers, motifLength)
)

setMethod("trueAlleles", "list",
          function(alleles, strings, markers, motifLength = 4)
              .list.trueAlleles(alleles, strings, markers)
)


.stringCoverageList.heterozygoteBalance <- function(stringCoverageListObject, trueAllelesObject, beliefCol) {
    if (class(trueAllelesObject) != "trueAllelesList")
        stop("'trueAllelesObject' needs to be an object of class 'trueAllelesList'. See the 'trueAlleles'-function for more information.")

    ans <- vector("list", length(stringCoverageListObject))
    for (i in seq_along(stringCoverageListObject)) {
        alleles <- stringCoverageListObject[[i]]

        trueAlleles_i <- trueAllelesObejct[[i]]$Alleles
        trueStrings_i <- trueAllelesObejct[[i]]$Strings

        if (is.null(trueStrings_i)) {
            trueAlleleIndex <- which((alleles$Allele %in% trueAlleles_i))
        } else {
            trueAlleleIndex <- which((alleles$Allele %in% trueAlleles_i) && (alleles$String %in% trueStrings_i))
        }


        heights <- alleles[trueAlleleIndex, beliefCol]
        calledAlleles <- alleles[trueAlleleIndex, "Allele"]

        ans[[i]] <- data.frame(Locus = names(stringCoverageListObject)[i],
                               Allele1Belief = ifelse(length(heights) >= 1L, heights[1L], NA),
                               Allele2Belief = ifelse(length(heights) >= 2L, heights[2L], NA),
                               HeterozygosityBalance = ifelse(length(heights) == 2L, heights[2L] / heights[1L], NA),
                               RepeatNumberDifference = ifelse(length(heights) == 2L, calledAlleles[2L] - calledAlleles[1L], 0))
    }

    ans <- do.call(rbind, ans)
    names(ans) <- names(stringCoverageListObject)
    class(ans) <- "heterozygoteBalanceList"
    return(ans)
}

#' Finds the heterozygote balance.
#'
#' @description Given either a stringCoverageListObject and a list of the class 'trueAllelesList' or 'stringCoverageGenotypeListObject', the function finds if possible the heterozygote balance for every provided marker (the elements of the first list).
#'
#' @param genotypeStringCoverageListObject either a stringCoverageList-object or a stringCoverageGenotypeList-object.
#' @param trueAllelesObject a trueAlleles-object created using the \link{trueAlleles}-function.
#' @param beliefCol the coloumn in which we believe.
#'
#' @return The function returns an object of class 'heterozygoteBalanceList' or 'heterozygoteBalanceUnknownList'.
setGeneric("heterozygoteBalance", signature = "stringCoverageGenotypeListObject",
           function(stringCoverageGenotypeListObject, trueAllelesObject, beliefCol = "Coverage")
               standardGeneric("heterozygoteBalance")
)

setMethod("heterozygoteBalance", "stringCoverageGenotypeList",
          function(stringCoverageGenotypeListObject, trueAllelesObject = NULL, beliefCol = "Coverage")
              .stringCoverageGenotypeList.heterozygoteBalance(stringCoverageGenotypeListObject, trueAllelesObject, beliefCol)
)

setMethod("heterozygoteBalance", "stringCoverageList",
          function(stringCoverageGenotypeListObject, trueAllelesObject, beliefCol = "Coverage")
              .stringCoverageList.heterozygoteBalance(stringCoverageGenotypeListObject, trueAllelesObject, beliefCol)
)

setClass("heterozygoteBalanceUnknownList")
setClass("heterozygoteBalanceList")
