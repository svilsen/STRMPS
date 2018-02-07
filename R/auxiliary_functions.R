## Used for suppressing dplyr/ggplot variable names...
globalVariables(c("AdjustedBasePairs", "AggregateQuality", "Allele", "BasePairs", "Coverage",
                  "ExpandedRegion", "ForwardFlank", "LUS", "Marker", "Motif", "MotifLength",
                  "NeighbourAllele", "Quality", "Region", "Repeats", "ReverseFlank", "Type"))

##
.appendExtractLists <- function (x, val, addRCIndex = TRUE, addRCIndexRef = "matchedSeq") {
    if (is.null(val)) {
        val <- structure(vector("list", length(x)), .Names = as.character(names(x)))
    }
    else if (is.null(x)) {
        x <- structure(vector("list", length(val)), .Names = as.character(names(val)))
    }

    xnames <- names(x)
    for (v in names(val)) {
        if (v == "name") {
            x[[v]] <- if (is.null(x[[v]])) as.character(val[[v]]) else  as.character(x[[v]])
        }
        else {
            if (is.null(x[[v]])) {
                x[[v]] <- val[[v]]
            }
            else if (is.null(val[[v]])) {
                x[[v]] <- x[[v]]
            }
            else {
                if (is_tibble(x[[v]]) & is_tibble(val[[v]])) {
                    x[[v]] <- bind_rows(x[[v]], val[[v]])
                }
                else {
                    x[[v]] <- append(x[[v]], val[[v]])
                }
            }
        }
    }
    if (addRCIndex & (addRCIndexRef %in% names(val))) {
        x[["ReverseComplement"]] <- data.frame(SeqID = c(x[[addRCIndexRef]]), ReverseComplement = rep(c(FALSE, TRUE), c(length(x[[addRCIndexRef]]) - length(val[[addRCIndexRef]]), length(val[[addRCIndexRef]]))))
    }

    return(x)
}

#' @title Block length of the missing motif.
#'
#' @description Given a motif length and a string it finds the blocks of the string.
#'
#' @param s a string of either class: 'character' or 'DNAString'.
#' @param motifLength the known motif length of the STR region.
#' @param returnType the type of return wanted. Takes three values 'numeric', 'string', or 'fullList' (or any other combination cased letters).
#'
#' @details If returnType is 'numeric', the function returns the numeric value of the LUS.
#' If returnType is instead chosen as 'string', the function returns "[AATG]x" i.e. the motif, AATG, is repeated 'x' times.
#' Lastly if the returnType is set to fullList, the function returns a list of data.frames containing every possible repeat structure their start and the numeric value of the repeat unit length.
#'
#' @return Depending on returnType it return an object of class 'numeric', 'string', or 'fulllist'.
BLMM <- function(s, motifLength = 4, returnType = "numeric") {
    motifLength <- if (!is.integer(motifLength)) as.integer(motifLength) else motifLength

    if (class(s) == "character") {
        sD <- DNAString(s)
    } else if (class(s) == "DNAString"){
        sD <- s
        s <- as.character(s)
    } else {
        stop("LUS only implemented for class 'character' or 'DNAString'")
    }

    typesOfMotifs <- oligonucleotideFrequency(sD, motifLength)
    typesOfMotifs <- typesOfMotifs[typesOfMotifs > 0]
    motifs <- names(typesOfMotifs)

    positionOfMotifs <- lapply(as.list(motifs), function(y, text = s) unlist(gregexpr2(y, text = text)))

    allRepeats <- structure(vector("list", length(positionOfMotifs)), .Names = motifs)
    for (i in seq_along(positionOfMotifs)) {
        y = positionOfMotifs[[i]]
        rleValues <- rle(y-(motifLength*(0:(length(y) - 1))))
        end <- y[cumsum(rleValues$lengths)] + motifLength

        startEndFrame <- data.frame(Start = end - rleValues$length*motifLength, End = end, Repeats = rleValues$length)
        j = 1
        while (j <= dim(startEndFrame)[1]) {
            whichExtending <- which(startEndFrame$Start == startEndFrame$End[j])
            if (length(whichExtending) > 0) {
                startEndFrame$End[j] <- startEndFrame$End[whichExtending]
                startEndFrame$Repeats[j] <- startEndFrame$Repeats[j] + startEndFrame$Repeats[whichExtending]
                startEndFrame <- startEndFrame[-whichExtending, ]
            }

            j = j + 1
        }

        allRepeats[[i]] <- startEndFrame
    }


    if (length(allRepeats) == 0) {
        allRepeats <- list("NA" = data.frame(Start = NA, End = NA, Repeats = NA))
    }
    allRepeats <- enframe(allRepeats, name = "Motif") %>% unnest()

    reducedRepeats <- allRepeats
    i = 1
    while ((i <= dim(reducedRepeats)[1]) & !is.na(reducedRepeats$Repeats[1])) {
        whichWithin <- rep(FALSE, dim(reducedRepeats)[1])
        whichWithin[-i] <- IRanges(reducedRepeats$Start[-i], reducedRepeats$End[-i]) %within% IRanges(reducedRepeats$Start[i], reducedRepeats$End[i])

        reducedRepeats <- reducedRepeats[!whichWithin, ]
        i = i + 1
    }

    lusOfMotifs <- reducedRepeats %>% group_by(Motif) %>% filter(Repeats == max(Repeats)) %>% ungroup()
    lus <- which.max(lusOfMotifs$Repeats)

    if (tolower(returnType) == "numeric") {
        numeric_format <- lusOfMotifs$Repeats[lus]
        return(numeric_format)
    }
    else if (tolower(returnType) == "string") {
        string_format <- paste("[", lusOfMotifs$Motif[lus], "]", lusOfMotifs$Repeats[lus], sep="")
        return(string_format)
    }
    else if (tolower(returnType) == "fulllist") {
        return(reducedRepeats)
    }
    else
        stop(paste(returnType, "is not valid. Please use 'numeric', 'string', or 'fullList'."))
}

.cyclicRotation <- function(x, y) {
    (nchar(y) == nchar(x)) && (grepl(y, strrep(x, 2), fixed = TRUE))
}

.loadRData <- function(fileName) {
    load(fileName)
    get(ls()[ls() != "fileName"])
}
