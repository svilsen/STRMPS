#' Phred Quality Score
#'
#' @description Calculates phred quality score given probability.
#'
#' @param p probability of error.
#'
#' @return Phred quality score.
phredQualityScore <- function(p) {
    q_p <- -10*log10(p)
    return(q_p)
}

#' Phred to Probability
#'
#' @description Converts Phred Quality Score to probability.
#'
#' @param q Phred Quality Score.
#'
#' @return Probability of error.
phredQualityProbability <- function(q) {
    p <- 10^(-q/10)
    return(p)
}

#' Solexa Quality Score
#'
#' @description Calculates Solexa Quality Score given probability.
#'
#' @param p probability of error.
#'
#' @return Solexa Quality Score.
solexaQualityScore <- function(p) {
    q_s <- -10*log10(p/(1 - p))
    return(q_s)
}

#' Solexa to Probability
#'
#' @description Converts Solexa Quality Score to probability.
#'
#' @param q Solexa Quality Score.
#'
#' @return Probability of error.
solexaQualityProbability <- function(q) {
    p <- 1/(10^(q/10) - 1)
    return(p)
}

.geometricMean <- function(x) {
    return(exp(1/length(x) * sum(log(x))))
}

.aggregateQuality <- function(q) {
    qS <- as(PhredQuality(as.character(q)), "matrix")
    qAvg <- PhredQuality(apply(qS, 2, function(x) .geometricMean(phredQualityProbability(x))))
    return(as.character(qAvg))
}

.createRanking <- function(x) {
    unique_x <- unique(x)
    rank = c()
    coverage = c()
    for (i in seq_along(unique_x)) {
        rank = c(rank, rep(i, length(which(x == unique_x[i]))))
        coverage = c(coverage, rep(length(which(x == unique_x[i])), length(which(x == unique_x[i]))))
    }

    final <- data.frame(Rank = rank, Coverage = coverage)
    return(final)
}

#' The weighted average probability of error for each base.
#'
#' @param averageQuality A string containing the phred encoded quality.
#' @param coverageString The coverage of the string assumed to be correct.
#' @param sumCoverage The sum of the coverage on the marker.
#' @param windowSize The size of the window used for weighting the probabilities.
#'
#' @return A vector containing the weighted probability of error for each base.
baseErrorProbability <- function(averageQuality, coverageString, sumCoverage, windowSize = 5) {
    averageQuality <- c(as(PhredQuality(averageQuality), "matrix"))
    rawProbabilityOfError <- phredQualityProbability(averageQuality)
    leftmostNeighbour <- ifelse(seq(1, length(rawProbabilityOfError)) - windowSize < 1, 1, seq(1, length(rawProbabilityOfError)) - windowSize)
    rightmostNeighbour <- ifelse(seq(1, length(rawProbabilityOfError)) + windowSize > length(rawProbabilityOfError), length(rawProbabilityOfError), seq(1, length(rawProbabilityOfError)) + windowSize)
    probabilityOfError <- sapply(seq_along(rawProbabilityOfError), function(i) {
        weightedBaseProbability <- rawProbabilityOfError / (abs(1:length(rawProbabilityOfError) - i) + 1)
        baseProbability <- max(weightedBaseProbability[leftmostNeighbour[i]:rightmostNeighbour[i]]) / sum(weightedBaseProbability)
        return(baseProbability * coverageString / sumCoverage)
    })

    return(probabilityOfError)
}
