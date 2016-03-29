#' Phred Quality Score
#'
#' @description Calculates phred quality score given probability.
#'
#' @param p probability of error.
#'
#' @return Phred quality score.
#'
#' @export
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
#'
#' @export
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
#'
#' @export
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
#'
#' @export
solexaQualityProbability <- function(q) {
    p <- 1/(10^(q/10) - 1)
    return(p)
}

createRanking <- function(x) {
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
