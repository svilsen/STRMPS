% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allele_genotype_functions.R
\name{stringCoverage,extractedReadsList-method}
\alias{stringCoverage,extractedReadsList-method}
\title{Get string coverage STR identified objects.}
\usage{
\S4method{stringCoverage}{extractedReadsList}(
  extractedReadsListObject,
  flankingRegions,
  control = stringCoverage.control()
)
}
\arguments{
\item{extractedReadsListObject}{an extractedReadsList-object, created using the \link{identifySTRRegions}-function.}

\item{flankingRegions}{containing marker ID/name, the directly adjacent forward and reverse flanking regions, used for identification.}

\item{control}{an \link{stringCoverage.control}-object.}
}
\value{
Returns a list, with an element for every marker in extractedReadsList-object, each element contains the string coverage of all unique strings of a given marker.
}
\description{
\code{stringCoverage} takes an extractedReadsList-object and finds the coverage of every unique string for every marker in the provided list.
}
\examples{
# Regions identified using 'identifySTRs()'
data("identifiedSTRs")

# Aggregate the strings
stringCoverage(
    extractedReadsListObject = identifiedSTRs,
    flankingRegions = flankingRegions,
    control = stringCoverage.control(
        numberOfThreads = 1,
        trace = FALSE,
        simpleReturn = TRUE
    )
)
}
