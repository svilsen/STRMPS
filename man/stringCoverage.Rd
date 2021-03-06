% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/allele_genotype_functions.R
\name{stringCoverage}
\alias{stringCoverage}
\title{Get string coverage STR identified objects.}
\usage{
stringCoverage(extractedReadsListObject, control = stringCoverage.control())
}
\arguments{
\item{extractedReadsListObject}{an extractedReadsList-object, created using the \link{identifySTRRegions}-function.}

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

# Limiting and restructuring
sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned),
                                function(m) which(m == flankingRegions$Marker))

# Aggregate the strings
stringCoverage(extractedReadsListObject = identifiedSTRs,
               control = stringCoverage.control(
                   motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
                   Type = flankingRegions$Type[sortedIncludedMarkers],
                   numberOfThreads = 1,
                   trace = FALSE,
                   simpleReturn = TRUE))
}
