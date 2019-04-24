## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- message = FALSE----------------------------------------------------
library("Biostrings")
library("ShortRead")

read_path <- system.file('extdata', 'sampleSequences.fastq', package = 'STRMPS')
sequences <- readFastq(read_path)

## ------------------------------------------------------------------------
sequences@sread

## ------------------------------------------------------------------------
library("STRMPS")
data("flankingRegions")

## ------------------------------------------------------------------------
head(flankingRegions, 5)

## ----iden, cache = TRUE--------------------------------------------------
identifiedSTRs <- identifySTRRegions(reads = sequences,
                                     flankingRegions = flankingRegions,
                                     numberOfMutation = 1,
                                     control = identifySTRRegions.control(
                                            numberOfThreads = 1,
                                            includeReverseComplement = FALSE
                                         )
                                     )

## ------------------------------------------------------------------------
names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned$CSF1PO)

## ------------------------------------------------------------------------
sortedIncludedMarkers <- sapply(names(identifiedSTRs$identifiedMarkersSequencesUniquelyAssigned),
                                function(m) which(m == flankingRegions$Marker))

## ----aggr, cache = TRUE--------------------------------------------------
stringCoverageList <- stringCoverage(extractedReadsListObject = identifiedSTRs,
                                     motifLength = flankingRegions$MotifLength[sortedIncludedMarkers],
                                     Type = flankingRegions$Type[sortedIncludedMarkers],
                                     control = stringCoverage.control(
                                            numberOfThreads = 1,
                                            trace = FALSE,
                                            simpleReturn = TRUE
                                         )
                                     )

## ------------------------------------------------------------------------
stringCoverageList$CSF1PO

## ------------------------------------------------------------------------
genotypeList <- getGenotype(stringCoverageList)

## ------------------------------------------------------------------------
genotypeList$CSF1PO

## ------------------------------------------------------------------------
noiseList <- identifyNoise(stringCoverageList, thresholdSignal = 0.03)

## ------------------------------------------------------------------------
noiseList$CSF1PO

## ------------------------------------------------------------------------
BLMM(stringCoverageList$CSF1PO$Region[1], 
     motifLength = stringCoverageList$CSF1PO$MotifLength[1], 
     returnType = "numeric")

## ------------------------------------------------------------------------
BLMM(stringCoverageList$CSF1PO$Region[1], 
     motifLength = stringCoverageList$CSF1PO$MotifLength[1], 
     returnType = "string")

## ------------------------------------------------------------------------
fullBLMM <- BLMM(stringCoverageList$CSF1PO$Region[1], 
                 motifLength = stringCoverageList$CSF1PO$MotifLength[1], 
                 returnType = "fulllist")
head(fullBLMM, 5)

## ------------------------------------------------------------------------
stringCoverageGenotypeList <- mergeGenotypeStringCoverage(stringCoverageList, genotypeList)

stutterList <- findStutter(stringCoverageGenotypeList)
stutterTibble <- subset(do.call("rbind", stutterList), !is.na(Genotype))
head(stutterTibble, 5)

## ---- eval = FALSE-------------------------------------------------------
#  STRMPSWorkflow(read_path,
#                 control = workflow.control(
#                      restrictType = "Autosomal",
#                      numberOfThreads = 1
#                     )
#                 )

