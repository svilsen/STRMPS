###### Playground ######
if (FALSE) {
    dilution_files <- list.files("~/NGS/Data/Illumina/ForenSeqDIlution/ForenSeq 2016-01-22 Sensitivity panel A/FastQ", full.names = T)
    dilution_files_explained <- read.csv("~/NGS/Data/Illumina/ForenSeqDIlution/ForenSeq 2016-01-22 Sensitivity panel A/ForenSeq_20160114_sample_explanation.csv", sep = ",")

    getDilutionFile <- function(sampleNumber, run = 1, files = dilution_files) {
        s <- paste0("S", sampleNumber, "_")
        r <- paste0("R", run, "_")

        file = files[which(grepl(r, files) & grepl(s, files))]
        return(file)
    }
    flankingRegionsForenSeq <- structure(read.csv("~/AAU/PhD/Articles/2015/DropOut/R/STRFlankingRegions/forenseq.csv", sep = "\t", header = FALSE)[, 1:4], .Names = c("Marker", "Type", "Forward", "Reverse"))

    firstFile <- readFastq(getDilutionFile(1, 1))
    firstCase <- identifySTRRegions(reads = firstFile, flankingRegions = flankingRegionsForenSeq, nrOfMutations = 1, control = identifySTRRegions.control(includeReverseComplement = FALSE))
    firstCoverage <- stringCoverage(extractedReadsListObject = firstCase, control = stringCoverage.control(trace = T))

    stringCoverageListMarker <- firstCoverage[["D16S539"]]

    plotSequence.ggplot(firstCoverage, "D16S539")

}

