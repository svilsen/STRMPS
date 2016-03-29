###### Playground ######
if (FALSE) {
    load("~/Dropbox/PhD/Rpackage/NewCode/noiseDiscLap.RData")
    read_path <- "~/NGS/Data/LT-IonTorrentPGM/STR-10plex/Fortyndinger/IonXpress_002_R_2013_11_20_07_47_42_user_SN2-44-Fred-HSM_2ng-50pg_STRpool1_Auto_user_SN2-44-Fred-HSM_2ng-50pg_STRpool1_89.fastq"
    #read_path <- "~/NGS/Data/LT-IonTorrentPGM/SomalierPopulationSTR10Plex/temp/rawlib.basecaller/PGM_318C_IonXpress_001_1.fastq"
    readfile <- readFastq(read_path)
    flankingRegions <- read.table("~/NGS/Scripts/R/STRflanks.txt", header=FALSE, col.names=c("Locus", "Forward.Flank", "Reverse.Flank"), stringsAsFactors = FALSE)
    str_10_plex <- c("TPOX", "CSF1PO", "D5S818", "D7S820", "D16S539", "D3S1358", "D8S1179", "vWA", "TH01")
    flankingRegions_10plex <- flankingRegions[which(flankingRegions$Locus %in% str_10_plex), ]

    nrThreads = 4L

    #identifiedFlanksObj <- STRMPS:::.identifyFlankingRegions(sread(reads), flankingRegions_10plex)
    identifiedSTRs <- identifySTRRegions(reads = read_path, flankingRegions = flankingRegions_10plex, nrOfMutations = 1, control = identifySTRRegions.control())
    #extractedReadsListObject <- identifiedSTRs

    stringCoverageListObject <- stringCoverage(extractedReadsListObject = identifiedSTRs, control = stringCoverage.control(trace = T))
    head(stringCoverageListObject[[1]], 50)
    noiseGenotypeIdentifiedListObject <- getGenotype(stringCoverageListObject)
    stringCoverageGenotypeListObject <- mergeGenotypeStringCoverage(stringCoverageListObject, noiseGenotypeIdentifiedListObject)

    head(stringCoverageGenotypeListObject[[1]], 15)
}

