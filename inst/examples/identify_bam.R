\donttest{
    library("Biostrings")
    library("ShortRead")
    library("Rsamtools")

    # Path to file
    bamReadPath <- system.file('extdata', "sampleSequences.bam", package = 'STRMPS')
    baiReadPath <- system.file('extdata', "sampleSequences.bam.bai", package = 'STRMPS')

    # Flanking regions
    data("flankingRegions")

    # Read the file into memory
    readFile <- scanBam(bamReadPath, baiReadPath)

    # Identify the STR's of the file, both readPath and readFile can be used.
    identifySTRRegions(reads = readFile,
                       flankingRegions = flankingRegions,
                       numberOfMutation = 1,
                       control = identifySTRRegions.control(numberOfThreads = 1)
    )
}
