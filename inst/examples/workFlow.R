read_path <-
    paste(system.file('sampleFiles', package = 'STRMPS'),
          "sampleSequences.fastq", sep = "/")

stringCoverageListTrimmed <- STRMPSWorkflow(read_path,
                                            control = workflow.control(
                                                restrictType = "Autosomal"
                                            ))
