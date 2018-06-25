readPath <- system.file('extdata', 'sampleSequences.fastq', package = 'STRMPS')

STRMPSWorkflow(readPath, control = workflow.control(restrictType = "Autosomal"))
