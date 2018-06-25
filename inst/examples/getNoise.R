# Strings aggregated by 'stringCoverage()'
data("stringCoverageList")

noiseList <- identifyNoise(stringCoverageList, thresholdSignal = 0.03)
