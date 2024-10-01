library(rmspc)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
print(args)

# Tec 0.001 0.01 0.05 10 outputDir inputs
replicateType = args[1]
stringencyThreshold = as.numeric(args[2])
weakThreshold = as.numeric(args[3])
alpha = as.numeric(args[4])
keep = TRUE
GRanges = FALSE
multipleIntersections = "Lowest"
degreeOfParallelism = as.numeric(args[5])
c = 2
outputBed = args[6]
outputPath = dirname(args[6])
inputs = args[7:length(args)]
print(inputs)

results <- mspc(
  input = inputs, 
  replicateType = replicateType,
  stringencyThreshold = stringencyThreshold,
  weakThreshold = weakThreshold, 
  keep = keep,
  GRanges = GRanges,
  multipleIntersections = multipleIntersections,
  degreeOfParallelism = degreeOfParallelism,
  outputPath = outputPath,
  c = c,
  alpha = alpha)

tab = fread(paste0(outputPath, "ConsensusPeaks.bed"), sep = "\t", col.names = c("chr","start","end","peak_id","pvalue"))
tab