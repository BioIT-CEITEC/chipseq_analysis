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

print("# Creating output directory")
dir.create(outputPath, recursive = T)

print("# calling MSPC")
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

print("# postprocessing output BED")
tab = fread(paste0(outputPath, "/ConsensusPeaks.bed"), sep = "\t", col.names = c("chr","start","end","peak_id","score"))
tab[,chr:=sub("^chr","",chr)]
tab[,strand:="."]
tab[,signal:=1]
tab[,pvalue:=score]
tab[,qvalue:=-1*log10(p.adjust(10^-pvalue,"BH"))]
tab[,summit:=floor((end - start)*0.5)]
tab[,score:=floor(10*qvalue)]
fwrite(tab, outputBed, sep = "\t", quote = F, row.names = F, col.names = F)