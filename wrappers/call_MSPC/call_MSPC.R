library(rmspc)
library(data.table)

args = commandArgs(trailingOnly=TRUE)
print(args)

# Tec 0.001 0.01 0.05 10 outputDir inputs
replicateType = args[1]
stringencyThreshold = as.numeric(args[2])
weakThreshold = as.numeric(args[3])
alpha = as.numeric(args[4])
keep = FALSE
GRanges = TRUE
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
  outputPath = NULL,
  c = c,
  alpha = alpha)

print("# MSPC finished, writinng down the results:")
fwrite(as.data.table(results$GRangesObjects$ConsensusPeaks),
       paste0(outputPath,"/ConsensusPeaks.bed"), 
       sep="\t", 
       row.names=F, 
       col.names=T)

print("# postprocessing output table")
tab = copy(as.data.table(results$GRangesObjects$ConsensusPeaks))[,.(
  #chr=sub("^chr","",seqnames),   # Not sure why did I set this hard-coded name trimming, might be usefull in some cases
  chr=seqnames,
  start, 
  end, 
  name, 
  score, 
  strand=".", 
  signal=1, 
  pvalue=score,
  qvalue=-1*log10(p.adjust(10^-score,"BH")), 
  summit=floor((end - start)*0.5)
  )]
tab[,score:=floor(10*qvalue)]
fwrite(tab, outputBed, sep = "\t", quote = F, row.names = F, col.names = F)