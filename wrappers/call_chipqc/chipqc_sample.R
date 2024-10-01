library(ChIPQC)
options("browser"="false")

args = commandArgs(trailingOnly=TRUE)
print(args)
outdir = args[1]
html_prefix = args[2]
sample_data = args[3]
peaks = args[4]
reads = args[5]

print("# Loading data into ChIPQCsample object")
sample = ChIPQCsample(reads, peaks)
print("# Saving ChIPQCsample object")
save(sample, file=sample_data)
print("# Producing ChIPQC report")
ChIPQCreport(sample, reportName=html_prefix, reportFolder=outdir)