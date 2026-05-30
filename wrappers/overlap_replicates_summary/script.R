#########################################
# wrapper for rule: overlap_replicates_summary
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[["run"]]
sink(logfile, append = T, type = "output")
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: overlap_replicates_summary \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(data.table)

beds = snakemake@input[['bed']]
tab = data.table()
for(bed in beds) {
  cat("## INFO: working on: ",bed,"\n")
  tab = rbind(tab, fread(bed))
}

cat("## printing summary table\n")
fwrite(tab, snakemake@output[["tab"]], sep='\t', row.names = F, col.names = T, quote = F)
