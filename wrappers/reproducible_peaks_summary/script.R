#########################################
# wrapper for rule: reproducible_peaks_summary
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[["run"]]
sink(logfile, append = T, type = "output")
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: reproducible_peaks_summary \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(data.table)

beds = snakemake@input[['bed']]
tab = data.table()
for(bed in beds) {
  cat("## INFO: working on: ",bed,"\n")
  path = dirname(bed)
  name = basename(bed)
  cond = tstrsplit(path,"/")[[3]]
  source = tstrsplit(path,"/")[[2]]
  source = substr(source,1,nchar(source)-6)
  rtype = ifelse(grepl('true_reps', name), 'true', 'pseudo')
  comp = ifelse(grepl('IDR', source), 'rep1_VS_rep2', 'all_reps')
  if(grepl('_VS_', name)) {
    comp = sub('.+\\.(rep[0-9]+_VS_rep[0-9]+)\\..+','\\1',name, perl=T)
  }
  tab = rbind(tab, data.table(condition=cond, comparison=comp, source=source, reps_type=rtype, N=fread(bed,sep='\t',header = F)[,.N]))
}
wtab = dcast(tab, condition + comparison + source ~ reps_type, value.var = 'N')
wtab[,ratio:=max(pseudo,true)/min(pseudo,true), by=seq_along(condition)]
wtab[,test:=fifelse(ratio>2, "failed", "passed")]

fwrite(wtab[,.(condition, comparison, source, N_true=true, N_pseudo=pseudo, rescue_ratio=ratio, reproducibility_test=test)],
       snakemake@output[["tab"]], sep='\t', row.names = F, col.names = T, quote = F)
