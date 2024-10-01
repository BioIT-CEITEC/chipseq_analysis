#########################################
# wrapper for rule: diff_summary
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[["run"]]
sink(logfile, append = T)
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: diff_summary \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(data.table)

beds = snakemake@input[['bed']]
tab = data.table()
for(bed in beds) {
  cat("## INFO: working on: ",bed,"\n")
  btab = fread(bed, sep='\t', header = T)
  path = dirname(bed)
  name = basename(bed)
  comp = tstrsplit(path,"/")[[3]]
  source = tstrsplit(path,"/")[[2]]
  if(source == 'overlapped_peaks') {
    source = ifelse(grepl('by_MACS', name), "MACS_overlap", "MSPC_overlap")
    total = as.numeric(btab[,2])
    up = as.numeric(btab[,6])
    down = as.numeric(btab[,7])
    sig = up+down
  }
  if(source == 'MACS_bdgdiff') {
    total = as.numeric(btab[,2])
    up = as.numeric(btab[,3])
    down = as.numeric(btab[,4])
    sig = up+down
  }
  if(source == "DiffBind") {
    tab = rbind(tab, data.table(comparison=comp, approach=paste0(source,"_DESeq2"), total=btab$total, significant=btab$deseq_sig, sig_in_cond_1=btab$deseq_up, sig_in_cond_2=btab$deseq_dn))
    tab = rbind(tab, data.table(comparison=comp, approach=paste0(source,"_edgeR"), total=btab$total, significant=btab$edger_sig, sig_in_cond_1=btab$edger_up, sig_in_cond_2=btab$edger_dn))
  } else {
    tab = rbind(tab, data.table(comparison=comp, approach=source, total=total, significant=sig, sig_in_cond_1=up, sig_in_cond_2=down))
  }
}
fwrite(tab, snakemake@output[["tab"]], sep='\t', row.names = F, col.names = T, quote = F)
