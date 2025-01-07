library(data.table)
library(stringr)

#setwd("/mnt/ssd/ssd_1/workspace/martin/massari_cutnrun_spike_test")
#args = c(
#  "results/SEACR_peaks/1Ex10_KO_STAT1_1_VS_NT2_IgG_1/1Ex10_KO_STAT1_1_VS_NT2_IgG_1.keep_dups.peaks.all.narrowPeak",
#  "results/SEACR_peaks/1Ex10_KO_STAT1_2_VS_NT2_IgG_2/1Ex10_KO_STAT1_2_VS_NT2_IgG_2.keep_dups.peaks.all.narrowPeak",
#  "results/overlapped_replicates/1Ex10_KO_STAT1/1Ex10_KO_STAT1.merged_peaks.keep_dups.by_SEACR.bed",
#  "results/overlapped_replicates/1Ex10_KO_STAT1/1Ex10_KO_STAT1.overlapped_peaks.keep_dups.by_SEACR.bed",
#  "results/overlapped_replicates/1Ex10_KO_STAT1/1Ex10_KO_STAT1.singletons_in_rep1.keep_dups.by_SEACR.tsv",
#  "results/overlapped_replicates/1Ex10_KO_STAT1/1Ex10_KO_STAT1.singletons_in_rep2.keep_dups.by_SEACR.tsv",
#  "results/overlapped_replicates/1Ex10_KO_STAT1/1Ex10_KO_STAT1.summary_table.keep_dups.by_SEACR.tsv",
#  "1Ex10_KO_STAT1",
#  "SEACR"
#)

args = commandArgs(trailingOnly = T)
tab_all = args[1]
tab_olap = args[2]
tab_s1 = args[3]
tab_s2 = args[4]
out_sum_tab= args[5]
sample_name = args[6]
tool = args[7]
rep1 = args[8]
rep2 = args[9]
# fdr_cutof = as.numeric(args[10])
# l2fc_cutof= as.numeric(args[11])


names = c("chr", "start", "end", "name", "score", "signal", "real_summit") # rsp = relative summit position
classes = c("character","integer","integer","character","numeric", "numeric", "character")

peaks1 = fread(rep1, sep = "\t", col.names = paste0(names,"_1"), key = c("chr_1","start_1","end_1"), colClasses = classes)
if(peaks1[,.N]==0) {
  peaks1 = data.table(matrix(ncol = length(names), nrow = 0))
  names_1 = paste0(names,'_1')
  colnames(peaks1) = names_1
  setkeyv(peaks1, names_1[1:3])
  for(i in seq_along(names_1)) {
    class(peaks1[[names_1[i]]]) = classes[i]
  }
  peaks1[,summit_1 := integer()]
  peaks1[,len_1 := integer()]
  peaks1[,tpm_1 := numeric()]
} else {
	peaks1[,name_1 := paste0(name_1,"_1")]
	peaks1[,summit_1:={
	  sum_bounds=str_split(sub('^.*:','',real_summit_1),'-')[[1]]
	  (as.numeric(sum_bounds[1])-start_1)+floor((as.numeric(sum_bounds[2])-as.numeric(sum_bounds[1]))/2)
	}, by=seq_along(chr_1)]
	peaks1[,len_1:=end_1-start_1]
	peaks1[,tpm_1:=1000*score_1/len_1]
}

peaks2 = fread(rep2, sep = "\t", col.names = paste0(names,"_2"), key = c("chr_2","start_2","end_2"), colClasses = classes)
if(peaks2[,.N]==0) {
  peaks2 = data.table(matrix(ncol = 10, nrow = 0))
  names_2 = paste0(names,'_2')
  colnames(peaks2) = names_2
  setkeyv(peaks2, names_2[1:3])
  for(i in seq_along(names_2)) {
    class(peaks2[[names_2[i]]]) = classes[i]
  }
  peaks2[,summit_2 := integer()]
  peaks2[,len_2 := integer()]
  peaks2[,tpm_2 := numeric()]
} else {
  peaks2[,name_2 := paste0(name_2,"_2")]
  peaks2[,summit_2:={
    sum_bounds=str_split(sub('^.*:','',real_summit_2),'-')[[1]]
    (as.numeric(sum_bounds[1])-start_2)+floor((as.numeric(sum_bounds[2])-as.numeric(sum_bounds[1]))/2)
  }, by=seq_along(chr_2)]
  peaks2[,len_2:=end_2-start_2]
  peaks2[,tpm_2:=1000*score_2/len_2]
}

# overlap two sets of identified peaks
overlapped1 = foverlaps(peaks1, peaks2, by.x=c("chr_1","start_1","end_1"), by.y=c("chr_2","start_2","end_2"), nomatch = NA)
setnames(overlapped1, "chr_1", "chr")
# write down the singletons
fwrite(overlapped1[
  is.na(name_2), 
  .(chr,start=start_1,end=end_1,name=name_1,score=tpm_1,strand='.',signal=signal_1,pval=0,qval=0,summit=summit_1,
    real_summit=real_summit_1,peak_len=len_1,total_cov=score_1)][order(chr)], 
       tab_s1,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

overlapped2 = foverlaps(peaks2, peaks1, by.y=c("chr_1","start_1","end_1"), by.x=c("chr_2","start_2","end_2"), nomatch = NA)
setnames(overlapped2, "chr_2", "chr")
fwrite(overlapped2[
  is.na(name_1), 
  .(chr,start=start_2,end=end_2,name=name_2,score=tpm_2,strand='.',signal=signal_2,pval=0,qval=0,summit=summit_2,
    real_summit=real_summit_2,peak_len=len_2,total_cov=score_2)][order(chr)], 
       tab_s2,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

merged = rbind(overlapped1[,.(chr,start=start_1,end=end_1,name=name_1,score=tpm_1,strand='.',signal=signal_1,pval=0,qval=0,summit=summit_1,
                              real_summit=real_summit_1,peak_len=len_1,total_cov=score_1)],
               overlapped2[is.na(name_1),
                           .(chr,start=start_2,end=end_2,name=name_2,score=tpm_2,strand='.',signal=signal_2,pval=0,qval=0,summit=summit_2,
                             real_summit=real_summit_2,peak_len=len_2,total_cov=score_2)])
merged = unique(merged)
fwrite(merged[order(chr)], 
       tab_all,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

overlapped = rbind(overlapped1[!is.na(name_1)&!is.na(name_2), 
                               .(chr,start=start_1,end=end_1,name=name_1,score=tpm_1,strand='.',signal=signal_1,pval=0,qval=0,summit=summit_1,
                                 real_summit=real_summit_1,peak_len=len_1,total_cov=score_1)],
                   overlapped2[!is.na(name_1)&!is.na(name_2), 
                               .(chr,start=start_2,end=end_2,name=name_2,score=tpm_2,strand='.',signal=signal_2,pval=0,qval=0,summit=summit_2,
                                 real_summit=real_summit_2,peak_len=len_2,total_cov=score_2)])
overlapped = unique(overlapped)
fwrite(overlapped[order(chr)], 
       tab_olap,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

sum_tab = data.table(sample = sample_name)
sum_tab[,tool := tool]
sum_tab[,rep1_peaks := peaks1[,.N]]
sum_tab[,rep2_peaks := peaks2[,.N]]
sum_tab[,merged_peaks := merged[,.N]]
sum_tab[,overlap_peaks := overlapped1[!is.na(name_1) & !is.na(name_2), .N]]
sum_tab[,rep1_overlap := overlapped1[!is.na(name_1) & !is.na(name_2), length(unique(name_1))]]
sum_tab[,rep2_overlap := overlapped2[!is.na(name_1) & !is.na(name_2), length(unique(name_2))]]
sum_tab[,rep1_only := overlapped1[is.na(name_2), .N]]
sum_tab[,rep2_only := overlapped2[is.na(name_1), .N]]
if(sum_tab$rep1_peaks==0) {sum_tab[,rep1_reprod := 0]
} else {sum_tab[,rep1_reprod := 100-(rep1_only/rep1_peaks)*100]}
if(sum_tab$rep2_peaks==0) {sum_tab[,rep2_reprod := 0]
} else {sum_tab[,rep2_reprod := 100-(rep2_only/rep2_peaks)*100]}
# sum_tab[,eval(paste0(tool,"_sig_up")) := overlapped1[is.na(name_2) & qval_1 > -log10(fdr_cutof), .N]]
# sum_tab[,eval(paste0(tool,"_sig_dn")) := overlapped2[is.na(name_1) & qval_2 > -log10(fdr_cutof), .N]]
# sum_tab[,eval(paste0(tool,"_sig_overlap")) := overlapped[overlap_FDR < fdr_cutof, .N]]
fwrite(sum_tab, out_sum_tab, sep = '\t', row.names = F, col.names = T)
