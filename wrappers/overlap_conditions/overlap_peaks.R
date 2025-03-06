library(data.table)

#setwd("/mnt/nfs/shared/S3acgt/sequia/15869__chipseq_analysis__ChIP_peak_calling__231214")
#c1 = "results/SEACR_peaks/pCSTF2_INH_VS_IgG/pCSTF2_INH_VS_IgG.no_dups.peaks.all.narrowPeak"
#c2 = "results/SEACR_peaks/pCSTF2_VS_IgG/pCSTF2_VS_IgG.no_dups.peaks.all.narrowPeak"
#tab_all = "results/overlapped_peaks/pCSTF2_INH_vs_pCSTF2/overlapped_peaks.no_dups.by_SEACR.bed"
#tab_hist = "results/overlapped_peaks/pCSTF2_INH_vs_pCSTF2/overlapped_peaks.no_dups.by_SEACR.hist.tsv"
#tab_s1 = "results/overlapped_peaks/pCSTF2_INH_vs_pCSTF2/singletons_in_pCSTF2_INH.no_dups.by_SEACR.bed"
#tab_s2 = "results/overlapped_peaks/pCSTF2_INH_vs_pCSTF2/singletons_in_pCSTF2.no_dups.by_SEACR.bed"
#out_sum_tab= "results/overlapped_peaks/pCSTF2_INH_vs_pCSTF2/summary_table.no_dups.by_SEACR.tsv"
#comparison = "pCSTF2_INH_vs_pCSTF2"
#tool = "SEACR"
#fdr_cutof = 0.05
#l2fc_cutof= 0

args = commandArgs(trailingOnly = T)
c1 = args[1]
c2 = args[2]
tab_all = args[3]
tab_hist = args[4]
tab_s1 = args[5]
tab_s2 = args[6]
out_sum_tab= args[7]
comparison = args[8]
tool = args[9]
fdr_cutof = as.numeric(args[10])
l2fc_cutof= as.numeric(args[11])

if(tool == "SEACR") {
  names = c("chr","start","end","name","score","summit_cov","summit_pos")
  classes = c("character","integer","integer","character","numeric","numeric","character")
} else {
  names = c("chr", "start", "end", "name", "score", "strand", "l2fc", "pval", "qval", "rel_summit_pos") # rsp = relative summit position
  classes = c("character","integer","integer","character","numeric","character","numeric","numeric","numeric","numeric")
}

peaks1 = fread(cmd = paste0("cut -f 1-",length(names)," ",c1), sep = "\t", col.names = paste0(names,"_1"), key = c("chr_1","start_1","end_1"), colClasses = classes)
if(peaks1[,.N]==0) {
  peaks1 = data.table(matrix(ncol = length(names), nrow = 0))
  names_1 = paste0(names,'_1')
  colnames(peaks1) = names_1
  setkeyv(peaks1, names_1[1:3])
  for(i in seq_along(names_1)) {
    class(peaks1[[names_1[i]]]) = classes[i]
  }
}
peaks1[,len_1:=end_1-start_1]
peaks2 = fread(cmd = paste0("cut -f 1-",length(names)," ",c2), sep = "\t", col.names = paste0(names,"_2"), key = c("chr_2","start_2","end_2"), colClasses = classes)
if(peaks2[,.N]==0) {
  peaks2 = data.table(matrix(ncol = length(names), nrow = 0))
  names_2 = paste0(names,'_2')
  colnames(peaks2) = names_2
  setkeyv(peaks2, names_2[1:3])
  for(i in seq_along(names_2)) {
    class(peaks2[[names_2[i]]]) = classes[i]
  }
}
peaks2[,len_2:=end_2-start_2]

# overlap two sets of identified peaks
overlapped1 = foverlaps(peaks1, peaks2, by.x=c("chr_1","start_1","end_1"), by.y=c("chr_2","start_2","end_2"), nomatch = NA)
setnames(overlapped1, "chr_1", "chr")
# write down the singletons and remove them from table
if(tool == "SEACR") {
  result = overlapped1[is.na(name_2),
    .(`#chr` = chr,
      start = start_1,
      end = end_1,
      name = name_1,
      score = score_1,
      summit_cov = summit_cov_1,
      summit = summit_pos_1)
  ][order(`#chr`)]
} else {
  result = overlapped1[is.na(name_2),
    .(`#chr` = chr,
      start = start_1,
      end = end_1,
      name = name_1,
      score = score_1,
      strand = strand_1,
      l2fc = l2fc_1,
      pval = pval_1,
      qval = qval_1,
      summit = rel_summit_pos_1)
  ][order(`#chr`)]
}
fwrite(result,
       tab_s1,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

overlapped2 = foverlaps(peaks2, peaks1, by.y=c("chr_1","start_1","end_1"), by.x=c("chr_2","start_2","end_2"), nomatch = NA)
setnames(overlapped2, "chr_2", "chr")
if(tool == "SEACR") {
  result = overlapped2[is.na(name_1),
    .(`#chr` = chr,
      start = start_2,
      end = end_2,
      name = name_2,
      score = score_2,
      summit_cov = summit_cov_2,
      summit = summit_pos_2)
  ][order(`#chr`)]
} else {
  result = overlapped2[is.na(name_1),
    .(`#chr` = chr,
      start = start_2,
      end = end_2,
      name = name_2,
      score = score_2,
      strand = strand_2,
      l2fc = l2fc_2,
      pval = pval_2,
      qval = qval_2,
      summit = rel_summit_pos_2)
  ][order(`#chr`)]
}
fwrite(result,
       tab_s2,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

sum_tab = data.table(comparison = comparison)
sum_tab[,tool := tool]
sum_tab[,merged_peaks := overlapped1[,.N]+overlapped2[is.na(name_1),.N]]
sum_tab[,cond1_unique := overlapped1[is.na(name_2), .N]]
sum_tab[,cond2_unique := overlapped2[is.na(name_1), .N]]
sum_tab[,overlap_peaks := overlapped1[!is.na(name_1) & !is.na(name_2), .N]]
if(tool == "SEACR") {
  sum_tab[,cond1_signif := overlapped1[is.na(name_2), .N]]
  sum_tab[,cond2_signif := overlapped2[is.na(name_1), .N]]
} else {
  sum_tab[,cond1_signif := overlapped1[is.na(name_2) & qval_1 > -log10(fdr_cutof), .N]]
  sum_tab[,cond2_signif := overlapped2[is.na(name_1) & qval_2 > -log10(fdr_cutof), .N]]
}

overlapped = overlapped1[!is.na(name_1)&!is.na(name_2)]
# compute overlap length and its percentage  against the combined length of both peaks divided by 2
overlapped[ ,overlap_len:=fifelse(end_1<end_2,end_1,end_2)-ifelse(start_1>start_2,start_1,start_2) ]
overlapped[ ,overlap_perc:=round(overlap_len/((len_1+len_2)*0.5),4) ]
# overlapped[ ,overlap_pval:=10^-sum(pval_1,pval_2)]
if(tool == "SEACR") {
#  overlapped[ ,overlap_l2fc:=fifelse(l2fc_1<l2fc_2, l2fc_2, l2fc_1)]
#  overlapped[ ,overlap_pval:=0, by=seq_along(chr)]
#  overlapped[ ,overlap_FDR:= 0, by=seq_along(chr)]
#  overlapped[ ,overlap_qval:=0, by=seq_along(chr)]
  overlapped[ ,overlap_score:=ifelse(score_1>score_2, score_1, score_2)]
  # add number of significant overlapped peaks into summary file
  sum_tab[,eval(paste0(tool,"_sig_overlap")) := overlapped[,.N]]
} else {
  overlapped[ ,overlap_l2fc:=fifelse(l2fc_1<l2fc_2, l2fc_2, l2fc_1)]
  overlapped[ ,overlap_pval:=pchisq(-2*sum(log(10^-c(pval_1,pval_2))), 4, lower.tail=FALSE), by=seq_along(chr)]
  overlapped[ ,overlap_FDR:=p.adjust(overlap_pval, method = "fdr")]
  overlapped[ ,overlap_qval:=-log10(overlap_FDR)]
  overlapped[ ,overlap_score:=floor(10*overlap_qval)]
  # add number of significant overlapped peaks into summary file
  sum_tab[,eval(paste0(tool,"_sig_overlap")) := overlapped[overlap_FDR < fdr_cutof, .N]]
}
fwrite(sum_tab, out_sum_tab, sep = '\t', row.names = F, col.names = T)

fwrite(overlapped[order(chr)], 
       tab_all,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

# sum up a histo table with count of overlapped peaks based on the minimum overlap length in percentage of both peaks or one of the peaks respectively
histo = data.table(perc=c(0:10)*10)
histo[,count_overlap_perc:=nrow(overlapped[overlap_perc*100>=perc]), by=perc]
histo[,count_over_p1:=nrow(overlapped[overlap_len/len_1*100>=perc]), by=perc]
histo[,count_over_p2:=nrow(overlapped[overlap_len/len_2*100>=perc]), by=perc]
histo[perc == 0] = list(perc=0,count_overlap_perc=nrow(overlapped),count_over_p1=nrow(peaks1),count_over_p2=nrow(peaks2))
# histo

fwrite(histo,
       tab_hist,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)
