library(data.table)

#setwd("/mnt/nfs/shared/S3acgt/sequia/15869__chipseq_analysis__ChIP_peak_calling__231214")
#c1 = "results/MACS_peaks/H7_M_VS_Ig_M/H7_M_VS_Ig_M.no_dups.peaks.all.narrowPeak"
#c2 = "results/MACS_peaks/BRD4_R_VS_Ig_R/BRD4_R_VS_Ig_R.no_dups.peaks.all.narrowPeak"

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


names = c("chr", "start", "end", "name", "score", "strand", "l2fc", "pval", "qval", "rel_summit_pos") # rsp = relative summit position
classes = c("character","integer","integer","character","numeric","character","numeric","numeric","numeric","numeric")

peaks1 = fread(c1, sep = "\t", col.names = paste0(names,"_1"), key = c("chr_1","start_1","end_1"), colClasses = classes)
if(peaks1[,.N]==0) {
  peaks1 = data.table(matrix(ncol = 10, nrow = 0))
  names_1 = paste0(names,'_1')
  colnames(peaks1) = names_1
  setkeyv(peaks1, names_1[1:3])
  for(i in seq_along(names_1)) {
    class(peaks1[[names_1[i]]]) = classes[i]
  }
}
peaks1[,len_1:=end_1-start_1]
peaks2 = fread(c2, sep = "\t", col.names = paste0(names,"_2"), key = c("chr_2","start_2","end_2"), colClasses = classes)
if(peaks2[,.N]==0) {
  peaks2 = data.table(matrix(ncol = 10, nrow = 0))
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
fwrite(overlapped1[is.na(name_2), .(chr,start=start_1,end=end_1,name=name_1,score=score_1,strand=strand_1,l2fc=l2fc_1,pval=pval_1,qval=qval_1,summit=rel_summit_pos_1)][order(chr)], 
       tab_s1,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

overlapped2 = foverlaps(peaks2, peaks1, by.y=c("chr_1","start_1","end_1"), by.x=c("chr_2","start_2","end_2"), nomatch = NA)
setnames(overlapped2, "chr_2", "chr")
fwrite(overlapped2[is.na(name_1), .(chr,start=start_2,end=end_2,name=name_2,score=score_2,strand=strand_2,l2fc=l2fc_2,pval=pval_2,qval=qval_2,summit=rel_summit_pos_2)][order(chr)], 
       tab_s2,
       col.names = T,
       row.names = F,
       sep = "\t",
       quote = F)

sum_tab = data.table(comparison = comparison)
sum_tab[,eval(paste0(tool,"_total")) := overlapped1[,.N]+overlapped2[is.na(name_1),.N]]
sum_tab[,eval(paste0(tool,"_up")) := overlapped1[is.na(name_2), .N]]
sum_tab[,eval(paste0(tool,"_dn")) := overlapped2[is.na(name_1), .N]]
sum_tab[,eval(paste0(tool,"_overlap")) := overlapped1[!is.na(name_1) & !is.na(name_2), .N]]
sum_tab[,eval(paste0(tool,"_sig_up")) := overlapped1[is.na(name_2) & qval_1 > -log10(fdr_cutof), .N]]
sum_tab[,eval(paste0(tool,"_sig_dn")) := overlapped2[is.na(name_1) & qval_2 > -log10(fdr_cutof), .N]]

overlapped = overlapped1[!is.na(name_1)&!is.na(name_2)]
# compute overlap length and its percentage  against the combined length of both peaks divided by 2
overlapped[ ,overlap_len:=ifelse(end_1<end_2,end_1,end_2)-ifelse(start_1>start_2,start_1,start_2) ]
overlapped[ ,overlap_perc:=round(overlap_len/((len_1+len_2)*0.5),4) ]
overlapped[ ,overlap_l2fc:=ifelse(l2fc_1<l2fc_2, l2fc_2, l2fc_1)]
# overlapped[ ,overlap_pval:=10^-sum(pval_1,pval_2)]
overlapped[ ,overlap_pval:=pchisq(-2*sum(log(10^-c(pval_1,pval_2))), 4, lower.tail=FALSE), by=seq_along(chr)]
overlapped[ ,overlap_FDR:=p.adjust(overlap_pval, method = "fdr")]
overlapped[ ,overlap_qval:=-log10(overlap_FDR)]
overlapped[ ,overlap_score:=floor(10*overlap_qval)]
# add number of sginificant overlapped peaks into summary file and write it down
sum_tab[,eval(paste0(tool,"_sig_overlap")) := overlapped[overlap_FDR < fdr_cutof, .N]]
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
