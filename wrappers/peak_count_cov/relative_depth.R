library(data.table)

# setwd("/mnt/ssd/ssd_1/snakemake")
# ref_file <- "stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_BY_WT_vs_ChiP_Drai1/peaks_from_reps.no_dups/merged_peaks.saf"
# depth_file <- "stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_BY_WT_vs_ChiP_Drai1/peaks_from_reps.no_dups/counts/ChiP_BY_WT_rep1.peaks_cov.bdg"
# out_file <- "stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/featureCounts/ChiP_BY_WT_rep1.ChiP_BY_WT_vs_ChiP_Drai1.reps.no_dups.featureCounts.tsv"

args <- commandArgs(trailingOnly = T)

ref_file <- args[1]
depth_file <- args[2]
out_file <- args[3]
metric <- args[4]

ref_file <- fread(ref_file, sep="\t")
bed_file <- ref_file[, .(chr=V1,start=as.numeric(V2),end=as.numeric(V3))]
setkeyv(bed_file,c("chr","start","end"))

depth_file <- fread(depth_file, sep="\t", col.names = c("chr","start","depth"), colClasses = c("character","numeric","numeric"))
depth_file[,end:=start]
setkeyv(depth_file, c("chr","start","end"))

overlapped <- foverlaps(depth_file, bed_file, by.x=c("chr","start","end"), by.y=c("chr","start","end"), nomatch = 0)[,.(max_depth=max(depth), sum_depth=sum(depth), mean_depth=mean(depth), median_depth=median(depth)),by=.(chr,start,end)]
overlapped[,depth_ratio := max_depth/mean_depth]

out <- ref_file[,.(Chr=V1, Start=V2, End=V3, Geneid=V5, Strand=V4, Length=V3-V2+1)]

if(metric == "ratio") {
  out <- out[overlapped[,.SD,.SDcols=c(1:3,8)], on=.(Chr=chr, Start=start, End=end), nomatch = 0]
} else {
  out <- out[overlapped[,.SD,.SDcols=c(1:3,4)], on=.(Chr=chr, Start=start, End=end), nomatch = 0]
}

fwrite(out, out_file, quote = F, sep = "\t", row.names = F, col.names = T)
