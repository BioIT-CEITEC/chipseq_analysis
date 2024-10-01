library(data.table)

fread_vector_of_files <- function(file_list,regex = NULL,add_column = "sample"){
  rbindlist(lapply(file_list,function(x){
    res <- fread(x, sep="\t", col.names = c("Chr","Start","End","Geneid","Strand","Length","Count"))
    if(is.null(regex)){
      res[,(add_column) := x]
    } else {
      res[,(add_column) := gsub(regex,"\\1",x)]
    }
    
  }))
}

run_all <- function(file_list, output_file){
  res_tab <- fread_vector_of_files(file_list, ".*\\/([^\\.]*).*\\.counts\\..*_based.tsv$")
  res_tab[,sample := make.names(sample)]
  res_tab <- dcast.data.table(res_tab, Chr + Start + End + Geneid + Strand + Length ~ sample, value.var = "Count")
  write.table(res_tab, file = output_file, quote = F, row.names = F, col.names = T, sep = "\t")
}

# develop and test 2
# setwd('/mnt/ssd/ssd_1/snakemake/')
# file_list <- c("stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_Dtrf4_vs_ChiP_BY_WT/peaks_from_pooled.no_dups/counts/ChiP_BY_WT_rep1.counts.max_based.tsv","stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_Dtrf4_vs_ChiP_BY_WT/peaks_from_pooled.no_dups/counts/ChiP_Dtrf4_rep1.counts.max_based.tsv","stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_Dtrf4_vs_ChiP_BY_WT/peaks_from_pooled.no_dups/counts/ChiP_BY_WT_rep2.counts.max_based.tsv","stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_Dtrf4_vs_ChiP_BY_WT/peaks_from_pooled.no_dups/counts/ChiP_Dtrf4_rep2.counts.max_based.tsv")
# output_file <- "stage239_Vanacova_ChIP-seq_Pol2.run1/ChIP-seq_profile_investigation/diff_peak_calling/ChiP_Dtrf4_vs_ChiP_BY_WT/peaks_from_pooled.no_dups/counts/merged_counts.max_based.tsv"

# run as Rscript
args <- commandArgs(trailingOnly = T)
file_list <- tail(args,-1)
output_file <- args[1]
run_all(file_list,output_file)