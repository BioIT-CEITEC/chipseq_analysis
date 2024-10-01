library(data.table)
library(rjson)

# setwd("/mnt/nfs/shared/S3acgt/sequia/16808__chipseq_analysis__peak_calling_DMSO_vs_IBU_test_new_alignment__240116")
# mtx <- "results/peaks_QC/peak_profiles/over_all_genes/POL_II_IBU_VS_control.no_dups.from_MACS.all.over_all_genes.mtx.tmp.gz"
# bed <- "gene_sets/all_genes.gene_set.bed"
# out <- "results/peaks_QC/peak_profiles/over_all_genes/POL_II_IBU_VS_control.no_dups.from_MACS.all.over_all_genes.mtx.gz"
#out <- "test.mtx.gz"

args = commandArgs(trailingOnly=TRUE)

mtx <- args[1]
bed <- args[2]
out <- args[3]

head <- suppressMessages(fread(cmd=paste0("zcat ", mtx, "|head -n 1"), sep=NULL, header = F, quote=""))
data <- suppressMessages(fread(cmd=paste0("zcat ", mtx, "|tail -n +2")))
setkeyv(data, c("V1","V2","V3","V4","V5","V6"))
ccols = c("V1","V4","V5","V6")
data[,(ccols):= lapply(.SD, as.character), .SDcols = ccols]
bed <- suppressMessages(fread(cmd=paste0("cat ", bed, "|cut -f 1-6"), header=F, key=c("V1","V2","V3","V4","V5","V6")))
bed[,(ccols):= lapply(.SD, as.character), .SDcols = ccols]

#res <- merge(data, bed, by = c("V1","V2","V3","V4","V5","V6"), all = T)
#res <- res[order(bed)]
res <- data[bed, on=.(V1,V2,V3,V4,V5,V6)]
res[is.na(res)] <- 0

# head <- fromJSON(sub("^@","",head[[1]]))
# head$group_boundaries <- c(0, nrow(bed))
# fwrite(as.list(paste0("@",toJSON(head))), out, col.names = F, sep="\t", quote = F, compress = "gzip")

fwrite(head[,.(sub("(\"group_boundaries\":\\[0,)[0-9]+(\\])",paste0("\\1",nrow(bed),"\\2"),V1))], out, col.names = F, sep="\t", quote = F, compress = "gzip")
fwrite(res,  out, append = T, col.names = F, sep="\t", quote = F, compress = "gzip")

