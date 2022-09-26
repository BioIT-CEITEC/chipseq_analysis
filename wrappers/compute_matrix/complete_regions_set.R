library(data.table)
library(rjson)

#setwd("/mnt/ssd/ssd_1/snakemake/")
#mtx <- "stage241_Smida_ChiP-Seq_pilot.RTXresistant/ChIP-seq_profile_investigation/peaks_profile/over_all_genes/RTXresistant_TFEB/pooled/coverage_matrix.no_dups.mtx"
#bed  <- "stage241_Smida_ChiP-Seq_pilot.RTXresistant/ChIP-seq_profile_investigation/gene_sets/all_genes.gene_set.bed12"
#out <- "stage241_Smida_ChiP-Seq_pilot.RTXresistant/ChIP-seq_profile_investigation/peaks_profile/over_all_genes/RTXresistant_TFEB/pooled/coverage_matrix.no_dups.mtx"
#out <- "test.mtx.gz"

args = commandArgs(trailingOnly=TRUE)

mtx <- args[1]
bed <- args[2]
out <- args[3]

head <- suppressMessages(fread(cmd=paste0("zcat ", mtx, "|head -n 1"), sep=NULL, header = F, quote=""))
data <- suppressMessages(fread(cmd=paste0("zcat ", mtx, "|tail -n +2"), key=c("V1","V2","V3","V4","V5","V6")))

bed <- suppressMessages(fread(cmd=paste0("cat ", bed, "|cut -f 1-6"), header = F, key=c("V1","V2","V3","V4","V5","V6")))

#res <- merge(data, bed, by = c("V1","V2","V3","V4","V5","V6"), all = T)
#res <- res[order(bed)]
res <- data[bed]
res[is.na(res)] <- 0

# head <- fromJSON(sub("^@","",head[[1]]))
# head$group_boundaries <- c(0, nrow(bed))
# fwrite(as.list(paste0("@",toJSON(head))), out, col.names = F, sep="\t", quote = F, compress = "gzip")

fwrite(head[,.(sub("(\"group_boundaries\":\\[0,)[0-9]+(\\])",paste0("\\1",nrow(bed),"\\2"),V1))], out, col.names = F, sep="\t", quote = F, compress = "gzip")
fwrite(res,  out, append = T, col.names = F, sep="\t", quote = F, compress = "gzip")

