library(data.table)

args = commandArgs(trailingOnly=TRUE)

mtx <- args[1]
out <- args[2]
smooth <- as.numeric(args[3])

#setwd("/mnt/ssd/ssd_3/temp/martin/Vanacova_chip-seq_playground/")
#mtx <- "ChiP_BY_WT/ChiP_BY_WT.MQ_over_20_treat_pileup.non_overlap_genes.b200_a200.mtx.gz"
#out <- sub(".mtx.gz",".rel_counts.mtx.gz",mtx)
#smooth <- 1

cat(mtx)
head <- suppressMessages(fread(paste0("zcat ",mtx,"|head -1"), sep=NULL, header = F, quote=""))
info <- suppressMessages(fread(paste0("zcat ",mtx,"|tail -n +2|cut -f 1-6")))
data <- suppressMessages(fread(paste0("zcat ",mtx,"|tail -n +2|cut -f 7-")))
for (j in names(data)) set(data,which(is.na(data[[j]])),j,0)
data <- data.table(sapply(data, function(x) {(x+smooth)/mean(x+smooth)}))
fwrite(head, out, col.names = F, sep="\t", quote = F)
fwrite(cbind(info,data), out, append = T, col.names = F, sep="\t", quote = F)

