library(data.table)
library(rtracklayer)
library(stringr)

args <- commandArgs(trailingOnly = T)
bed_file <- args[1]
gtf_file <- args[2]
out_bed  <- args[3]
feat_type <- args[4]
annotate_by<- str_split(args[5], ",")[[1]]  # must be comma separated list of valid GTF TAGS (e.g., gene_name, gene_id)

# bed_file <- "/mnt/ssd/ssd_1/snakemake/stage233_Dragana_Clip-Seq_4lib/CLIP-seq_general_analysis/macs2/Adar1_CLASH_rep2/Adar1_CLASH_rep2.uniq_reads.no_dups.peaks.narrowPeak"
# # bed_file <- "/mnt/ssd/ssd_1/snakemake/stage233_Dragana_Clip-Seq_4lib/CLIP-seq_general_analysis/pureClip/Adar1_CLASH_rep2/Adar1_CLASH_rep2.uniq_reads.no_dups.binding_regions.bed"
# gtf_file <- "/mnt/ssd/ssd_3/references/homsap/GRCh38-p10/annot/GRCh38-p10.gtf"
# feat_type <- "gene"
# annotate_by<- c("gene_name")

comp_tab <- fread(cmd=paste0("cut -f 1-3 ",bed_file), sep = "\t", col.names = c("chr","start","end"), colClasses = c("character","numeric","numeric"), key = c("chr","start","end"))
comp_tab[,start:=start+1]

ref <- as.data.table(rtracklayer::import(gtf_file))[type == feat_type, c("seqnames","start","end",annotate_by), with=F]
setkeyv(ref, c("seqnames","start","end"))

#overlapped <- comp_tab[ref,on=c("chr==seqnames"),allow.cartesian=T][start+1<=i.end & end>=i.start]
overlapped <- foverlaps(comp_tab, ref, by.x=c("chr","start","end"), by.y=c("seqnames","start","end"), nomatch = 0)
overlapped[,c("start","end","i.start","i.end") := list(i.start, i.end, NULL, NULL)]
res <- overlapped[, paste0(unique(get(annotate_by[1])),collapse = ","), by=.(chr,start,end)]
setnames(res, "V1", annotate_by[1])
for(idf in annotate_by[-1]) {
  res <- merge(res,overlapped[, paste0(unique(get(idf)),collapse = ","), by=.(chr,start,end)], by=c("chr","start","end"))
  setnames(res, "V1", idf)
}

bed_file <- fread(bed_file, sep = "\t")
setnames(bed_file, 1:3, c("chr","start","end"))
bed_file$chr <- as.character(bed_file$chr)
res[,start:=start-1]
res <- merge(bed_file, res, by.y = c("chr","start","end"), by.x = c("chr","start","end"), all.x = T)
res[is.na(res)] <- ""

fwrite(res, out_bed, sep="\t", row.names = F, col.names = T, quote = F)
