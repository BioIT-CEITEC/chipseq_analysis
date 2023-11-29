#########################################
# wrapper for rule: call_DiffBind
#########################################
shell = function(cmd) {
  cat(system(cmd, intern = T), sep = '\n')
}

logfile = snakemake@log[["run"]]
sink(logfile, append = T, type = "output")
sink(stdout(), append = T, type = "message")

cat("##\n## RULE: call_DiffBind \n##\n")
cat("## CONDA:\n")
shell("conda list 2>&1")

library(data.table)
library(DiffBind)
# library(ChIPQC)
# library(ggplot2)
# library(DESeq2)
# library(edgeR)
options("browser"="false")

## Doesn't work very well
# do_plot = function(toPlot, toName) {
#   png(filename = gsub(".pdf$",".png",toName), width = 1080, height = 1080, pointsize = 20)
#   toPlot
#   dev.off()
#   ggsave(filename = gsub(".png$",".pdf",toName), plot = toPlot, width = 1080, height = 1080, units = "px", pointsize = 15, device = "pdf")
# }



# setwd("/mnt/ssd/ssd_1/sequia/5232__chipseq_analysis__peak_calling_test__220810")
# design_tab = "results/DiffBind/p53R273H_360_vs_p53R273H_393/design_table.no_dups.tsv"

design_tab = snakemake@input[['tab']]
out_deseq_tab = snakemake@output[['deseq_tab']]
out_edger_tab = snakemake@output[['edger_tab']]
out_sum_tab = snakemake@output[['sum_tab']]
info_tab = snakemake@output[['info_tab']]
cor_heat_1 = snakemake@output[['cor_heat_1']]
cor_heat_2 = snakemake@output[['cor_heat_2']]
cor_heat_deseq = snakemake@output[['cor_heat_deseq']]
cor_heat_edger = snakemake@output[['cor_heat_edger']]
pca_plot_1 = snakemake@output[['pca_plot_1']]
venn_plot = snakemake@output[['venn_plot']]
PCA_plot_deseq = snakemake@output[['pca_plot_deseq']]
PCA_plot_edger = snakemake@output[['pca_plot_edger']]
MA_plot_deseq = snakemake@output[['ma_plot_deseq']]
MA_plot_edger = snakemake@output[['ma_plot_edger']]
XY_plot_deseq = snakemake@output[['xy_plot_deseq']]
XY_plot_edger = snakemake@output[['xy_plot_edger']]
volcano_plot_deseq = snakemake@output[['volcano_plot_deseq']]
volcano_plot_edger = snakemake@output[['volcano_plot_edger']]
DBA_heat_deseq = snakemake@output[['DBA_heat_deseq']]
DBA_heat_edger = snakemake@output[['DBA_heat_edger']]
fdr_cutof = snakemake@params[['fdr_cutof']]
l2fc_cutof= snakemake@params[['l2fc_cutof']]
comparison= snakemake@params[['comparison']]
top_num = snakemake@params[['top']]

cat("# reading design samplesheet ",design_tab," as DBA object\n")
design_tab = read.csv(design_tab, sep = '\t')
print(design_tab)
data = dba(sampleSheet = design_tab)
# do_plot(dba.plotHeatmap(data, margin = 20, main = "Correlation heatmap (using peak scores)"), "test_plot.png")
png(filename = gsub(".pdf$",".png",cor_heat_1), width = 1080, height = 1080, pointsize = 20)
dba.plotHeatmap(data, margin = 20, main = "Correlation heatmap (using peak scores)")
dev.off()
pdf(file = gsub(".png$",".pdf",cor_heat_1), width = 12, height = 12, pointsize = 15)
dba.plotHeatmap(data, margin = 20, main = "Correlation heatmap (using peak scores)")
dev.off()

cat('# re-counting coverage and summits\n')
data = dba.count(data)
png(filename = gsub(".pdf$",".png",cor_heat_2), width = 1080, height = 1080, pointsize = 20)
dba.plotHeatmap(data, margin = 20, main = "Correlation heatmap (using read counts)")
dev.off()
pdf(file = gsub(".png$",".pdf",cor_heat_2), width = 12, height = 12, pointsize = 15)
dba.plotHeatmap(data, margin = 20, main = "Correlation heatmap (using read counts)")
dev.off()
info = dba.show(data)
print(info)
if("Reads" %in% colnames(info) && "FRiP" %in% colnames(info)) info = cbind(info, PeakReads=round(info$Reads * info$FRiP))

cat('# normalizing data\n')
data = dba.normalize(data, method = DBA_ALL_METHODS)
norm = dba.normalize(data, method = DBA_ALL_METHODS, bRetrieve = T)
info = cbind(info, DESeqNormFactor=norm$DESeq2$norm.factors, DESeqNormReads=round(info$Reads/norm$DESeq2$norm.factors))
info = cbind(info, edgeRNormFactor=norm$edgeR$norm.factors, edgeRNormReads=round(info$Reads/norm$edgeR$norm.factors))
fwrite(info, info_tab, sep = '\t', row.names = F, col.names = T)
data = dba.contrast(data, minMembers = 2)

cat('# doing differential analysis\n')
data = dba.analyze(data, method = DBA_ALL_METHODS, bBlacklist = F, bGreylist = F)
# dba.analyze(data, method = DBA_ALL_METHODS, bRetrieveAnalysis = T)

edger_tab = as.data.table(dba.report(data, method = DBA_EDGER, th = 1, bCounts = T, DataType = DBA_DATA_FRAME))
fwrite(edger_tab,
       out_edger_tab,
       sep = '\t',
       row.names = F,
       col.names = T,
       quote = F)
deseq_tab = as.data.table(dba.report(data, method = DBA_DESEQ2, th = 1, bCounts = T, DataType = DBA_DATA_FRAME))
fwrite(deseq_tab,
       out_deseq_tab,
       sep = '\t',
       row.names = F,
       col.names = T,
       quote = F)

sum_tab = data.table(
  comparison = comparison,
  total = deseq_tab[,.N],
  deseq_sig= deseq_tab[FDR < fdr_cutof,.N],
  deseq_up = deseq_tab[FDR < fdr_cutof & Fold > l2fc_cutof,.N],
  deseq_dn = deseq_tab[FDR < fdr_cutof & Fold < -l2fc_cutof,.N],
  edger_sig= edger_tab[FDR < fdr_cutof,.N],
  edger_up = edger_tab[FDR < fdr_cutof & Fold > l2fc_cutof,.N],
  edger_dn = edger_tab[FDR < fdr_cutof & Fold < -l2fc_cutof,.N]
)
fwrite(sum_tab, out_sum_tab, sep = '\t', row.names = F, col.names = T)

if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.edgeR']]) > 1) {
  png(filename = gsub(".pdf$",".png",cor_heat_edger), width = 1080, height = 1080, pointsize = 20)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, method = DBA_EDGER, margin = 20, main = "Correlation heatmap (using diff. bound sites)")
  dev.off()
  pdf(file = gsub(".png$",".pdf",cor_heat_edger), width = 12, height = 12, pointsize = 15)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, method = DBA_EDGER, margin = 20, main = "Correlation heatmap (using diff. bound sites)")
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",cor_heat_edger), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("Not enough results from edgeR with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",cor_heat_edger), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("Not enough results from edgeR with FDR below ",fdr_cutof))
  dev.off()
}
if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.DESeq2']]) > 1) {
  png(filename = gsub(".pdf$",".png",cor_heat_deseq), width = 1080, height = 1080, pointsize = 20)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, method = DBA_DESEQ2, margin = 20, main = "Correlation heatmap (using diff. bound sites)")
  dev.off()
  pdf(file = gsub(".png$",".pdf",cor_heat_deseq), width = 12, height = 12, pointsize = 15)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, method = DBA_DESEQ2, margin = 20, main = "Correlation heatmap (using diff. bound sites)")
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",cor_heat_deseq), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("Not enough results from DESeq2 with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",cor_heat_deseq), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("Not enough results from DESeq2 with FDR below ",fdr_cutof))
  dev.off()
}

if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.edgeR']]) > 0 && as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.DESeq2']]) > 0) {
  png(filename = gsub(".pdf$",".png",venn_plot), width = 1080, height = 1080, pointsize = 20)
  dba.plotVenn(data, contrast = 1, bDB = T, method = DBA_ALL_METHODS, th = fdr_cutof)
  dev.off()
  pdf(file = gsub(".png$",".pdf",venn_plot), width = 12, height = 12, pointsize = 15)
  dba.plotVenn(data, contrast = 1, bDB = T, method = DBA_ALL_METHODS, th = fdr_cutof)
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",venn_plot), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("No results from edgeR or DESeq2 with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",venn_plot), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("No results from edgeR or DESeq2 with FDR below ",fdr_cutof))
  dev.off()
}

png(filename = gsub(".pdf$",".png",pca_plot_1), width = 1080, height = 1080, pointsize = 20)
dba.plotPCA(data, DBA_CONDITION, label = DBA_REPLICATE)
dev.off()
pdf(file = gsub(".png$",".pdf",pca_plot_1), width = 12, height = 12, pointsize = 15)
dba.plotPCA(data, DBA_CONDITION, label = DBA_REPLICATE)
dev.off()

if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.edgeR']]) >= nrow(design_tab)) {
  png(filename = gsub(".pdf$",".png",PCA_plot_edger), width = 1080, height = 1080, pointsize = 20)
  dba.plotPCA(data, contrast = 1, th=fdr_cutof, label = DBA_REPLICATE, method = DBA_EDGER)
  dev.off()
  pdf(file = gsub(".png$",".pdf",PCA_plot_edger), width = 12, height = 12, pointsize = 15)
  dba.plotPCA(data, contrast = 1, th=fdr_cutof, label = DBA_REPLICATE, method = DBA_EDGER)
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",PCA_plot_edger), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("Too few results from edgeR with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",PCA_plot_edger), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("Too few results from edgeR with FDR below ",fdr_cutof))
  dev.off()
}
if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.DESeq2']]) >= nrow(design_tab)) {
  png(filename = gsub(".pdf$",".png",PCA_plot_deseq), width = 1080, height = 1080, pointsize = 20)
  dba.plotPCA(data, contrast = 1, th=fdr_cutof, label = DBA_REPLICATE, method = DBA_DESEQ2)
  dev.off()
  pdf(file = gsub(".png$",".pdf",PCA_plot_deseq), width = 12, height = 12, pointsize = 15)
  dba.plotPCA(data, contrast = 1, th=fdr_cutof, label = DBA_REPLICATE, method = DBA_DESEQ2)
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",PCA_plot_deseq), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("Too few results from DESeq2 with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",PCA_plot_deseq), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("Too few results from DESeq2 with FDR below ",fdr_cutof))
  dev.off()
}

png(filename = gsub(".pdf$",".png",MA_plot_deseq), width = 1080, height = 1080, pointsize = 20)
dba.plotMA(data, contrast = 1, method = DBA_DESEQ2, th = fdr_cutof)
dev.off()
pdf(file = gsub(".png$",".pdf",MA_plot_deseq), width = 12, height = 12, pointsize = 15)
dba.plotMA(data, contrast = 1, method = DBA_DESEQ2, th = fdr_cutof)
dev.off()
png(filename = gsub(".pdf$",".png",MA_plot_edger), width = 1080, height = 1080, pointsize = 20)
dba.plotMA(data, contrast = 1, method = DBA_EDGER, th = fdr_cutof)
dev.off()
pdf(file = gsub(".png$",".pdf",MA_plot_edger), width = 12, height = 12, pointsize = 15)
dba.plotMA(data, contrast = 1, method = DBA_EDGER, th = fdr_cutof)
dev.off()

png(filename = gsub(".pdf$",".png",XY_plot_deseq), width = 1080, height = 1080, pointsize = 20)
dba.plotMA(data, contrast = 1, method = DBA_DESEQ2, th = fdr_cutof, bXY = T)
dev.off()
pdf(file = gsub(".png$",".pdf",XY_plot_deseq), width = 12, height = 12, pointsize = 15)
dba.plotMA(data, contrast = 1, method = DBA_DESEQ2, th = fdr_cutof, bXY = T)
dev.off()
png(filename = gsub(".pdf$",".png",XY_plot_edger), width = 1080, height = 1080, pointsize = 20)
dba.plotMA(data, contrast = 1, method = DBA_EDGER, th = fdr_cutof, bXY = T)
dev.off()
pdf(file = gsub(".png$",".pdf",XY_plot_edger), width = 12, height = 12, pointsize = 15)
dba.plotMA(data, contrast = 1, method = DBA_EDGER, th = fdr_cutof, bXY = T)
dev.off()

png(filename = gsub(".pdf$",".png",volcano_plot_deseq), width = 1080, height = 1080, pointsize = 20)
dba.plotVolcano(data, contrast = 1, method = DBA_DESEQ2, th = fdr_cutof, bLabels = T, maxLabels = top_num)
dev.off()
pdf(file = gsub(".png$",".pdf",volcano_plot_deseq), width = 12, height = 12, pointsize = 15)
dba.plotVolcano(data, contrast = 1, method = DBA_DESEQ2, th = fdr_cutof, bLabels = T, maxLabels = top_num)
dev.off()
png(filename = gsub(".pdf$",".png",volcano_plot_edger), width = 1080, height = 1080, pointsize = 20)
dba.plotVolcano(data, contrast = 1, method = DBA_EDGER, th = fdr_cutof, bLabels = T, maxLabels = top_num)
dev.off()
pdf(file = gsub(".png$",".pdf",volcano_plot_edger), width = 12, height = 12, pointsize = 15)
dba.plotVolcano(data, contrast = 1, method = DBA_EDGER, th = fdr_cutof, bLabels = T, maxLabels = top_num)
dev.off()

hmap <- colorRampPalette(c("red", "black", "green"))(n = 13)
if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.edgeR']]) > 1) {
  png(filename = gsub(".pdf$",".png",DBA_heat_edger), width = 1080, height = 1080, pointsize = 20)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, correlations = F, maxSites = top_num, method = DBA_EDGER, colScheme = hmap, margin = 20)
  dev.off()
  pdf(file = gsub(".png$",".pdf",DBA_heat_edger), width = 12, height = 12, pointsize = 15)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, correlations = F, maxSites = top_num, method = DBA_EDGER, colScheme = hmap, margin = 20)
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",DBA_heat_edger), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("Not enough results from edgeR with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",DBA_heat_edger), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("Not enough results from edgeR with FDR below ",fdr_cutof))
  dev.off()
}
if(as.numeric(dba.show(data, bContrasts=TRUE, th=fdr_cutof)[['DB.DESeq2']]) > 1) {
  png(filename = gsub(".pdf$",".png",DBA_heat_deseq), width = 1080, height = 1080, pointsize = 20)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, correlations = F, maxSites = top_num, method = DBA_DESEQ2, colScheme = hmap, margin = 20)
  dev.off()
  pdf(file = gsub(".png$",".pdf",DBA_heat_deseq), width = 12, height = 12, pointsize = 15)
  dba.plotHeatmap(data, contrast = 1, th=fdr_cutof, correlations = F, maxSites = top_num, method = DBA_DESEQ2, colScheme = hmap, margin = 20)
  dev.off()
} else {
  png(filename = gsub(".pdf$",".png",DBA_heat_deseq), width = 1080, height = 1080, pointsize = 20)
  plot.new()
  title(main=paste0("Not enough results from DESeq2 with FDR below ",fdr_cutof))
  dev.off()
  pdf(file = gsub(".png$",".pdf",DBA_heat_deseq), width = 12, height = 12, pointsize = 15)
  plot.new()
  title(main=paste0("Not enough results from DESeq2 with FDR below ",fdr_cutof))
  dev.off()
}