#########################################
# wrapper for rule: plot_profile_and_heatmap
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: plot_profile_and_heatmap \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# ## Add part with combining mtx files when there are more of them using deeptools computeMatrixOperations cbind
# if len(snakemake.input.mtx) > 1:
#   input_mtx = snakemake.params.mtx
#   
#   command = "computeMatrixOperations cbind -m "+" ".join(snakemake.input.mtx)+" -o "+input_mtx+" >> "+snakemake.log.run+" 2>&1"
#   f = open(snakemake.log.run, 'at')
#   f.write("## COMMAND: "+command+"\n")
#   f.close()
#   shell(command)
# else:
#   input_mtx = snakemake.input.mtx[0]

input_mtx = snakemake.input.mtx
title = snakemake.params.title+" ("+snakemake.wildcards.filt+")"

## PLOT ONLY HEATMAP
# plotHeatmap -m ${i} -o ${i%.mtx.gz}.heatmap.pdf --plotTitle "$(basename $i)" --dpi 300 --heatmapWidth 22 --startLabel start --endLabel end --whatToShow "heatmap and colorbar"
command = "plotHeatmap -m "+input_mtx+" -o "+snakemake.output.heatmap+" --outFileNameMatrix "+snakemake.params.heatmtx+" --regionsLabel "+snakemake.wildcards.gs+" --plotTitle \""+title+"\" --dpi "+str(snakemake.params.dpi)+" --heatmapWidth 22 --startLabel start --endLabel end --whatToShow \"heatmap and colorbar\" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

## PLOT ONLY PROFILE
# plotProfile -m ${i} -o ${i%.mtx.gz}.profile.pdf --plotTitle "$(basename $i)" --dpi 300 --averageType mean --plotType std --plotHeight 14 --plotWidth 22 --startLabel start --endLabel end
command = "plotProfile -m "+input_mtx+" -o "+snakemake.output.profile+" --outFileNameData "+snakemake.params.data+" --regionsLabel "+snakemake.wildcards.gs+" --perGroup --plotTitle \""+title+"\" --dpi "+str(snakemake.params.dpi)+" --averageType "+snakemake.params.profile_avrg+" --plotType "+snakemake.params.profile_type+" --plotHeight 14 --plotWidth 22 --startLabel start --endLabel end >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

## PLOT HEATMAP AND PROFILE
# plotHeatmap -m ${i} -o ${i%.mtx.gz}.combined.pdf --plotTitle "$(basename $i)" --dpi 300 --averageType mean --plotType std --heatmapWidth 12 --startLabel start --endLabel end
command = "plotHeatmap -m "+input_mtx+" -o "+snakemake.output.combined+" --regionsLabel "+snakemake.wildcards.gs+" --plotTitle \""+title+"\" --dpi "+str(snakemake.params.dpi)+" --averageTypeSummaryPlot "+snakemake.params.profile_avrg+" --plotType "+snakemake.params.profile_type+" --heatmapWidth 12 --startLabel start --endLabel end >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# ## Add part with cleaning tmp files when there are more input matrix files
# command = "rm -f "+snakemake.params.mtx+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
