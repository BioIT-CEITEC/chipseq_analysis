#########################################
# wrapper for rule: plot_average_peak_profile
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
#import pdfkit
from fpdf import FPDF

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: plot_average_peak_profile \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) == 0:
    f = open(snakemake.log.run, 'at')
    f.write("## WARNING: Input file "+snakemake.input.bed+" is missing or empty\n")
    f.close()
    #pdfkit.from_string('Input file '+snakemake.input.bed+' is missing or empty!', snakemake.output.plot)
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("helvetica", size=16)
    pdf.write(txt='Input file '+snakemake.input.bed+' is missing or empty!')
    pdf.output(snakemake.output.plot)
    
else:
    extra = "--sortRegions keep"
    if "skip_zeros" in snakemake.params and snakemake.params.skip_zeros:
        extra = extra+" --skipZeros"
    title = snakemake.params.title+" ("+snakemake.wildcards.filt+")"
    
    command = "computeMatrix reference-point -S "+snakemake.input.bwg+" -R "+snakemake.input.bed+" -o "+snakemake.params.mtx+" "+extra+" --samplesLabel "+snakemake.params.sample_name+" -a "+str(snakemake.params.after)+" -b "+str(snakemake.params.before)+" -p "+str(snakemake.threads)+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    # ## PLOT ONLY HEATMAP
    # # plotHeatmap -m ${i} -o ${i%.mtx.gz}.heatmap.pdf --plotTitle "$(basename $i)" --dpi 300 --heatmapWidth 22 --startLabel start --endLabel end --whatToShow "heatmap and colorbar"
    # command = "plotHeatmap -m "+snakemake.params.data+" -o "+snakemake.output.heatmap+" --regionsLabel peaks --plotTitle "+title+" --dpi "+str(snakemake.params.dpi)+" --heatmapWidth 22 --refPointLabel peak_summit --whatToShow \"heatmap and colorbar\" >> "+snakemake.log.run+" 2>&1"
    # f = open(snakemake.log.run, 'at')
    # f.write("## COMMAND: "+command+"\n")
    # f.close()
    # shell(command)

    ## PLOT ONLY PROFILE
    # plotProfile -m ${i} -o ${i%.mtx.gz}.profile.pdf --plotTitle "$(basename $i)" --dpi 300 --averageType mean --plotType std --plotHeight 14 --plotWidth 22 --startLabel start --endLabel end
    command = "plotProfile -m "+snakemake.params.mtx+\
              " -o "+snakemake.output.profile+\
              " --outFileNameData "+snakemake.output.profile_data+\
              " --regionsLabel peaks"+\
              " --perGroup"+\
              " --plotTitle \""+title+"\""+\
              " --dpi "+str(snakemake.params.dpi)+\
              " --averageType "+snakemake.params.profile_avrg+\
              " --plotType "+snakemake.params.profile_type+\
              " --plotHeight 14 --plotWidth 22"+\
              " --refPointLabel peak_summit"+\
              " >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    ## PLOT HEATMAP AND PROFILE
    # plotHeatmap -m ${i} -o ${i%.mtx.gz}.combined.pdf --plotTitle "$(basename $i)" --dpi 300 --averageType mean --plotType std --heatmapWidth 12 --startLabel start --endLabel end
    command = "plotHeatmap -m "+snakemake.params.mtx+\
              " -o "+snakemake.output.plot+\
              " --outFileNameMatrix "+snakemake.params.data+\
              " --regionsLabel peaks"+\
              " --plotTitle \""+title+"\""+\
              " --dpi "+str(snakemake.params.dpi)+\
              " --averageType "+snakemake.params.profile_avrg+\
              " --plotType "+snakemake.params.profile_type+\
              " --heatmapWidth 12"+\
              " --refPointLabel peak_summit"+\
              " >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    command = "rm -f "+snakemake.params.mtx+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
