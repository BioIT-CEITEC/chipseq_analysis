#########################################
# wrapper for rule: plot_selected_peak_profiles
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
from PyPDF2 import PdfFileMerger, PdfFileReader
import pdfkit

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: plot_selected_peak_profiles \n##\n")
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
    pdfkit.from_string('Input file '+snakemake.input.bed+' is missing or empty!', snakemake.output.plot)
    
else:
    os.makedirs(snakemake.params.prefix, exist_ok=True)
    f = open(snakemake.log.run, 'at')
    f.write("## INFO: creating directory "+snakemake.params.prefix+" if doesn't exist already\n\n")
    f.close()
    
    top_limit = snakemake.wildcards.top
    if int(top_limit) < 1:
      top_limit = 20
    elif int(top_limit) > 200:
      top_limit = 200
      
    # plot_title = "Top "+str(top_limit)+" peaks by score(=-10*log10(Q-value))"
    
    command = "cat "+snakemake.input.bed+" | tail -n +2 | sort -k5,5nr | awk 'NR <= "+str(top_limit)+"' > "+snakemake.params.bed+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    extra = "--sortRegions keep"
    if "skip_zeros" in snakemake.params and snakemake.params.skip_zeros:
        extra += " --skipZeros"
    
    bed_list = []
    with open(snakemake.params.bed, 'r') as f:
      for line in f.readlines():
        sline = line.split('\t')
        peak_ID = sline[3].split('_')[-1]
        pos = str(sline[0])+':'+str(sline[1])+'-'+str(sline[2])
        name = pos+'_peak_'+str(peak_ID)
        plot_title = "Peak "+str(peak_ID)+" at position "+pos
        region_label = 'score:'+str(sline[4])+' logFC:'+str(sline[6])
        bed_list.append(snakemake.params.prefix+name)
        with open(bed_list[-1]+'.bed', 'w') as o:
          o.write(line+os.linesep)
    
        command = "computeMatrix scale-regions -S "+snakemake.input.bwg+" -R "+bed_list[-1]+'.bed'+" -o "+snakemake.params.mtx+" "+extra+" --smartLabels -a "+str(snakemake.params.after)+" -b "+str(snakemake.params.before)+" -p "+str(snakemake.threads)+" >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    
        command = "plotProfile -m "+snakemake.params.mtx+" -o "+bed_list[-1]+'.profile.pdf'+" --outFileNameData "+bed_list[-1]+'.profile.data'+" --perGroup --regionsLabel \""+region_label+"\" --plotTitle \""+plot_title+"\" --samplesLabel \""+snakemake.params.sample_name+"\" --dpi "+str(snakemake.params.dpi)+" --averageType "+snakemake.params.profile_avrg+" --plotType "+snakemake.params.profile_type+" --plotHeight 14 --plotWidth 22 --startLabel start --endLabel end >> "+snakemake.log.run+" 2>&1"
        # command = "plotProfile -m "+snakemake.params.mtx+" -o "+bed_list[-1]+'.pdf'+" --perGroup --regionsLabel \""+region_label+"\" --plotTitle \""+plot_title+"\" --samplesLabel \""+snakemake.params.sample_name+"\" --dpi "+str(snakemake.params.dpi)+" --averageType "+snakemake.params.profile_avrg+" --plotType "+snakemake.params.profile_type+" --plotHeight 14 --plotWidth 22 --startLabel peak_start --endLabel peak_end >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    
    # Call the PdfFileMerger
    mergedObject = PdfFileMerger()
    
    # Loop over all files you want to merge
    for bed in bed_list:
      mergedObject.append(PdfFileReader(bed+'.profile.pdf', 'rb'))
    
    # Define final PDF file
    mergedObject.write(snakemake.output.plot)
    
    command = "rm -f "+snakemake.params.mtx+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

