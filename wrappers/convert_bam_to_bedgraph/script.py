#########################################
# wrapper for rule: convert_bam_to_bedgraph
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

log_filename = str(snakemake.log)
shell.executable("/bin/bash")

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: convert_bam_to_bedgraph \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

paired = True
command = 'samtools view '+snakemake.input.bam+' 2>> '+log_filename+' | head -1 2>> '+log_filename+' | cut -f 2 2>> '+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
flag = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f.write("## FLAG:"+flag+"\n")
f.write("## INFO: file "+snakemake.input.bam+" is "+("paired-end" if int(flag)%2==1 else "single-end")+"\n")
f.close()
if int(flag)%2==0:
    paired = False
  
scaling_infix = ""  
if snakemake.params.spikein:    
  if paired:
    command = "$(which time) --verbose samtools view -uh -e 'flag.proper_pair' -@ "+str(snakemake.threads)+" "+snakemake.input.sbam+" 2>> "+log_filename+\
              " | $(which time) --verbose samtools sort -n -@ "+str(snakemake.threads)+" 2>> "+log_filename+\
              " | $(which time) --verbose bedtools bamtobed -bedpe -i stdin 2>> "+log_filename+\
              " | awk '$1==$4' 2>> "+log_filename+\
              " | wc -l 2>> "+log_filename+\
              " | cut -f 1 -d ' ' 2>> "+log_filename
  else: 
    command = "$(which time) --verbose samtools view -c -@ "+str(snakemake.threads)+" "+snakemake.input.sbam+" 2>> "+log_filename
  f = open(log_filename, 'at')
  f.write("## COMMAND: "+command+"\n")
  spike_frags = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
  f.write("## INFO: Spike-in fragments: "+str(spike_frags)+"\n")
  scale_factor = round(float(snakemake.params.scalefac)/int(spike_frags), 8)
  f.write("## INFO: Spike-in scale factor: "+str(scale_factor)+"\n")
  f.close()
  
  scaling_infix = " -scale "+str(scale_factor)
  
if paired:
  command = "$(which time) --verbose samtools view -uh -e 'flag.proper_pair' -@ "+str(snakemake.threads)+" "+snakemake.input.bam+" 2>> "+log_filename+\
            " | $(which time) --verbose samtools sort -n -@ "+str(snakemake.threads)+" 2>> "+log_filename+\
            " | $(which time) --verbose bedtools bamtobed -bedpe -i stdin 2>> "+log_filename+\
            " | awk '$1==$4' 2>> "+log_filename+\
            " | cut -f 1,2,6 | $(which time) --verbose sort -k1,1 -k2,2n -k3,3n 2>> "+log_filename+\
            " | $(which time) --verbose bedtools genomecov -bg"+scaling_infix+" -i stdin -g "+snakemake.input.ref+\
            " > "+snakemake.output.bdg+" 2>> "+log_filename
else: 
  command = "$(which time) --verbose bedtools genomecov -bg"+scaling_infix+" -ibam "+snakemake.input.bam+" > "+snakemake.output.bdg+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
