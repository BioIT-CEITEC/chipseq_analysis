#########################################
# wrapper for rule: call_IDR
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_IDR \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# idr --samples <(sort -k9,9nr -T /mnt/ssd/ssd_1/tmp/ macs2/${name}1/${name}1.no_dups.peaks.all.narrowPeak) 
#    <(sort -k9,9nr -T /mnt/ssd/ssd_1/tmp/ macs2/${name}2/${name}2.no_dups.peaks.all.narrowPeak) 
#    --input-file-type narrowPeak --rank q.value --output-file IDR/${name}s/${name}s.no_dups --plot 2>&1 | tee IDR/${name}s/${name}s.IDR_run.no_dups.log
# mv IDR/${name}s/${name}s.no_dups IDR/${name}s/${name}s.no_dups.peaks.all.bed
# awk '$5 >= 540' OFS='\t' IDR/${name}s/${name}s.no_dups.peaks.all.bed > IDR/${name}s/${name}s.no_dups.peaks.bed
command = "$(which time) idr --samples "+\
          " <(sort -k9,9nr -T "+snakemake.params.tmpd+" "+snakemake.input.peaks[0]+")"+\
          " <(sort -k9,9nr -T "+snakemake.params.tmpd+" "+snakemake.input.peaks[1]+")"+\
          " --input-file-type narrowPeak"+\
          "  --rank q.value"+\
          " --output-file "+snakemake.output.all_bed+\
          " --plot "+\
          " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) awk '$5 >= "+str(-125*math.log2(snakemake.params.cutof))+"' OFS='\t' FS='\t' "+snakemake.output.all_bed+\
          " > "+snakemake.output.bed+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

