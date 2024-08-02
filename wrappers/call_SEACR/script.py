#########################################
# wrapper for rule: call_SEACR
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
import pandas

log_filename = str(snakemake.log)
shell.executable("/bin/bash")

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: call_SEACR \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

input_line = snakemake.input.trt
if hasattr(snakemake.input, 'ctl'):
  input_line += " "+snakemake.input.ctl
else:
  input_line += " "+str(snakemake.params.ctl_thr)
            
command = "$(which time) SEACR_1.3.sh "+input_line+" "+snakemake.params.norm+\
          " relaxed "+snakemake.params.prefix+\
          " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

#command = "mv "+snakemake.params.prefix+".relaxed.bed "+snakemake.output.tab_all+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## INFO: Adding peak_id column and renaming '"+snakemake.params.prefix+".relaxed.bed' to '"+snakemake.output.tab_all+"'\n")
f.close()
tab = pandas.read_csv(snakemake.params.prefix+".relaxed.bed", sep="\t", header=None, names=["chr","start","end","score","summit_cov","summit_pos"])
tab.insert(3, "peak_id", ["peak_"+str(i) for i in range(1,len(tab)+1)], True)
tab.to_csv(snakemake.output.tab_all, sep="\t", header=False, index=False)

command = "$(which time) SEACR_1.3.sh "+input_line+" "+snakemake.params.norm+\
          " stringent "+snakemake.params.prefix+\
          " >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "mv "+snakemake.params.prefix+".stringent.bed "+snakemake.output.tab+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## INFO: Adding peak_id column and renaming '"+snakemake.params.prefix+".stringent.bed' to '"+snakemake.output.tab+"'\n")
f.close()
tab = pandas.read_csv(snakemake.params.prefix+".stringent.bed", sep="\t", header=None, names=["chr","start","end","score","summit_cov","summit_pos"])
tab.insert(3, "peak_id", ["peak_"+str(i) for i in range(1,len(tab)+1)], True)
tab.to_csv(snakemake.output.tab, sep="\t", header=False, index=False)

command = "$(which time) bedtools genomecov -bg -i "+snakemake.output.tab_all+" -g "+snakemake.input.ref+\
            " > "+snakemake.output.bdg_all+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) bedtools genomecov -bg -i "+snakemake.output.tab+" -g "+snakemake.input.ref+\
            " > "+snakemake.output.bdg+" 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# bedGraphToBigWig ${i} /mnt/ssd/ssd_3/references/saccharomyces_cerevisiae/R64-1-1.100/seq/chrom.sizes ${i%.bdg}.bigWig
command = "$(which time) bedGraphToBigWig "+snakemake.output.bdg_all+" "+snakemake.input.ref+" "+snakemake.output.bwg_all+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) bedGraphToBigWig "+snakemake.output.bdg+" "+snakemake.input.ref+" "+snakemake.output.bwg+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
