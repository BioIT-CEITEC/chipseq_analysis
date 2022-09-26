#########################################
# wrapper for rule: call_MSPC
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_MSPC \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

# Tec 0.001 0.01 0.05 10 outputDir inputs
command = "(time Rscript "+snakemake.params.rscript+\
          " "+snakemake.params.rep_type+\
          " "+str(snakemake.params.strong)+\
          " "+str(snakemake.params.weak)+\
          " "+str(snakemake.params.fdr)+\
          " "+str(snakemake.threads)+\
          " "+snakemake.output.bed+\
          " "+snakemake.input+\
          " ) >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.trt_bdg+" "+snakemake.output.trt_bdg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.ctl_bdg+" "+snakemake.output.ctl_bdg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# rename original output file
command = "mv "+snakemake.params.xls_tab+" "+snakemake.params.xls_tab_all+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# rename original output file with no cutof
command = "mv "+snakemake.params.sum_tab+" "+snakemake.params.sum_tab_all+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# rename original output file with no cutof
command = "mv "+snakemake.params.nar_tab+" "+snakemake.params.nar_tab_all+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# use Qvalue cutof
command = "(time awk -v OFS='\t' -F '\t' '$9 > -log("+str(snakemake.params.qval_cutof)+")/log(10)' "+snakemake.params.nar_tab_all+" > "+snakemake.output.nar_tab+") 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# use Qvalue cutof
command = "(time awk -v OFS='\t' -F '\t' '$5 > -log("+str(snakemake.params.qval_cutof)+")/log(10)' "+snakemake.params.sum_tab_all+" > "+snakemake.output.sum_tab+") 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# bedGraphToBigWig ${i} /mnt/ssd/ssd_3/references/saccharomyces_cerevisiae/R64-1-1.100/seq/chrom.sizes ${i%.bdg}.bigWig
command = "(time bedGraphToBigWig "+snakemake.output.trt_bdg+" "+snakemake.input.ref+" "+snakemake.output.trt_bwg+") >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "(time bedGraphToBigWig "+snakemake.output.ctl_bdg+" "+snakemake.input.ref+" "+snakemake.output.ctl_bwg+") >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
