#########################################
# wrapper for rule: call_macs2
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_macs2 \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

keep_dups = "all"

if snakemake.params.frag_len == "unk":
  nomodel = ""
else:
  nomodel = "--nomodel --extsize "+str(snakemake.params.frag_len)

input_line = "-t "+" ".join(snakemake.input.trt)
if hasattr(snakemake.input, 'ctl'):
  input_line += " -c "+" ".join(snakemake.input.ctl)
  
# if int(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0])
if hasattr(snakemake.params, 'paired'):
  input_line += " -f BAMPE" if snakemake.params.paired else " -f BAM"
            
command = "(time macs2 callpeak "+input_line+\
          " --keep-dup "+keep_dups+\
          " -g "+str(snakemake.params.effective_GS)+\
          " --outdir "+snakemake.params.dir+\
          " --name "+snakemake.params.name+\
          " "+nomodel+\
          " --bdg"+\
          " --tempdir "+snakemake.params.temp+\
          " -q 1 ) >> "+snakemake.log.run+" 2>&1"
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
command = "mv "+snakemake.params.xls_tab+" "+snakemake.output.xls_tab_all+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# rename original output file with no cutof
command = "mv "+snakemake.params.sum_tab+" "+snakemake.output.sum_tab_all+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# rename original output file with no cutof
command = "mv "+snakemake.params.nar_tab+" "+snakemake.output.nar_tab_all+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# use Qvalue cutof
command = "(time awk -v OFS='\t' -F '\t' '$9 > -log("+str(snakemake.params.qval_cutof)+")/log(10)' "+snakemake.output.nar_tab_all+" > "+snakemake.output.nar_tab+") 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# use Qvalue cutof
command = "(time awk -v OFS='\t' -F '\t' '$5 > -log("+str(snakemake.params.qval_cutof)+")/log(10)' "+snakemake.output.sum_tab_all+" > "+snakemake.output.sum_tab+") 2>> "+snakemake.log.run
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
