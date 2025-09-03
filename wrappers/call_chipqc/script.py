#########################################
# wrapper for rule: call_chipqc
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
import pandas

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_chipqc \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

# The ChipQC needs more than 1000 reads to count mean read length and it will throw an error otherwise!
reads = str(subprocess.Popen("samtools view -c "+snakemake.input.reads+" 2>> "+snakemake.log.run, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: samtools view -c "+snakemake.input.reads+" 2>> "+snakemake.log.run+"\n")
f.write("## INFO: There are "+str(reads)+" of reads in "+snakemake.input.reads+"\n")
f.close()
if float(reads) < 1000:
  command = "touch "+(" ".join(snakemake.output))+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## INFO: Not enough reads to run ChipQC...touching an empty output files!\n")
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  exit()

if os.path.isfile(snakemake.input.peaks) and sum(1 for line in open(snakemake.input.peaks, 'r') if not (line.startswith('track') or line.startswith('#'))) > 0:
  command = "grep -vP '^(#|track)' "+snakemake.input.peaks+"|awk '{{$2=$2+$10;$3=$2+1;print}}' OFS='\\t' > "+snakemake.params.input_peaks+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

  command = "$(which time) --verbose Rscript "+snakemake.params.rscript+\
            " "+snakemake.params.odir+\
            " "+snakemake.params.prefix+\
            " "+snakemake.output.Rsam+\
            " "+snakemake.params.input_peaks+\
            " "+snakemake.input.reads+\
            " >> "+snakemake.log.run+" 2>&1"
else:
  command = "touch "+(" ".join(snakemake.output))+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

