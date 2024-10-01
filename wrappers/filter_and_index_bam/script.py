#########################################
# wrapper for rule: filter_and_index_bam
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: filter_and_index_bam \n##\n")
f.close()

shell.executable("/bin/bash")

if str(snakemake.params.quality_cutof) == "nan":
  raise Exception("quality_cutoff je na picu")

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "samtools view -@ "+str(snakemake.threads)+" -q "+str(snakemake.params.quality_cutof)+" -b "+snakemake.input.bam+" > "+snakemake.output.bam+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "samtools index -@ "+str(snakemake.threads)+" "+snakemake.output.bam+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
