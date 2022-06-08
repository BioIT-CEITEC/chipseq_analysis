#########################################
# wrapper for rule: dedup_bam
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: dedup_bam \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()


command = "(time samtools view -b -h -F 1024 -@ "+str(snakemake.threads)+" "+snakemake.input.bam+" > "+snakemake.output.bam+" ) 2>> "+snakemake.log.run+" && (time samtools index -@ "+str(snakemake.threads)+" "+snakemake.output.bam+" ) >> "+snakemake.log.run+" 2>&1"
#command = "bamCoverage -b "+snakemake.input.bam+" -o "+snakemake.output.bw+" -p "+str(snakemake.threads)+" "+ignore_dups+" "+extra+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
