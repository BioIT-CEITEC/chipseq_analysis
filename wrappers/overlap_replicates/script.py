#########################################
# wrapper for rule: overlap_replicates
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

SAMTOOLS = "samtools"

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: overlap_replicates \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "mkdir -p "+snakemake.params.odir+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "(time Rscript "+snakemake.params.rscript+\
          " "+snakemake.output.tall+\
          " "+snakemake.output.tol+\
          " "+snakemake.output.s1+\
          " "+snakemake.output.s2+\
          " "+snakemake.output.smr+\
          " "+snakemake.wildcards.name+\
          " "+snakemake.wildcards.tool+\
          " "+" ".join(snakemake.input.reps)+\
          ") >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

