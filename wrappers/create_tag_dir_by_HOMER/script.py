#########################################
# wrapper for rule: create_tag_dir_by_HOMER
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
import pandas

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: create_tag_dir_by_HOMER \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

command = "(time makeTagDirectory"+\
          " "+snakemake.output.tagdir+\
          " "+" ".join(snakemake.input.bam)+\
          " -format sam"+\
          " -totalReads all"+\
          " -genome "+snakemake.input.fa+\
          " -checkGC"+\
          ") >> "+snakemake.log.run+" 2>&1 && touch "+snakemake.output.ok+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

