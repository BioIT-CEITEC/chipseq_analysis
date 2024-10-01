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
f.write("## CONDA: "+version+"\n")
f.close()

command = "(time chipr --input "+" ".join(snakemake.input.peaks)+\
          " --minentries "+str(snakemake.params.minentries)+\
          " --size "+str(snakemake.params.minsize)+\
          " --rankmethod qvalue"+\
          " --output "+snakemake.params.prefix+\
          " ) >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "awk '$9 <= "+str(snakemake.params.cutof)+"' "+snakemake.params.prefix+"_all.bed > "+snakemake.params.prefix+".bed 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.prefix+"_all.bed "+snakemake.params.prefix+".all.bed >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "mv "+snakemake.params.prefix+"_optimal.bed "+snakemake.params.prefix+".optimal.bed >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
