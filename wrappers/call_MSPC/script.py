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
f.write("## CONDA:\n"+version+"\n")
f.close()

# Tec 0.001 0.01 0.05 10 outputDir inputs
command = "$(which time) Rscript "+snakemake.params.rscript+\
          " "+snakemake.params.rep_type+\
          " "+str(snakemake.params.strong)+\
          " "+str(snakemake.params.weak)+\
          " "+str(snakemake.params.fdr)+\
          " "+str(snakemake.threads)+\
          " "+snakemake.output.bed+\
          " "+" ".join(snakemake.input)+\
          " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "mv "+snakemake.params.trt_bdg+" "+snakemake.output.trt_bdg+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

