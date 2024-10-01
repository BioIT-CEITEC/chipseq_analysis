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

if os.path.isfile(snakemake.input.peaks) and os.path.getsize(snakemake.input.peaks) > 0:
  command = "(time Rscript "+snakemake.params.rscript+\
            " "+snakemake.params.odir+\
            " "+snakemake.params.prefix+\
            " "+snakemake.output.Rsam+\
            " "+snakemake.input.peaks+\
            " "+snakemake.input.reads+\
            ") >> "+snakemake.log.run+" 2>&1"
else:
  command = "touch "+(" ".join(snakemake.output))+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

