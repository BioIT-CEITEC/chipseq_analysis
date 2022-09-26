#########################################
# wrapper for rule: compute_matrix_relative
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: compute_matrix_relative \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("R --version 2>&1 | grep 'R version' ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "Rscript "+snakemake.params.rscript+" "+snakemake.input.mtx+" "+snakemake.params.mtx+" "+str(snakemake.params.smooth)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "gzip "+snakemake.params.mtx+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
