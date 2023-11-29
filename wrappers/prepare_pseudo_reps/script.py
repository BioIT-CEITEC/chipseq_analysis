#########################################
# wrapper for rule: prepare_pseudo_reps
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: prepare_pseudo_reps \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

with open(snakemake.log.run, 'at') as f:
  f.write("## INFO: Processing following tags: "+", ".join(snakemake.params.tags)+"\n")

command = "(time samtools merge -f -u -o "+snakemake.params.merged+" "+" ".join(snakemake.input.bam)+\
          " ) >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "(time samtools view -H "+snakemake.params.merged+" > "+snakemake.params.header+") 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

num_reps = str(snakemake.params.num_of_reps)
nlines = str(subprocess.Popen("samtools view -c "+snakemake.params.merged, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## INFO: Total number of lines in "+snakemake.params.merged+" file: "+nlines+"\n")
nlines = str(math.ceil(int(nlines)/int(num_reps)))
f.write("## INFO: Number of lines in each of "+num_reps+" pseudo replicates: "+nlines+"\n")
f.close()

command = "(time samtools view "+snakemake.params.merged+" 2>> "+snakemake.log.run+\
          " | shuf 2>> "+snakemake.log.run+\
          " | split -d -a 1 -l "+nlines+" --numeric-suffixes=1 --additional-suffix='."+snakemake.wildcards.dups+".sam'"+\
          " - "+snakemake.params.prefix+\
          " ) >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

prefix_counter = 1
for rep in ['rep1', 'rep2', 'rep3', 'rep4']:
  if rep in snakemake.params.tags:
      command = "cat "+snakemake.params.header+" "+snakemake.params.prefix+str(prefix_counter)+"."+snakemake.wildcards.dups+".sam"+\
                " | samtools view -b - > "+snakemake.output[rep]+\
                " 2>> "+snakemake.log.run
      prefix_counter += 1
  else:
      command = "touch "+snakemake.output[rep]+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

# if 'rep2' in snakemake.params.tags:
#     command = "cat "+snakemake.params.header+" "+snakemake.params.prefix+"2."+snakemake.wildcards.dups+".sam"+\
#               " | samtools view -b - > "+snakemake.output.r2+\
#               " 2>> "+snakemake.log.run
# else:
#     command = "touch "+snakemake.output.r2+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)
# 
# if 'rep3' in snakemake.params.tags:
#     command = "cat "+snakemake.params.header+" "+snakemake.params.prefix+"3."+snakemake.wildcards.dups+".sam"+\
#               " | samtools view -b - > "+snakemake.output.r3+\
#               " 2>> "+snakemake.log.run
# else:
#     command = "touch "+snakemake.output.r3+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

command = "rm -rf "+snakemake.params.merged+" "+snakemake.params.header+" "+snakemake.params.prefix+"*sam >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
