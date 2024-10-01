#########################################
# wrapper for rule: compute_matrix
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: compute_matrix \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

extra = "--missingDataAsZero"
# if snakemake.params.skip_zeros:
#     extra = "--skipZeros"

if not isinstance(snakemake.input.bwg, list):
    # computeMatrix scale-regions -S $i -R ${REF} -o ${i%.bigWig}.${TAG}.mtx.gz --samplesLabel ${i%%.*} -a $AFTER -b $BEFORE -p 1
    command = "computeMatrix scale-regions -S "+snakemake.input.bwg+" -R "+snakemake.input.ref+" -o "+snakemake.params.mtx+" "+extra+" --samplesLabel "+snakemake.params.sample_name+" -a "+str(snakemake.params.after)+" -b "+str(snakemake.params.before)+" -p "+str(snakemake.threads)+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND:\n"+command+"\n")
    f.close()
    shell(command)
else:
    # computeMatrix scale-regions -S ${s}_rep1.MQ_over_20/${s}_rep1.${i} ${s}_rep2.MQ_over_20/${s}_rep2.${i} -R ${REF} -o "${s}/${s}.${i%.bigWig}.${TAG}.mtx.gz" --samplesLabel "${s}_rep1" "${s}_rep2" -a $AFTER -b $BEFORE -p 1
    names = [os.path.basename(x).replace(".bigWig", "") for x in snakemake.input.bwg]
    command = "computeMatrix scale-regions -S "+" ".join(snakemake.input.bwg)+" -R "+snakemake.input.ref+" -o "+snakemake.params.mtx+" "+extra+" --samplesLabel "+" ".join(names)+" -a "+str(snakemake.params.after)+" -b "+str(snakemake.params.before)+" -p "+str(snakemake.threads)+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

command = "Rscript "+snakemake.params.rscript+" "+snakemake.params.mtx+" "+snakemake.input.ref+" "+snakemake.output.mtx+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# command = "gzip "+snakemake.params.mtx+" >> "+snakemake.log.run+" 2>&1"
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

command = "rm -f "+snakemake.params.mtx+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
    
