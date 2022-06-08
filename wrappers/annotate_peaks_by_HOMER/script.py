#########################################
# wrapper for rule: annotate_peaks_by_HOMER
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
import pandas

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: annotate_peaks_by_HOMER \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA: "+version+"\n")
f.close()

if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) > 0:
    command = "(time annotatePeaks.pl"+\
              " "+snakemake.input.bed+\
              " "+snakemake.input.fa+\
              " -gtf "+snakemake.input.gtf+\
              " -annStats "+snakemake.params.annstats+\
              " -size given"+\
              " -d "+" ".join(snakemake.input.tagdir)+\
              " > "+snakemake.output.tsv+") 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
    with open(snakemake.log.run, 'at') as f:
        f.write("## NOTE: merging "+snakemake.output.tsv+" with "+snakemake.input.bed+"\n")
    orig = pandas.read_csv(snakemake.input.bed, sep="\t", header=None, names=["chr","start","end","peak_id","score","strand","signal","pvalue","qvalue","summit"])
    new = pandas.read_csv(snakemake.output.tsv, sep="\t", header=0)
    new.rename(columns={new.columns[0]:"peak_id"}, inplace=True)
    new.rename(columns=lambda s:s.replace(" ","_"), inplace=True)
    out = pandas.merge(orig, new, on=["peak_id"])
    out.drop(columns=['Chr','Start','End','Strand','Peak_Score'], inplace=True)
    out.to_csv(snakemake.output.tsv, sep="\t", header=True, index=False)
    
    command = "(time Rscript "+snakemake.params.rscript+\
              " "+snakemake.output.tsv+\
              " "+str(snakemake.params.fdr_cutof)+\
              " "+str(snakemake.params.best)+\
              ") >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
else:
    with open(snakemake.log.run, 'at') as f:
        f.write("## NOTE: "+snakemake.input.bed+" is empty\n")
  
    command = "touch "+snakemake.output.tsv+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

