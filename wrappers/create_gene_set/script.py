#########################################
# wrapper for rule: create_gene_set
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: create_gene_set \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

gene_set = snakemake.input.gset
if ',' in gene_set:
    f=open(snakemake.params.list,'w')
    for el in gene_set.split(','):
        f.write(el+'\n')
    f.close()
    
    # grep -w -f list_of_genes.Vanacova_ChIP-seq R64-1-1.100.gtf | gffread - --bed -o R64-1-1.100.Vanacova_ChIP-seq.bed
    # command = "grep -w -i -f "+snakemake.params.list+" <(awk '$3 == \"gene\"' "+snakemake.input.ref+") | gffread - --bed -o "+snakemake.output.ref+" >> "+snakemake.log.run+" 2>&1"
    # command = "grep -w -i -f "+snakemake.params.list+" "+snakemake.input.ref+" | gffread - --bed -o "+snakemake.output.ref+" >> "+snakemake.log.run+" 2>&1"
    command = "grep -w -i -f "+snakemake.params.list+" <(awk '$3==\"gene\"' "+snakemake.input.ref+") | awk '{{print $1,$4,$5,gensub(/.*gene_name \"([^\"]+)\";.*/,\"\\\\1\",\"g\",$9),$6,$7}}' FS='\t' OFS='\t' > "+snakemake.output.ref+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
else:
    if gene_set.endswith(".gtf"):
        # command = "gffread "+gene_set+" --bed -o "+snakemake.output.ref+" >> "+snakemake.log.run+" 2>&1"
        command = "awk '$3==\"gene\" {{print $1,$4,$5,gensub(/.*gene_id \"([^\"]+)\";.*/,\"\\\\1\",\"g\",$9),$6,$7}}' FS='\t' OFS='\t' "+gene_set+\
                  " > "+snakemake.output.ref+" 2>> "+snakemake.log.run
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    elif gene_set.endswith(".gff") or gene_set.endswith(".gff3"):
        command = "awk '$3==\"gene\" {{print $1,$4,$5,gensub(/.*ID=([^;]+);.*/,\"\\\\1\",\"g\",$9),$6,$7}}' FS='\t' OFS='\t' "+gene_set+\
                  " > "+snakemake.output.ref+" 2>> "+snakemake.log.run
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    elif gene_set.endswith(".bed") or gene_set.endswith(".bed12"):
        command = "cp "+gene_set+" "+snakemake.output.ref+" >> "+snakemake.log.run+" 2>&1"
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    elif gene_set.endswith(".list"):
        command = "grep -w -i -f "+gene_set+" <(awk '$3==\"gene\"' "+snakemake.input.ref+")"+\
                  " | awk '{{print $1,$4,$5,gensub(/.*gene_name \"([^\"]+)\";.*/,\"\\\\1\",\"g\",$9),$6,$7}}' FS='\t' OFS='\t'"+\
                  " > "+snakemake.output.ref+" 2>> "+snakemake.log.run
        f = open(snakemake.log.run, 'at')
        f.write("## COMMAND: "+command+"\n")
        f.close()
        shell(command)
    else:
        raise ValueError(gene_set+" has unsupported extension! Use one of [gtf, gff3, gff, bed, bed12, list].")
