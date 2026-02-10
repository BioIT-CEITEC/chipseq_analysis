#########################################
# wrapper for rule: filter_bam
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")
log_filename = str(snakemake.log)

f = open(log_filename, 'a+')
f.write("\n##\n## RULE: filter_bam \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(log_filename, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

bad_tags = 2828 # read_unmapped + mate_unmapped + not_primary_alignment + read_fails_platform/vendor_quality_checks + supplementary_alignment
if not snakemake.params.keep_dups:
  bad_tags += 1024 # read_is_PCR_or_optical_duplicate

## Myslim ze by vsetky filtre mali byt prevedene na -e formu a musi sa pridat kontrola na PE a vetveni kvoli kontrole tlen
command = "$(which time) samtools view"+\
          " -@ "+str(snakemake.threads)+\
          " -q "+str(snakemake.params.min_mapq)+\
          " -e '(!flag.paired && tlen == 0) || (flag.paired && tlen != 0 && ((tlen <= "+str(snakemake.params.max_tlen)+" && tlen >= "+str(snakemake.params.min_tlen)+\
           ") || (tlen >= -"+str(snakemake.params.max_tlen)+" && tlen <= -"+str(snakemake.params.min_tlen)+")))'"+\
          " -F "+str(bad_tags)+\
          " -U "+snakemake.params.bam_fail+\
          " -b -h "+snakemake.input.bam+\
          " 2>> "+log_filename+\
          " | "+\
          "$(which time) samtools view"+\
          " -@ "+str(snakemake.threads)+\
          " -L "+snakemake.input.bed+\
          " -U "+snakemake.output.bam+\
          " -b -h -"+\
          " >> "+snakemake.params.bam_fail+\
          " 2>> "+log_filename
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) samtools index -@ "+str(snakemake.threads)+" "+snakemake.output.bam+" >> "+log_filename+" 2>&1"
f = open(log_filename, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
