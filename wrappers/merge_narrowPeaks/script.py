#############################################################
# wrapper for rule: merge_narrowPeaks
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: merge_narrowPeaks \n##\n")
f.close()

command = "rm -f "+snakemake.params.bed+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

command = "for f in "+" ".join(snakemake.input.peaks)+"; do echo \"processing $f\"; cut -f 1-5,7-9 $f >> "+snakemake.params.bed+"; done >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)

# command = "sort -k1,1 -k2,2n "+snakemake.params.bed+" | bedtools merge -i - -c 4,5,6,7,8 -o collapse | awk -v OFS=\"\\t\" -F\"\\t\" 'BEGIN{{counter=0; print \"GeneID\",\"Chr\",\"Start\",\"End\",\"Strand\",\"Peaks\",\"Scores\",\"FoldChange\",\"-Log10pval\",\"-Log10qval\"}}{{counter++; print \"peak\"counter,$1,$2,$3,\".\",$4,$5,$6,$7,$8}}' > "+snakemake.output.ref+" 2>> "+snakemake.log.run
command = "sort -k1,1 -k2,2n "+snakemake.params.bed+" | bedtools merge -i - -c 4,5,6,7,8 -o collapse | awk -v OFS=\"\\t\" -F\"\\t\" 'BEGIN{{counter=0}}{{counter++; print $1,$2,$3,\".\",\"peak\"counter,$4,$5,$6,$7,$8}}' > "+snakemake.output.ref+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)


# command = " Rscript "+os.path.abspath(os.path.dirname(__file__))+"/combine_featureCount_tabs.R "+\
#             snakemake.output.table+ " " +\
#             " ".join(snakemake.input.featureCounts)
# 
# f = open(snakemake.log.run, 'a+')
# f.write("## COMMAND: "+command+"\n")
# shell(command)
