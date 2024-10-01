#############################################################
# wrapper for rule: merge_counts_table/
#############################################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: merge_counts_table/ \n##\n")
f.close()

command = " Rscript "+snakemake.params.rscript+" "+\
            snakemake.output.table+ " " +\
            " ".join(snakemake.input.counts)

f = open(snakemake.log.run, 'a+')
f.write("## COMMAND: "+command+"\n")
shell(command)
