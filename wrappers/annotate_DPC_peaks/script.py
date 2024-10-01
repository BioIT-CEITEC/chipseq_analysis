#########################################
# wrapper for rule: annotate_peaks
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: annotate_peaks \n##\n")
f.close()

shell.executable("/bin/bash")

conda = str(subprocess.Popen("conda list 2>&1", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+conda+"\n")
f.close()

if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) > 0:
    command = "(time Rscript "+snakemake.params.rscript+\
              " "+snakemake.input.bed+\
              " "+snakemake.input.gtf+\
              " "+snakemake.output.bed+\
              " "+snakemake.params.feat_type+\
              " "+snakemake.params.annotate_by+\
              ") >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

    command = "sed -i '1 s/^\([^\t]\+\t\)\{{5\}}/chr\tstart\tend\tlogLR\tname\t/' "+snakemake.output.bed+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

else:
    with open(snakemake.log.run, 'at') as f:
        f.write("## NOTE: "+snakemake.input.bed+" is empty\n")

    command = "touch "+snakemake.output.bed+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

