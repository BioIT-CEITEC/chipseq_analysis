#########################################
# wrapper for rule: call_bdgdiff
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_bdgdiff \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

command = "cat "+snakemake.input.xls1+" | grep -qF '# Paired-End mode is on'"+\
          " && echo $(cat "+snakemake.input.xls1+"|grep -F '# total fragments in treatment:'|sed 's/ //g'|cut -f2 -d':')"+\
          " || echo $(cat "+snakemake.input.xls1+"|grep -F '# total tags in treatment:'|sed 's/ //g'|cut -f2 -d':')"
# command = "cat "+snakemake.input.xls1+" | grep -F '# total "+tags+" in treatment:' | sed 's/ //g' | cut -f 2 -d ':'"
depth = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## RESULT (Cond1 depth): "+depth+"\n")
f.close()
d1 = int(depth)

command = "cat "+snakemake.input.xls2+" | grep -qF '# Paired-End mode is on'"+\
          " && echo $(cat "+snakemake.input.xls2+"|grep -F '# total fragments in treatment:'|sed 's/ //g'|cut -f2 -d':')"+\
          " || echo $(cat "+snakemake.input.xls2+"|grep -F '# total tags in treatment:'|sed 's/ //g'|cut -f2 -d':')"
# command = "cat "+snakemake.input.xls2+" | grep -F '# total "+tags+" in treatment:' | sed 's/ //g' | cut -f 2 -d ':'"
depth = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.write("## RESULT (Cond2 depth): "+depth+"\n")
f.close()
d2 = int(depth)


command = "$(which time) macs2 bdgdiff"+\
          " --t1 "+snakemake.input.trt1+\
          " --t2 "+snakemake.input.trt2+\
          " --c1 "+snakemake.input.ctl1+\
          " --c2 "+snakemake.input.ctl2+\
          " -C "+str(snakemake.params.logLR_cutoff)+\
          " -l "+str(snakemake.params.min_len)+\
          " -g "+str(snakemake.params.max_gap)+\
          " --d1 "+str(d1)+\
          " --d2 "+str(d2)+\
          " -o "+snakemake.params.cnd1+" "+snakemake.params.cnd2+" "+snakemake.params.comn+\
          " >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cat "+snakemake.params.cnd1+" | tail -n +2 | awk -v OFS='\t' '{{print $1,$2,$3,$5,$4}}' > "+snakemake.output.cnd1+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cat "+snakemake.params.cnd2+" | tail -n +2 | awk -v OFS='\t' '{{print $1,$2,$3,$5,$4}}' > "+snakemake.output.cnd2+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "cat "+snakemake.params.comn+" | tail -n +2 | awk -v OFS='\t' '{{print $1,$2,$3,$5,$4}}' > "+snakemake.output.comn+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -f "+snakemake.params.cnd1+" "+snakemake.params.cnd2+" "+snakemake.params.comn+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
