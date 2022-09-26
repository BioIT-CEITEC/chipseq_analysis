#########################################
# wrapper for rule: multiqc_report
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: multiqc_report \n##\n")
f.close()

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# command = "mkdir -p "+os.path.dirname(snakemake.output.report) + " >> "+snakemake.log.run+" 2>&1 "
# f = open(snakemake.log.run, 'at')
# f.write("## COMMAND: "+command+"\n")
# f.close()
# shell(command)

cl_conf = "--cl_config \""
cl_conf+= "sp:{{ "
cl_conf+= "macs2:{{fn:\'*.peaks.all.xls\'}},"
cl_conf+= "phantompeakqualtools/out:{{fn:\'*.phantom_peak_calling_stats.tsv\'}},"
cl_conf+= "deeptools/plotProfile:{{fn:\'*.average_peak_profile.data\'}},"
cl_conf+= " }}\""
cl_conf+= " --cl_config \"extra_fn_clean_exts:['.peaks.all.xls','.phantom_peak_calling_stats.tsv','.average_peak_profile.data']\""
cl_conf+= " --cl_config \"fn_ignore_files:['*.profile.data']\""

multiqc_search_paths = "./"
multiqc_search_paths+= " --ignore \""+snakemake.params.repdir+"\""
command = "multiqc -f -n "+snakemake.output.report+" "+cl_conf+" "+multiqc_search_paths+" >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "zip -r "+snakemake.output.zip+" "+snakemake.params.repdir+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

