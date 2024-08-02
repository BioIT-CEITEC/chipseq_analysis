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

cl_conf = "--cl-config \""
cl_conf+= "sp:{{ "
cl_conf+= "macs2:{{fn:\'*.peaks.all.xls\'}},"
cl_conf+= "phantompeakqualtools/out:{{fn:\'*.phantom_peak_calling_stats.tsv\'}},"
cl_conf+= "deeptools/plotProfile:{{fn:\'*.average_peak_profile.data\'}},"
cl_conf+= " }}\""
cl_conf+= " --cl-config \"extra_fn_clean_exts:['.peaks.all.xls','.phantom_peak_calling_stats.tsv','.average_peak_profile.data']\""
cl_conf+= " --cl-config \"fn_ignore_files:['*.profile.data']\""

multiqc_search_paths = "./"
multiqc_search_paths+= " --ignore '"+snakemake.params.repdir+"' --ignore 'tmp' --ignore '.snakemake' --ignore '.git' --ignore '.taskrunner'"
command = "multiqc -f -z -n "+snakemake.output.report+" "+cl_conf+" "+multiqc_search_paths+" >> "+snakemake.log.run+" 2>&1 "
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
