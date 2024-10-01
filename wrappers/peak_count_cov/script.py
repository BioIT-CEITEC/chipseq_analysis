######################################
# wrapper for rule: peak_count_cov
######################################
import os
import sys
import math
import subprocess
from snakemake.shell import shell

shell.executable("/bin/bash")

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: peak_count_cov \n##\n")
f.close()

version = str(subprocess.Popen("samtools --version 2>&1 | grep samtools",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

version = str(subprocess.Popen("Rscript --version 2>&1 ",shell=True,stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## VERSION: "+version+"\n")
f.close()

command = "mkdir -p "+os.path.dirname(snakemake.output.counts)+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

input_bam = snakemake.input.bam
command = "samtools depth -b "+snakemake.input.ref+" "+input_bam+" > "+snakemake.params.bdg+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

metric = snakemake.wildcards.metric
if metric == "total_depth":
    command = "samtools view -c "+input_bam
    depth = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    f = open(snakemake.log.run, 'at')
    f.write("## Library depth: "+depth+" (counted as: "+command+")\n")
    f.close()
else:
    depth = ""

command = "Rscript "+snakemake.params.rscript+" "+snakemake.input.ref+" "+snakemake.params.bdg+" "+snakemake.output.counts+" "+metric+" "+depth+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "rm -f "+snakemake.params.bdg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)


'''
html_report("""
---
title: featureCounts report
---

Report from **featureCounts** tool about reads counting. For more information, please, follow
[this link](http://bioinf.wehi.edu.au/featureCounts/ "reference from: 24. August 2017").

### Input files:

"""+
unfold_list('* [BAM file (aligned to genome)](../mapped/',snakemake.input.bam,')')+
unfold_list('* [genomic reference info file](../',snakemake.input.ref_info,')')+
"""\n
### Process run:

"""+
unfold_list('* [log file](',snakemake.log.run,')')+
"""\n
### Output files:

"""+
unfold_list('* [counts - full (text file)](',snakemake.output.txt,')')+
unfold_list('* [counts - summary (text file)](',snakemake.output.txt,'.summary)')+
"""

----
[visual report by **MultiQC**](../sample_final_reports/"""+snakemake.wildcards.sample+""".RNA.multiqc_report.html#featurecounts)

Return to [start page](../"""+snakemake.wildcards.run_name+""".final_report.html)
""",snakemake.output.html)


# OLD STUFF:
#     run:
#         shell(" {FEATURECOUNTS} -t {params.feature} -g gene_id {params.paired} -s {params.strand} -T {threads} -a {input.ref} -o {output.txt} {input.bam} > {log.run} 2>&1 ")
#         html_report("""
# ---
# title: featureCounts report
# ---
#
# Report from **featureCounts** tool about reads counting. For more information, please, follow
# [this link](http://bioinf.wehi.edu.au/featureCounts/ "reference from: 24. August 2017").
#
# ### Input files:
#
# """+
# unfold_list('* [BAM file (aligned to genome)](../mapped/',input.bam,')')+
# unfold_list('* [genomic reference info file](../',input.ref_info,')')+
# """\n
# ### Process run:
#
# """+
# unfold_list('* [log file](',log.run,')')+
# """\n
# ### Output files:
#
# """+
# unfold_list('* [counts - full (text file)](',output.txt,')')+
# unfold_list('* [counts - summary (text file)](',output.txt,'.summary)')+
# """
#
# ----
# [visual report by **MultiQC**](../sample_final_reports/"""+wildcards.sample+""".RNA.multiqc_report.html#featurecounts)
#
# Return to [start page](../"""+wildcards.run_name+""".final_report.html)
# """,output.html)
'''
