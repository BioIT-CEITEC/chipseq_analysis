#########################################
# wrapper for rule: plot_FRiP
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell
import pdfkit

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: plot_FRiP\n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) == 0:
    f = open(snakemake.log.run, 'at')
    f.write("## WARNING: There are no peaks in "+snakemake.input.bed+"!\n")
    f.close()
    pdfkit.from_string('Input file '+snakemake.input.bed+' is missing or empty!', snakemake.output.plot)
    
    command = "touch "+snakemake.output.counts+" >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
    
else:
    extra = ""
    command = "$(which time) plotEnrichment -b "+" ".join(snakemake.input.bam)+\
              " --BED "+snakemake.input.bed+" --outRawCounts "+snakemake.output.counts+\
              " -o "+snakemake.output.plot+" --smartLabels --plotTitle 'Fraction of reads in filtered peaks' "+extra+\
              " >> "+snakemake.log.run+" 2>&1"
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)

# if os.path.isfile(snakemake.input.bed_all) and os.path.getsize(snakemake.input.bed_all) == 0:
#     f = open(snakemake.log.run, 'at')
#     f.write("## WARNING: There are no peaks in "+snakemake.input.bed_all+"!\n")
#     f.close()
#     pdfkit.from_string('Input file '+snakemake.input.bed_all+' is missing or empty!', snakemake.output.plot_all)
#     
# else:
#     extra = ""
#     command = "plotEnrichment -b "+" ".join(snakemake.input.bam)+" --BED "+snakemake.input.bed_all+" --outRawCounts "+snakemake.output.counts_all+" -o "+snakemake.output.plot_all+" --smartLabels --plotTitle 'Fraction of reads in all peaks' "+extra+" >> "+snakemake.log.run+" 2>&1"
#     f = open(snakemake.log.run, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
