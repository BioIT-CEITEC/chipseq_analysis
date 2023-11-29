#########################################
# wrapper for rule: plot_bdgdiff_venn
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

import matplotlib.pyplot as plt
from matplotlib_venn import venn2

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: plot_bdgdiff_venn \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

# -1 because of the header
nl_c1 = max(0,sum(1 for line in open(snakemake.input.c1))-1)
nl_c2 = max(0,sum(1 for line in open(snakemake.input.c2))-1)
nl_bt = max(0,sum(1 for line in open(snakemake.input.bt))-1)

sw = snakemake.wildcards
venn2(subsets = (nl_c1, nl_c2, nl_bt), set_labels = (sw.c1, sw.c2))
plt.title("{c1} vs {c2}.{dups}".format(c1=sw.c1, c2=sw.c2, dups=sw.dups))
plt.savefig(snakemake.output.venn)

f = open(snakemake.log.run, 'at')
f.write("## INFO: saving table "+snakemake.output.tab+"\n")
f.close()
with open(snakemake.output.tab, 'w') as f:
    f.write('comparison\ttotal\tcond1\tcond2\tboth\n')
    f.write('{c1}_vs_{c2}\t{total}\t{up}\t{down}\t{both}'.format(c1=sw.c1, c2=sw.c2, total=nl_c1+nl_c2+nl_bt, up=nl_c1, down=nl_c2, both=nl_bt))
