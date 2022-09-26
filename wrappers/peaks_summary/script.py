#########################################
# wrapper for rule: peaks_summary
#########################################
import os
import sys
import subprocess
import re
from snakemake.shell import shell
import pandas as pd
import math
import sys
import statistics

sys.stdout = open(snakemake.log.run, 'a+')

print("\n##\n## RULE: peaks_summary \n##\n")

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA: "+version+"\n")

cols = {
  'name': [],
  'source': [],
  'count_all': [],
  'count_0.1': [],
  'count_0.05': [],
  'count_0.01': [],
  'count_0.001': [],
  'mean_FDR': [],
  'median_FDR': [],
  'min_FDR': [],
  'max_FDR': []
  }
# for ann in snakemake.params.annot.split(","):
#     cols['count_all_annot_by_'+ann] = []
#     cols['count_0.1_annot_by_'+ann] = []
#     cols['count_0.05_annot_by_'+ann] = []
#     cols['count_0.01_annot_by_'+ann] = []
#     cols['count_0.001_annot_by_'+ann] = []

for bed in snakemake.input.bed:
    print("doing: "+bed)
    name = os.path.split(os.path.split(bed)[0])[1]
    # Following works only in case of paths like 'something/SOURCE_peaks/something...'
    source = bed.split('/')[1][:-6]
    if os.path.isfile(bed) and os.path.getsize(bed) > 0:
        tab = pd.read_table(bed, sep="\t")
        cols['name'].append(name)
        cols['source'].append(source)
        cols['count_all'].append(len(tab))
        cols['count_0.1'].append(len(tab.loc[tab.qvalue > -math.log10(0.1)]))
        cols['count_0.05'].append(len(tab.loc[tab.qvalue > -math.log10(0.05)]))
        cols['count_0.01'].append(len(tab.loc[tab.qvalue > -math.log10(0.01)]))
        cols['count_0.001'].append(len(tab.loc[tab.qvalue > -math.log10(0.001)]))
        fdr = [10 ** (-1*x) for x in tab.qvalue]
        cols['mean_FDR'].append(statistics.mean(fdr))
        cols['median_FDR'].append(statistics.median(fdr))
        cols['min_FDR'].append(min(fdr))
        cols['max_FDR'].append(max(fdr))
        # for ann in snakemake.params.annot.split(","):
        #     cols['count_all_annot_by_'+ann].append(len(tab.loc[(tab.qvalue > -math.log10(1)) & (tab[ann])]))
        #     cols['count_0.1_annot_by_'+ann].append(len(tab.loc[(tab.qvalue > -math.log10(0.1)) & (tab[ann])]))
        #     cols['count_0.05_annot_by_'+ann].append(len(tab.loc[(tab.qvalue > -math.log10(0.05)) & (tab[ann])]))
        #     cols['count_0.01_annot_by_'+ann].append(len(tab.loc[(tab.qvalue > -math.log10(0.01)) & (tab[ann])]))
        #     cols['count_0.001_annot_by_'+ann].append(len(tab.loc[(tab.qvalue > -math.log10(0.001)) & (tab[ann])]))
    else:
        cols['name'].append(name)
        cols['source'].append(source)
        cols['count_all'].append(0)
        cols['count_0.1'].append(0)
        cols['count_0.05'].append(0)
        cols['count_0.01'].append(0)
        cols['count_0.001'].append(0)
        cols['mean_FDR'].append(1)
        cols['median_FDR'].append(1)
        cols['min_FDR'].append(1)
        cols['max_FDR'].append(1)
        # for ann in snakemake.params.annot.split(","):
        #     cols['count_all_annot_by_'+ann].append(0)
        #     cols['count_0.1_annot_by_'+ann].append(0)
        #     cols['count_0.05_annot_by_'+ann].append(0)
        #     cols['count_0.01_annot_by_'+ann].append(0)
        #     cols['count_0.001_annot_by_'+ann].append(0)
        
print("\n")
# print(cols)

df = pd.DataFrame(data=cols).sort_values(by = 'name')
# sorted_cols = sorted(list(df.columns))
# sorted_cols.remove('name')
# sorted_cols = ['name']+sorted_cols
# print(sorted_cols)
# df[sorted_cols].to_csv(snakemake.output.tab, sep="\t", header = True, index = False)
df.to_csv(snakemake.output.tab, sep="\t", header = True, index = False)
        

# if os.path.isfile(snakemake.input.bed) and os.path.getsize(snakemake.input.bed) > 0:
#     command = "Rscript "+snakemake.params.rscript+" "+snakemake.input.bed+" "+snakemake.input.gtf+" "+snakemake.output.bed+" "+snakemake.params.feat_type+" "+snakemake.params.annotate_by+" >> "+snakemake.log.run+" 2>&1"
#     f = open(snakemake.log.run, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
#     
#     command = "sed -i '1 s/^\([^\t]\+\t\)\{{10\}}/chr\tstart\tend\tname\tscore\tstrand\tlog2FC\tPvalue\tQvalue\trel_peak_pos\t/' "+snakemake.output.bed+" >> "+snakemake.log.run+" 2>&1"
#     f = open(snakemake.log.run, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)
#     
# else:
#     with open(snakemake.log.run, 'at') as f:
#         f.write("## NOTE: "+snakemake.input.bed+" is empty\n")
#   
#     command = "touch "+snakemake.output.bed+" >> "+snakemake.log.run+" 2>&1"
#     f = open(snakemake.log.run, 'at')
#     f.write("## COMMAND: "+command+"\n")
#     f.close()
#     shell(command)

