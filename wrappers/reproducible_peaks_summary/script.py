#########################################
# wrapper for rule: reproducible_peaks_summary
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

print("\n##\n## RULE: reproducible_peaks_summary \n##\n")

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
print("## CONDA: "+version+"\n")

cols = {
  'condition': [],
  'comparison': [],
  'source': [],
  'reps_type': [],
  'N': []
  }
# for ann in snakemake.params.annot.split(","):
#     cols['count_all_annot_by_'+ann] = []
#     cols['count_0.1_annot_by_'+ann] = []
#     cols['count_0.05_annot_by_'+ann] = []
#     cols['count_0.01_annot_by_'+ann] = []
#     cols['count_0.001_annot_by_'+ann] = []

for bed in snakemake.input.bed:
    print("doing: "+bed)
    name = os.path.split(bed)
    cond = os.path.split(name[0])[1]
    rtype = 'true' if 'true_reps' in name[1] else 'pseudo'
    # Following works only in case of paths like 'something/SOURCE_peaks/something...'
    source = bed.split('/')[1][:-6]
    comp = 'rep1_VS_rep2' if source == 'IDR' else 'all_reps'
    if '_VS_' in name[1]:
      comp = re.findall('\w+\.(rep\d+_VS_rep\d+)\.\w+',name[1])
      # print(comp)
    if os.path.isfile(bed) and os.path.getsize(bed) > 0:
        tab = pd.read_table(bed, sep="\t")
        cols['condition'].append(cond)
        cols['comparison'].append(comp)
        cols['source'].append(source)
        cols['reps_type'].append(rtype)
        cols['N'].append(len(tab))
    else:
        cols['condition'].append(cond)
        cols['comparison'].append(comp)
        cols['source'].append(source)
        cols['reps_type'].append(rtype)
        cols['N'].append(0)
        
print("\n")
# print(cols)

df = pd.DataFrame(data=cols).sort_values(by = 'condition')
# sorted_cols = sorted(list(df.columns))
# sorted_cols.remove('condition')
# sorted_cols = ['condition']+sorted_cols
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

