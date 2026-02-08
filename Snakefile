import os
import pandas as pd
import json
import numpy as np
import shutil
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = config["globalResources"]
GLOBAL_TMPD_PATH = config["globalTmpdPath"]
#GLOBAL_REF_PATH = "/mnt/references/"
#GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

##### BioRoot utilities #####
module BR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="bioroots_utilities.smk",branch="master")
    config: config

use rule * from BR as BR_*

##### Config processing #####

config = BR.load_organism()

sample_tab = BR.load_sample()
sample_tab['num_of_reps'] = sample_tab.groupby("condition")["tag"].transform('nunique')
sample_tab = pd.DataFrame(sample_tab.control.str.split(';').tolist(), index=sample_tab.sample_name).stack().reset_index([0, 'sample_name']).merge(sample_tab.drop('control',axis=1), on='sample_name')
sample_tab.rename(columns = {0:'control'}, inplace = True)
sample_tab['control'] = sample_tab.apply(lambda row: 'no_control' if (not row.is_control) & (not row.control) else row.control, axis=1)
sample_tab['name'] = sample_tab.apply(lambda row: row.condition if row.tag == '' else "_".join([row.condition, row.tag]), axis=1)
sample_tab['peaks_name'] = sample_tab.apply(lambda row: '' if row['is_control'] else '_VS_'.join([row['name'],row['control']]), axis=1)
print(sample_tab)

# ChIP-seq parameters processing
#
if not 'rel_profile' in config:
  config['rel_profile'] = 'no' # [yes, no, only]
if not "keep_duplicates" in config:
  config["keep_duplicates"] = False
config["dups"] = "keep_dups" if config["keep_duplicates"] else "no_dups"
  
if config['spikein']:
  config['seacr_normalisation'] = "non"

if not 'ignore_regions' in config:
  config['ignore_regions'] = "blacklist,MT,X,Y,not_chr"

## Process various possibilities of ignore_regions config param
reserved_words = ['blacklist','blcklist','blacklst','not_chr','non_chr']
filter_regions_bed = "mapped/filter_regions.bed"
## add all blacklisted regions from ENCODE's list into ignore_regions bed if asked for
if 'blacklist' in config['ignore_regions'] or 'blcklist' in config['ignore_regions'] or 'blacklst' in config['ignore_regions']:
  config['bam_remove_blacklisted'] = True
  blck_bed = config['reference_dir']+"/others/ChIP-seq/blacklist.v2.bed"
  if os.path.isfile(blck_bed):
    print("## INFO: Adding all blacklisted regions from ENCODE's blacklist ("+blck_bed+") into ignore_regions BED file.")
    shutil.copy(blck_bed, filter_regions_bed)
  else:
    print("## INFO: The ENCODE's blacklist BED was not found. Creating an empty ignore_regions BED file.")
    open(filter_regions_bed, 'x').close()
else:
  print("## INFO: Creating an empty ignore_regions BED file.")
  config['bam_remove_blacklisted'] = False
  open(filter_regions_bed, 'x').close()
## add all non-main chromosomes into ignore_regions bed if asked for
if 'not_chr' in config['ignore_regions'] or 'non_chr' in config['ignore_regions']:
  if os.path.isfile(config['organism_chr_sizes']):
    print("## INFO: Adding all alternative (non-main) chromosomes and contigs into ignore_regions BED file.")
    tab = pd.read_table(config['organism_chr_sizes'], header=None, names=['chr','end'])
    ## filter out all main chromosomes starting with number or chr[0-9]+ or X, Y, MT
    ftab = tab[~tab['chr'].str.contains("^([0-9]+|chr([0-9]+|X|Y|M)|X|Y|MT)$")]
    ## append 2nd column with zeros for start and append the table to the ignore_regions bed
    ftab['start'] = 0
    ftab[['chr','start','end']].to_csv(filter_regions_bed, header=False, index=False, sep='\t', mode='a')
## add into ignore_regions bed everything else what user specified with start=0 and end=99999999999
with open(filter_regions_bed, 'a') as b:
  print("## INFO: Adding all user-specified genomic regions into ignore_regions BED file.")
  for reg in config['ignore_regions'].split(','):
    if not reg in reserved_words:
      b.write(reg+'\t0\t99999999999\n')

## Process tlen_range and define proper min_tlen and max_tlen in the config
if not 'tlen_range' in config:
  config['tlen_range'] = "0,1000"
if not ',' in config['tlen_range']:
  val = int(float(config['tlen_range'].replace(" ", "")))
  config['min_tlen'] = 0
  config['max_tlen'] = abs(val)
else:
  val = [abs(int(float(v))) for v in config['tlen_range'].replace(" ", "").split(',') ]
  config['min_tlen'] = min(val)
  config['max_tlen'] = max(val)
print("## INFO: Valid range of reads template length is: ["+str(config['min_tlen'])+":"+str(config['max_tlen'])+"]")

#### Setting up the reference gene set ####
default_reference = config["organism_gtf"]
if 'gene_sets' in config:
    gene_sets = {y[0]:y[1] if y[0] != "all_genes" else default_reference for y in [x.split(':') for x in config['gene_sets'].split(';')]}
else:
    gene_sets = {'all_genes': default_reference}
print(gene_sets)

#### Setting up wildcard constraints ####
wildcard_constraints:
    # sample = "|".join(sample_tab.sample_name) + "|all_samples",
    dups="no_dups|keep_dups"

##### Target rules #####

rule all:
    input:  "final_report.html"

##### Modules #####

include: "rules/peak_calling.smk"

##### BioRoot utilities - prepare reference #####
module PR:
    snakefile: github("BioIT-CEITEC/bioroots_utilities", path="prepare_reference.smk",branch="master")
    config: config

use rule * from PR as PR_*
