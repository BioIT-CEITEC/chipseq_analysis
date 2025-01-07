import os
import pandas as pd
import json
import numpy as np
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

use rule * from BR as other_*

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

use rule * from PR as other_*
