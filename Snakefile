import os
import pandas as pd
import json
import numpy as np
from snakemake.utils import min_version

min_version("5.18.0")

configfile: "config.json"

GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

# RNA parameteres processing
#
# if not "strandness" in config:
#     config["strandness"] = "unstr"

# ChIP-seq parameters processing
#
# config['keep_dups'] = ['yes']
# config['rel_profile'] = ['yes'] # [yes, no, only]
# config['gene_sets'] = ['all_genes;non_overlap_500:/mnt/ssd/ssd_3/references/saccharomyces_cerevisiae/R64-1-1.100/annot/R64-1-1.100.non-overlapping_around_500.bed12;group6:TUB1,PMA1,CDC60,HTD2,RPL25,RPP1B']
# config['profile_avrg'] = ['median']
# config['profile_type'] = ['lines']
# config['gene_sets'] = ['all_genes']
# config['top_peaks'] = ['20'] # DONE: should be a pipeline parameter [0 < X < 100:Integer]
# config['top_peaks_neighborhood'] = ['200']
# config['peak_profiles_radius'] = ['500']
# config['bdgdiff_logLR_cutoff'] = ['3']
# config['bdgdiff_min_len'] = ['200']
# config['bdgdiff_max_gap'] = ['100']
# 
# # Peaks annotation part
# config['feat_type'] = ['gene']
# config['annotate_by'] = ['gene_name,gene_id']
# 
# # Differential peak calling part
# config['conditions_to_compare'] = ['ChiP_BY_WT:ChiP_Dtrf4,ChiP_BY_WT:ChiP_Drrp6,ChiP_BY_WT:ChiP_Drai1']
# config['dpc_metric'] = ['ratio,total_depth'] # combination of ["max", "ratio", "total_depth"] separated by comma

if not 'macs_padj_filter' in config:
  config['macs_padj_filter'] = '0.05'
if not 'filter_blacklist' in config:
  config['filter_blacklist'] = 'yes'
if not 'top_peaks_neighborhood' in config:
  config['top_peaks_neighborhood'] = 200
if not 'top_peaks' in config:
  config['top_peaks'] = 20
if not 'peak_profiles_radius' in config:
  config['peak_profiles_radius'] = 500
if not 'rel_profile' in config:
  config['rel_profile'] = 'no' # [yes, no, only]
if not 'gene_sets' in config:
  config['gene_sets'] = 'all_genes'
if not 'profile_avrg' in config:
  config['profile_avrg'] = 'median'
if not 'profile_type' in config:
  config['profile_type'] = 'lines'
if not 'after_genes' in config:
  config['after_genes'] = 500
if not 'before_genes' in config:
  config['before_genes'] = 500
  
if not 'bdgdiff_logLR_cutoff' in config:
  config['bdgdiff_logLR_cutoff'] = '3'
if not 'bdgdiff_min_len' in config:
  config['bdgdiff_min_len'] = '200'
if not "bdgdiff_max_gap" in config:
  config["bdgdiff_max_gap"] = '100'
if not "effective_genome_size" in config:
  config["effective_genome_size"] = "unk"
if not "fragment_length" in config:
  config["fragment_length"] = "unk"
if not "summary_correlation_method" in config:
  config["summary_correlation_method"] = "spearman"
if not "feat_type" in config:
  config['feat_type'] = 'gene'
if not "annotate_by" in config:
  config['annotate_by'] = 'gene_name,gene_id'
if not "UMI" in config:
  config["UMI"] = "no"
if not "dups" in config:
  config["dups"] = "no_dups"
# if not "use_groups" in config:
#   config["use_groups"] = False
if not "macs_broad_peaks" in config:
  config["macs_broad_peaks"] = False
if not "macs_broad_cutof" in config:
  config["macs_broad_cutof"] = 0.1
if not "idr_cutof" in config:
  config['idr_cutof'] = 0.05

#### Reference processing ####
# setting organism from reference
f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","reference.json"),)
reference_dict = json.load(f)
f.close()
config["organism"] = [organism_name.lower().replace(" ","_") for organism_name in reference_dict.keys() if isinstance(reference_dict[organism_name],dict) and config["reference"] in reference_dict[organism_name].keys()][0]

##### Config processing #####
# Folders
#
reference_directory = os.path.join(GLOBAL_REF_PATH,config["organism"],config["reference"])

# Samples
#
sample_tab = pd.DataFrame.from_dict(config["samples"],orient="index")
sample_tab['num_of_reps'] = sample_tab.groupby("condition")["tag"].transform('nunique')
sample_tab = pd.DataFrame(sample_tab.control.str.split(';').tolist(), index=sample_tab.sample_name).stack().reset_index([0, 'sample_name']).merge(sample_tab.drop('control',axis=1), on='sample_name')
sample_tab.rename(columns = {0:'control'}, inplace = True)
sample_tab['control'] = sample_tab.apply(lambda row: 'no_control' if (not row.is_control) & (not row.control) else row.control, axis=1)
sample_tab['name'] = sample_tab.apply(lambda row: row.condition if row.tag == '' else "_".join([row.condition, row.tag]), axis=1)
sample_tab['peaks_name'] = sample_tab.apply(lambda row: '' if row['is_control'] else '_VS_'.join([row['name'],row['control']]), axis=1)
# config["use_groups"] = False
# sample_tab = pd.DataFrame({'sample_name': ['karpas', 'pfeiffer', 'pfeiffer_ctr'],
#                             'condition': ['chip', 'chip2', 'ctrl'],
#                             'tag': ['','','']})
# if config["use_groups"] == False:
#   sample_tab['group'] = ""
# sample_tab['is_control'] = np.where(sample_tab['condition'].isin(["control","ctrl"]), True, False)
# sample_tab['name'] = sample_tab[["group", "condition"]].apply("_".join, axis=1) if config["use_groups"] else sample_tab['condition']
# # sample_tab = sample_tab[sample_tab.is_control == False]
# # sample_tab = sample_tab[sample_tab.group == "karpas"]
# # print(sample_tab)
# sample_tab['control'] = sample_tab.apply(lambda row: sample_tab.loc[(sample_tab.group == row.group) & (sample_tab.is_control == True), 'name'].unique(), axis=1)
print(sample_tab)

#### Setting up the reference gene set ####
default_reference = GLOBAL_REF_PATH+config['organism']+"/"+config['reference']+"/annot/"+config['reference']+".gtf"
if 'gene_sets' in config:
    gene_sets = {y[0]:y[1] if y[0] != "all_genes" else default_reference for y in [x.split(':') for x in config['gene_sets'].split(';')]}
else:
    gene_sets = {'all_genes': default_reference}
print(gene_sets)

#### Setting up wildcard constraints ####
wildcard_constraints:
     # sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     dups="no_dups|keep_dups"

##### Target rules #####

rule all:
    # input:  "final_report.html", "dummy"
    # input: expand("mapped/pseudo/{cond}_rep1.{dups}.bam", cond=sample_tab.condition.unique(), dups=config['dups'])
    # input: expand("results/MACS_peaks/{name}/{name}.{dups}.peaks.all.narrowPeak", name=[i for i in sample_tab.peaks_name.unique() if i != ''], dups=config['dups'])
    input: "dummy"

##### Modules #####

include: "rules/peak_calling.smk"
