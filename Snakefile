import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

# RNA parameteres processing
#
if not "strandness" in config:
    config["strandness"] = "unstr"

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

if not 'filter_blacklist' in config:
  config['filter_blacklist'] = ['yes']
if not 'top_peaks_neighborhood' in config:
  config['top_peaks_neighborhood'] = ['200']
if not 'top_peaks' in config:
  config['top_peaks'] = ['20']
if not 'peak_profiles_radius' in config:
  config['peak_profiles_radius'] = ['500']
if not 'bdgdiff_logLR_cutoff' in config:
  config['bdgdiff_logLR_cutoff'] = ['3']
if not 'bdgdiff_min_len' in config:
  config['bdgdiff_min_len'] = ['200']
if not "bdgdiff_max_gap" in config:
  config["bdgdiff_max_gap"] = ['100']
if not "effective_genome_size" in config:
    config["effective_genome_size"] = "unk"
if not "fragment_length" in config:
    config["fragment_length"] = "unk"
if not "summary_correlation_method" in config:
    config["summary_correlation_method"] = "spearman"

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

if not config["is_paired"]:
    read_pair_tags = [""]
    paired = "SE"
else:
    read_pair_tags = ["_R1","_R2"]
    paired = "PE"

wildcard_constraints:
     sample = "|".join(sample_tab.sample_name) + "|all_samples",
     lib_name="[^\.\/]+",
     read_pair_tag = "(_R.)?"

##### Target rules #####

rule all:
    input:  "final_report.html"

##### Modules #####

include: "rules/peak_calling.smk"
