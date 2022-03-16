import os
import pandas as pd
import json
from snakemake.utils import min_version

min_version("5.18.0")

GLOBAL_REF_PATH = "/mnt/references/"
GLOBAL_TMPD_PATH = "./tmp/"

os.makedirs(GLOBAL_TMPD_PATH, exist_ok=True)

# DNA parameteres processing
#
if not "lib_ROI" in config:
    config["lib_ROI"] = "wgs"

# RNA parameteres processing
#
if not "strandness" in config:
    config["strandness"] = "unstr"

# if not "count_over" in config:
#     config["count_over"] = "exon"

# if not "max_mapped_reads_to_run_biobloom" in config:
#     config["max_mapped_reads_to_run_biobloom"] = 100

# if not "feature_count" in config:
#     config["feature_count"] = False

# if not "RSEM" in config:
#     config["RSEM"] = False
#     
# ChIP-seq parameters processing
#
if not "effective_genome_size" in config:
    config["effective_genome_size"] = "unk"

if not "fragment_length" in config:
    config["fragment_length"] = "unk"
    
if not "summary_correlation_method" in config:
    config["summary_correlation_method"] = "spearman"

# Reference processing
#
if config["lib_ROI"] != "wgs":
    # setting reference from lib_ROI
    f = open(os.path.join(GLOBAL_REF_PATH,"reference_info","lib_ROI.json"))
    lib_ROI_dict = json.load(f)
    f.close()
    config["reference"] = [ref_name for ref_name in lib_ROI_dict.keys() if isinstance(lib_ROI_dict[ref_name],dict) and config["lib_ROI"] in lib_ROI_dict[ref_name].keys()][0]

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
    input:  "qc_reports/final_alignment_report.html"

##### Modules #####

include: "rules/quality_control.smk"
include: "rules/cross_sample_correlation.smk"
include: "rules/chipseq_specific_qc.smk"
include: "rules/sample_report.smk"
