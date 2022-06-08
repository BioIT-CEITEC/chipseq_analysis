import math
import subprocess
import json
import re
import os.path
import glob
import fnmatch
import pandas as pd
import itertools
from snakemake.utils import R
from snakemake.utils import report
from os.path import split


# rule final_report:
#     input:  report = rules.multiqc_report.output.report,
#     output: report = "final_report.html"
#     params: pname = PROJECT_NAME,
#             config = json.dumps(config),
#             macs_padj = config.macs_padj_filter.min(),
#             quality_cutoff = config.quality_cutoff.min(),
#             peak_profiles_radius = config.peak_profiles_radius.min(),
#             top_peaks = config.top_peaks.min(),
#             top_peaks_neighborhood = config.top_peaks_neighborhood.min(),
#             corr_method = config.corr_method.min(),
#             feat_type=config.feat_type.min(),
#             annotate_by = config.annotate_by.min(),
#             logLR_cutoff = config.bdgdiff_logLR_cutoff.min(),
#     conda:  "../wrappers/final_report/env.yaml"   
#     script: "../wrappers/final_report/ChIP-seq_analysis_report_template.Rmd" 
# 
# def specify_inputs_for_multiqc_report(wc):
#     inputs = [ 
#         "bam_QC/correlation_heatmap.no_dups.pdf",
#         expand("mapped/{full_name}.bam_cov.no_dups.bigWig", full_name=set(config["full_name"].tolist())),
#         expand("bam_QC/{name}.{dups}.cross-correlation.pdf", name=set(config.loc[config.condition != 'control', "full_name"].tolist()), dups="no_dups"),
#         expand("bam_QC/{name}.{dups}.cross-correlation.pdf", name=set(config.loc[config.condition != 'control', "full_name"].tolist()), dups="keep_dups"),
#     ]
#     # These results are meaningful for the pipeline assesment regardless of the keep_dups parameter value
#     inputs.append("bam_QC/correlation_heatmap.keep_dups.pdf")
#     inputs.append(expand("mapped/{full_name}.bam_cov.keep_dups.bigWig", full_name=config["full_name"].tolist()))
#       
#     # get FRiPs as output
#     dups = "no_dups"
#     if 'keep_dups' in config and config.keep_dups.tolist()[0] == "yes":
#         dups = [dups, "keep_dups"]
#     inputs.append(expand("peaks_QC/fraction_of_reads_in_peaks/{name}.{dups}.FRiP.pdf", name=samples_set, dups=dups))
#     inputs.append(expand("peaks_QC/average_peak_profile/{name}.{dups}.average_peak_profile.pdf", name=samples_set, dups=dups))
#     inputs.append(expand("peaks_QC/top_{top}_peak_profiles_by_score/{name}.{dups}.top_peak_profiles.pdf", top=config.top_peaks.min(), name=samples_set, dups=dups))
#     inputs.append(expand("results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.p-value_score.bdg", name=samples_set, dups=dups))
#     inputs.append(expand("results/peaks_by_macs2/{name}/{name}.{dups}.peaks.annotated.bed", name=samples_set, dups=dups))
#     inputs.append(expand("results/peaks_by_macs2/{name}/{name}.{dups}.peaks.annotated_by_HOMER.tsv", name=samples_set, dups=dups))
#     inputs.append(expand("results/peaks_by_macs2/enriched_peaks_summary.{dups}.tsv", dups=dups))
#     
#     sset = []
#     for sample in set(config.loc[config.condition != "control", "condition"].tolist()):
#         tags = config.loc[config.condition == sample, "tag"].tolist()
#         if len(set(tags)) > 1:
#             tags = set(tags+["pooled", "all"])
#         for tag in tags:
#             tag = tag if tag != "" else "norep"
#             sset.append([sample,tag])
#     tags = config.loc[config.condition != "control", "tag"].tolist()
#     if len(set(tags)) > 1:
#         tags = set(tags+["pooled"])
#     for tag in tags:
#         tag = tag if tag != "" else "norep"
#         sset.append(["all_samples",tag])
#     # print(sset)
#     for pa in config.profile_avrg.min().split(";"):
#         for pt in config.profile_type.min().split(";"):
#             for gs in gene_sets:
#                 if 'rel_profile' in config and config['rel_profile'].tolist()[0] != "only":
#                     inputs.append(["peaks_profile/over_"+gs+"/"+x[0]+"/"+x[1]+"/profile_average_"+pa+".profile_type_"+pt+"/"+x[0]+"."+x[1]+".no_dups.heatmap.pdf" for x in sset])
#                 if 'rel_profile' in config and config['rel_profile'].tolist()[0] in ["yes", "only"]:
#                     inputs.append(["peaks_profile/over_"+gs+"/"+x[0]+"/"+x[1]+"/profile_average_"+pa+".profile_type_"+pt+"/"+x[0]+"."+x[1]+".no_dups.rel_counts.heatmap.pdf" for x in sset])
#                 if 'keep_dups' in config and config.keep_dups.tolist()[0] == "yes":
#                     if 'rel_profile' in config and config['rel_profile'].tolist()[0] != "only":
#                         inputs.append(["peaks_profile/over_"+gs+"/"+x[0]+"/"+x[1]+"/profile_average_"+pa+".profile_type_"+pt+"/"+x[0]+"."+x[1]+".keep_dups.heatmap.pdf" for x in sset])
#                     if 'rel_profile' in config and config['rel_profile'].tolist()[0] in ["yes", "only"]:
#                         inputs.append(["peaks_profile/over_"+gs+"/"+x[0]+"/"+x[1]+"/profile_average_"+pa+".profile_type_"+pt+"/"+x[0]+"."+x[1]+".keep_dups.rel_counts.heatmap.pdf" for x in sset])
#       
#     if 'conditions_to_compare' in config and config.conditions_to_compare.min() != "":
#         if config.conditions_to_compare.min() == "all":
#             # create all pairs of conditions
#             conditions = itertools.combinations(set(config.condition.tolist()),2)
#         else:
#             # use only specified conditions
#             conditions = [x.split(":") for x in config.conditions_to_compare.min().split(",")]
#             
#         for cond in conditions:
#             if len(cond) > 2:
#                 raise NotImplementedError("Only comparison of two conditions is supported!")
#             elif len(cond) < 2:
#                 raise ValueError("Not enough conditions to compare, exactly 2 conditions needed!")
#             else:
#                 c0_tags = config.loc[config.condition == cond[0], "tag"].tolist()
#                 if len(c0_tags)>1:
#                     c0_tags.append("pooled")
#                 c1_tags = config.loc[config.condition == cond[1], "tag"].tolist()
#                 if len(c1_tags)>1:
#                     c1_tags.append("pooled")
#                 for t0 in c0_tags:
#                     for t1 in c1_tags:
#                         if t0 == '' or t1 == '' or t0 == t1:
#                             t0 = t0 if t0 != '' else 'norep'
#                             t1 = t1 if t1 != '' else 'norep'
#                             
#                             # # NOTE: DPC_tables are commented out intentionally, because they are computed out of raw reads coverage from BAMs without extension to fragment length
#                             # #       so these results don't correspond to reality as good as data from results/peaks_by_macs2
#                             # for metric in config['dpc_metric'].min().split(","):
#                             #     inputs.append("diff_peak_calling/"+cond[0]+"_vs_"+cond[1]+"/"+t0+"_vs_"+t1+".no_dups/DPC_table."+metric+"_based.annotated.tsv")
#                             #     if 'keep_dups' in config and config.keep_dups.tolist()[0] == "yes":
#                             #         inputs.append("diff_peak_calling/"+cond[0]+"_vs_"+cond[1]+"/"+t0+"_vs_"+t1+".keep_dups/DPC_table."+metric+"_based.annotated.tsv")
#                                     
#                             inputs.append("diff_peak_calling/"+cond[0]+"_vs_"+cond[1]+"/"+t0+"_vs_"+t1+".no_dups/bdgdiff/enriched_peaks_venn.pdf")
#                             inputs.append("diff_peak_calling/"+cond[0]+"_vs_"+cond[1]+"/"+t0+"_vs_"+t1+".no_dups/overlap/overlapped_peaks.tsv")
#                             if 'keep_dups' in config and config.keep_dups.tolist()[0] == "yes":
#                                 inputs.append("diff_peak_calling/"+cond[0]+"_vs_"+cond[1]+"/"+t0+"_vs_"+t1+".keep_dups/bdgdiff/enriched_peaks_venn.pdf")
#                                 inputs.append("diff_peak_calling/"+cond[0]+"_vs_"+cond[1]+"/"+t0+"_vs_"+t1+".keep_dups/overlap/overlapped_peaks.tsv")
#     return inputs
#     
rule multiqc_report:
    input:  expand("results/peaks_by_macs2/enriched_peaks_summary.{dups}.tsv", dups=["no_dups"]),
    output: report = "final_report.html",
    log:    run = "logs/multiqc_report.log",
    conda:  "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"

    
#########################################
# Peak calling
#

rule peaks_summary:
    input:  bed = lambda wcs: expand("results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.annotated_by_HOMER.tsv", name=sample_tab.loc[sample_tab.is_control==False, "name"].unique(), dups=wcs.dups),
    output: tab = "results/peaks_by_macs2/enriched_peaks_summary.{dups}.tsv",
    log:    run = "logs/peaks_summary.{dups}.log",
    params: annot = config["annotate_by"],
    conda:  "../wrappers/peaks_summary/env.yaml"
    script: "../wrappers/peaks_summary/script.py"
    

# rule calculate_score_tracks:
#     input:  trt = "results/peaks_by_macs2/{name}/{name}.{dups}.bdg",
#             ctl = "results/peaks_by_macs2/{name}/{name}.{dups}.control.bdg",
#     output: ppois = "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.p-value_score.bdg",
#             qpois = "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.q-value_score.bdg",
#             maxval= "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.maximum_value_score.bdg",
#             subtr = "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.subtract_from_treatment_score.bdg",
#             FE    = "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.linear_fold_enrichment_score.bdg",
#             logFE = "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.log10_fold_enrichment_score.bdg",
#             logLR = "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.log10_likelihood_ratio_score.bdg",
#             # slogLR= "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}.symmetric_log10_likelihood_ratio_score.bdg", # probably not usefull, as it should be run on two ChiP enriched tracks
#     log:    run = "logs/{name}/calculate_score_tracks.{dups}.log",
#     params: prefix= "results/peaks_by_macs2/{name}/score_tracks/{name}.{dups}",
#             pseudo= 0.00001,
#             stats = "ppois qpois max subtract FE logFE logLR"
#     conda:  "../wrappers/calculate_score_tracks/env.yaml"
#     script: "../wrappers/calculate_score_tracks/script.py"
    

# TODO: be sure to set this rule correctly in SeqUIa to be able to use tag dir of all replicates and/or pooled sample
rule annotate_peaks_by_HOMER:
    input:  bed = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.narrowPeak",
            gtf = expand(reference_directory+"/annot/{ref}.gtf", ref=config['reference'])[0],
            fa  = expand(reference_directory+"/seq/{ref}.fa", ref=config['reference'])[0],
            tagdir = lambda wc: expand("results/HOMER_tag_dirs/{sample}.{dups}/", 
                      sample=sample_tab.loc[sample_tab.group == sample_tab.loc[sample_tab.name == wc.name, 'group'].min(), 'sample_name'],
                      dups=wc.dups),
    output: tsv = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.annotated_by_HOMER.tsv",
            # annstats = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.ann_stats.tsv",
    log:    run = "logs/{name}/annotate_peaks_by_HOMER.{dups}.log",
    resources: mem = 5
    params: rscript = workflow.basedir+"/wrappers/annotate_peaks_by_HOMER/plots_and_stats.R",
            fdr_cutof = config["macs_padj_filter"],
            best = config["top_peaks"],
            annstats = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.ann_stats.tsv",
    conda:  "../wrappers/annotate_peaks_by_HOMER/env.yaml"
    script: "../wrappers/annotate_peaks_by_HOMER/script.py"
    
    
# TODO: be sure to set this rule correctly in SeqUIa to be able to do tag dir for all replicates 
rule create_tag_dir_by_HOMER:
    input:  bam = "mapped/{sample}.{dups}.bam",
            fa  = expand(reference_directory+"/seq/{ref}.fa", ref=config['reference'])[0],
    output: ok = "results/HOMER_tag_dirs/{sample}.{dups}/make_tag_dir_done",
            tagdir = directory("results/HOMER_tag_dirs/{sample}.{dups}/"),
    log:    run = "logs/{sample}.{dups}/create_tag_dir_by_HOMER.log",
    resources: mem = 5,
    conda:  "../wrappers/create_tag_dir_by_HOMER/env.yaml"
    script: "../wrappers/create_tag_dir_by_HOMER/script.py"


# rule annotate_peaks:
#     input:  bed = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.narrowPeak",
#             gtf = expand(reference_directory+"/annot/{ref}.gtf", ref=config['reference'])[0]
#     output: bed = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.annotated.bed",
#     log:    run = "logs/{name}/annotate_peaks.{dups}.log",
#     resources: mem=5
#     params: rscript = workflow.basedir+"/../scripts/annotate_bed_file.R",
#             feat_type=config["feat_type"],
#             annotate_by = config["annotate_by"],
#     conda:  "../wrappers/annotate_peaks/env.yaml"
#     script: "../wrappers/annotate_peaks/script.py"
    

def call_macs2_inputs(wc):
    # dups = "keep_dups"
    # if config["UMI"] != "no" and wc.dups != "keep_dups":
    #     dups = "no_dups"
    dups = wc.dups
    inputs = {
      "ref": expand(reference_directory+"/seq/{ref}.chrom.sizes", ref=config['reference'])[0],
      "trt": expand("mapped/{sample}.{dups}.bam", sample=sample_tab.loc[sample_tab.name == wc.name, "sample_name"].tolist(), dups=dups)
    }
    if sample_tab.loc[sample_tab.name == wc.name, "control"].min():
      control_name = sample_tab.loc[sample_tab.name == wc.name, "control"].min()
      inputs['ctl'] = expand("mapped/{sample}.{dups}.bam", sample=sample_tab.loc[sample_tab.name == control_name[0], 'sample_name'].tolist(), dups=dups)
      
    # print(inputs)
    return inputs
        
rule call_macs2:
    input:  unpack(call_macs2_inputs),
    output: trt_bdg = "results/peaks_by_macs2/{name}/{name}.{dups}.bdg",
            trt_bwg = "results/peaks_by_macs2/{name}/{name}.{dups}.bigWig",
            ctl_bdg = "results/peaks_by_macs2/{name}/{name}.{dups}.control.bdg",
            ctl_bwg = "results/peaks_by_macs2/{name}/{name}.{dups}.control.bigWig",
            xls_tab = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.xls",
            sum_tab = "results/peaks_by_macs2/{name}/{name}.{dups}.summits.bed",
            nar_tab = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.narrowPeak",
            nar_tab_all = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.narrowPeak",
    log:    run = "logs/{name}/call_macs2.{dups}.log",
    threads: 1
    params: trt_bdg = "results/peaks_by_macs2/{name}/{name}.{dups}_treat_pileup.bdg",
            trt_bwg = "results/peaks_by_macs2/{name}/{name}.{dups}_treat_pileup.bigWig",
            ctl_bdg = "results/peaks_by_macs2/{name}/{name}.{dups}_control_lambda.bdg",
            ctl_bwg = "results/peaks_by_macs2/{name}/{name}.{dups}_control_lambda.bigWig",
            xls_tab = "results/peaks_by_macs2/{name}/{name}.{dups}_peaks.xls",
            sum_tab = "results/peaks_by_macs2/{name}/{name}.{dups}_summits.bed",
            nar_tab = "results/peaks_by_macs2/{name}/{name}.{dups}_peaks.narrowPeak",
            xls_tab_all = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.xls",
            sum_tab_all = "results/peaks_by_macs2/{name}/{name}.{dups}.summits.all.bed",
            nar_tab_all = "results/peaks_by_macs2/{name}/{name}.{dups}.peaks.all.narrowPeak",
            effective_GS = config["effective_genome_size"],
            frag_len = config["fragment_length"],
            qval_cutof = config["macs_padj_filter"],
            umi = config["UMI"],
            dir = "results/peaks_by_macs2/{name}/",
            name= "{name}.{dups}",
            temp= "/mnt/ssd/ssd_1/tmp/",
    conda:  "../wrappers/call_macs2/env.yaml"
    script: "../wrappers/call_macs2/script.py"
    
##########################################
# PREPARE INPUT
#

# rule dedup_bam:
#     input:  bam = "mapped/{name}.keep_dups.bam",
#     output: bam = "mapped/{name}.no_dups.bam",
#     log:    run = "logs/{name}/dedup_bam.log"
#     threads: 5,
#     conda:  "../wrappers/dedup_bam/env.yaml"
#     script: "../wrappers/dedup_bam/script.py"
# 
# def filter_and_index_bam_inputs(wcs):
#     if 'filter_blacklist' in config and config["filter_blacklist"] == 'yes':
#         return "mapped/{name}.no_contam.bam"
#     else:
#         return "mapped/{name}.bam"
# 
# rule filter_and_index_bam:
#     input:  bam = filter_and_index_bam_inputs,
#     output: bam = "mapped/{name}.keep_dups.bam",
#             bai = "mapped/{name}.keep_dups.bam.bai",
#     log:    run = "logs/{name}/filter_and_index_bam.log"
#     threads: 5,
#     params: quality_cutof = config["quality_cutoff"],
#     conda:  "../wrappers/filter_and_index_bam/env.yaml"
#     script: "../wrappers/filter_and_index_bam/script.py"
    
