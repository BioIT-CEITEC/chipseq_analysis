import math
import subprocess
import json
import re
import os.path
import glob
import fnmatch
import pandas as pd
from pandas.api.types import CategoricalDtype
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
#     inputs.append(expand("results/MACS_peaks/{name}/score_tracks/{name}.{dups}.p-value_score.bdg", name=samples_set, dups=dups))
#     inputs.append(expand("results/MACS_peaks/{name}/{name}.{dups}.peaks.annotated.tsv", name=samples_set, dups=dups))
#     inputs.append(expand("results/MACS_peaks/enriched_peaks_summary.{dups}.tsv", dups=dups))
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
#     if 'conds_to_compare' in config and config['conds_to_compare'] != "":
#         if config['conds_to_compare'] == "all":
#             # create all pairs of conditions
#             conditions = itertools.combinations(set(config.condition.tolist()),2)
#         else:
#             # use only specified conditions
#             conditions = [x.split(":") for x in config['conds_to_compare'].split(",")]
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


def multiqc_report_inputs(wc):
    inputs = list()
    inputs+= expand("results/enriched_peaks_summary.{dups}.tsv", dups=config["dups"])
    if any(i > 1 for i in sample_tab['num_of_reps']):
        inputs+= expand("results/reproducible_peaks_summary.{dups}.tsv", dups=config["dups"])
        samples = [*sample_tab.loc[sample_tab.is_control==False, "peaks_name"].unique(), 
                      *sample_tab.loc[(sample_tab.is_control==False)&(sample_tab.num_of_reps>1), "condition"].unique()]
    else:
        print("No replicates!")
        samples = [*sample_tab.loc[sample_tab.is_control==False, "peaks_name"].unique()]
    if 'conds_to_compare' in config and config['conds_to_compare'] != "" and sample_tab.shape[0] > 1:
        inputs+= expand("results/differential_peaks_summary.{dups}.tsv", dups=config['dups'])
    inputs+= expand("results/ChIPQC/{sample}/{sample}.{dups}.report.html",
                sample=sample_tab.loc[sample_tab.is_control==False, "peaks_name"].unique(), #+sample_tab.loc[sample_tab.is_control==False, "condition"].unique(),
                dups=config["dups"])
    inputs+= expand("results/peaks_QC/peak_profiles/over_peaks/{name}.{dups}.from_{tool}.{filt}.average_peak_profile.pdf",
                name=samples,
                dups=config["dups"],
                tool=["MACS"],
                filt=["filtered","all"])
    inputs+= expand("results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.combined.pdf",
                name=samples,
                dups=config["dups"],
                tool=["MACS"],
                filt=["filtered","all"],
                gs=gene_sets.keys())
    inputs+= expand("results/peaks_QC/top_{top}_peak_profiles_by_score/{name}.{dups}.top_peak_profiles.pdf",
                name=samples,
                top=config["top_peaks"],
                dups=config["dups"])
    return inputs

rule multiqc_report:
    input:  multiqc_report_inputs,
    output: report = "final_report.html",
            zip = "final_report_data.zip",
    params: repdir = "final_report_data",
    resources: tmpdir=GLOBAL_TMPD_PATH,
    log:    run = "logs/multiqc_report.log",
    conda:  "../wrappers/multiqc_report/env.yaml"
    script: "../wrappers/multiqc_report/script.py"
    
#########################################
# Diff peak calling
#

def diff_summary_inputs(wc):
    inputs = list()
    if 'conds_to_compare' in config and config['conds_to_compare'] != "":
        if config['conds_to_compare'] == "all":
            # create all pairs of conditions
            conditions = list(itertools.combinations(sample_tab.loc[sample_tab.is_control==False,'condition'].unique(),2))
        else:
            # use only specified conditions
            conditions = [x.split(":") for x in config['conds_to_compare'].split(",")]
        
        # print([i for i in conditions])
        for cond in conditions:
            if len(cond) > 2:
                raise NotImplementedError("Only comparison of two conditions is supported!")
            elif len(cond) < 2:
                raise ValueError("Not enough conditions to compare, exactly 2 conditions needed!")
            elif not cond[0] in sample_tab['condition'].unique():
                raise ValueError(f"No such condition as {cond[0]} in sample conditions!")
            elif not cond[1] in sample_tab['condition'].unique():
                raise ValueError(f"No such condition as {cond[1]} in sample conditions!")
            else:
                inputs.append(f"results/MACS_bdgdiff/{cond[0]}_vs_{cond[1]}/enriched_peaks_venn.{wc.dups}.tsv")
                inputs.append(f"results/overlapped_peaks/{cond[0]}_vs_{cond[1]}/summary_table.{wc.dups}.by_MACS.tsv")
                if all(i > 1 for i in sample_tab.loc[(sample_tab.condition==cond[0])|(sample_tab.condition==cond[1]), 'num_of_reps'].unique()):
                    inputs.append(f"results/overlapped_peaks/{cond[0]}_vs_{cond[1]}/summary_table.{wc.dups}.by_MSPC.tsv")
                    inputs.append(f"results/DiffBind/{cond[0]}_vs_{cond[1]}/summary_table.{wc.dups}.tsv")
    return inputs

rule diff_summary:
    input:  bed = diff_summary_inputs
    output: tab = "results/differential_peaks_summary.{dups}.tsv",
    log:    run = "logs/diff_summary.{dups}.log",
    conda:  "../wrappers/diff_summary/env.yaml"
    script: "../wrappers/diff_summary/script.R"


def call_diffbind_inputs(wc):
  inputs = dict()
  design = f"results/DiffBind/{wc.c1}_vs_{wc.c2}/design_table.{wc.dups}.tsv"
  os.makedirs(os.path.dirname(design), exist_ok=True)
  design_tab = sample_tab.loc[sample_tab['condition'].isin([wc.c1, wc.c2]), ("name","condition","tag","control")].set_axis(['SampleID', 'Condition', 'Replicate', 'ControlID'], axis='columns').drop_duplicates()
  design_tab['bamReads'] = design_tab.apply(lambda row: "mapped/{sample}.{dups}.bam".format(
      sample = sample_tab.loc[sample_tab['name'] == row.SampleID, 'sample_name'].min(), 
      dups = wc.dups
  ), axis=1)
  design_tab['bamControl'] = design_tab.apply(lambda row: "mapped/{sample}.{dups}.bam".format(
      sample = sample_tab.loc[sample_tab['name'] == row.ControlID, 'sample_name'].min(), 
      dups = wc.dups
  ), axis=1)
  design_tab['Peaks'] = design_tab.apply(lambda row: "results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.narrowPeak".format(
      tool = "MACS",
      name = row.SampleID+"_VS_"+row.ControlID,
      dups = wc.dups
  ), axis=1)
  design_tab['PeakCaller'] = "narrow"
  # catgorical order of conditions must be the opposite of normal meaning as DiffBind assumes the first one to be a control condition
  condition_order = CategoricalDtype([wc.c2, wc.c1], ordered=True)
  design_tab['Condition'] = design_tab['Condition'].astype(condition_order)
  design_tab = design_tab.sort_values(['Condition', 'Replicate'])
  design_tab.to_csv(design, sep="\t", index=False)
  
  inputs['bed'] = expand("{peaks}", peaks=design_tab['Peaks'].unique())
  inputs['bam'] = expand("{reads}", reads=[*design_tab['bamReads'].unique(),*design_tab['bamControl'].unique()])
  inputs['tab'] = ancient(design)
  # print(inputs)
  return inputs

rule call_diffbind:
    input:  unpack(call_diffbind_inputs)
    output: deseq_tab = "results/DiffBind/{c1}_vs_{c2}/DESeq2_table.{dups}.full.tsv",
            edger_tab = "results/DiffBind/{c1}_vs_{c2}/edgeR_table.{dups}.full.tsv",
            sum_tab   = "results/DiffBind/{c1}_vs_{c2}/summary_table.{dups}.tsv",
            info_tab  = "results/DiffBind/{c1}_vs_{c2}/info_table.{dups}.tsv",
            cor_heat_1 = "results/DiffBind/{c1}_vs_{c2}/corr_heatmap_score_based.{dups}.png",
            cor_heat_2 = "results/DiffBind/{c1}_vs_{c2}/corr_heatmap_counts_based.{dups}.png",
            cor_heat_deseq = "results/DiffBind/{c1}_vs_{c2}/corr_heatmap_DBA_based_from_DESeq2.{dups}.png",
            cor_heat_edger = "results/DiffBind/{c1}_vs_{c2}/corr_heatmap_DBA_based_from_edgeR.{dups}.png",
            pca_plot_1 = "results/DiffBind/{c1}_vs_{c2}/PCA_counts_based.{dups}.png",
            venn_plot = "results/DiffBind/{c1}_vs_{c2}/Venn_DESeq2_edgeR.{dups}.png",
            pca_plot_deseq = "results/DiffBind/{c1}_vs_{c2}/PCA_DBA_based_from_DESeq2.{dups}.png",
            pca_plot_edger = "results/DiffBind/{c1}_vs_{c2}/PCA_DBA_based_from_edgeR.{dups}.png",
            ma_plot_deseq = "results/DiffBind/{c1}_vs_{c2}/MA_from_DESeq2.{dups}.png",
            ma_plot_edger = "results/DiffBind/{c1}_vs_{c2}/MA_from_edgeR.{dups}.png",
            xy_plot_deseq = "results/DiffBind/{c1}_vs_{c2}/XY_scatter_from_DESeq2.{dups}.png",
            xy_plot_edger = "results/DiffBind/{c1}_vs_{c2}/XY_scatter_from_edgeR.{dups}.png",
            volcano_plot_deseq = "results/DiffBind/{c1}_vs_{c2}/volcano_from_DESeq2.{dups}.png",
            volcano_plot_edger = "results/DiffBind/{c1}_vs_{c2}/volcano_from_edgeR.{dups}.png",
            DBA_heat_deseq = "results/DiffBind/{c1}_vs_{c2}/DBA_heatmap_from_DESeq2.{dups}.png",
            DBA_heat_edger = "results/DiffBind/{c1}_vs_{c2}/DBA_heatmap_from_edgeR.{dups}.png",
    log:    run = "logs/{c1}_vs_{c2}/call_diffbind.{dups}.log",
    params: fdr_cutof = config['diff_fdr_cutof'],
            l2fc_cutof= config['diff_l2fc_cutof'],
            comparison= "{c1}_vs_{c2}",
            top = config['diff_top_dbs'],
    conda:  "../wrappers/call_DiffBind/env.yaml"
    script: "../wrappers/call_DiffBind/script.R"

    
# rule annotate_DPC_peaks:
#     input:  bed = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/DPC_table.{metric}_based.tsv",
#             gtf = expand(REF_DIR+"/{organism}/{ref}/annot/{ref}.gtf", organism=config['organism'].tolist()[0], ref=config['reference'].tolist()[0])[0]
#     output: bed = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/DPC_table.{metric}_based.annotated.tsv",
#     log:    run = ADIR+"/logs/{comparison}/{subset}.{dups}/annotate_DPC_peaks.{metric}_based.log",
#     resources: mem=5
#     params: rscript = workflow.basedir+"/../scripts/annotate_bed_file.R",
#             feat_type=config.feat_type.min(),
#             annotate_by = config.annotate_by.min(),
#     conda:  "../wrappers/annotate_peaks/env.yaml"
#     script: "../wrappers/annotate_peaks/script.py"
# 
#       
# rule DPC_computation:
#     input:  expression_tab = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/counts/merged_counts.{metric}_based.tsv",
#             cfg_tab = ADIR+"/"+PROJECT_NAME+"."+config['analysis_class'].tolist()[0]+".config.json",
#     output: table = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/DPC_table.{metric}_based.tsv",
#     log:    run = ADIR+"/logs/{comparison}/{subset}.{dups}/DPC_computation.{metric}_based.log"
#     params: rscript = workflow.basedir+"/../wrappers/DPC_computation/mrna_de_counts.R",
#     conda:  "../wrappers/DPC_computation/env.yaml"
#     script: "../wrappers/DPC_computation/script.py"
# 
# 
# def merge_counts_table_inputs(wcs):
#     samples = list()
#     if wcs.t1 == "norep":
#         samples.append(config.loc[config.condition == wcs.c1, "full_name"].tolist())
#     elif wcs.t1 == "pooled":
#         for sa in config.loc[config.condition == wcs.c1, "full_name"].tolist():
#             samples.append(sa)
#     else:
#         samples.append(wcs.c1+"_"+wcs.t1)
#     if wcs.t2 == "norep":
#         samples.append(config.loc[config.condition == wcs.c2, "full_name"].tolist())
#     elif wcs.t2 == "pooled":
#         for sa in config.loc[config.condition == wcs.c2, "full_name"].tolist():
#             samples.append(sa)
#     else:
#         samples.append(wcs.c2+"_"+wcs.t2)
#     return(expand(ADIR+"/diff_peak_calling/{c1}_vs_{c2}/{t1}_vs_{t2}.{dups}/counts/{sample}.counts.{metric}_based.tsv", sample=set(samples), c1=wcs.c1, c2=wcs.c2, t1=wcs.t1, t2=wcs.t2, dups=wcs.dups, metric=wcs.metric))
#     
# rule merge_counts_table:
#     input:  counts = merge_counts_table_inputs
#     output: table = ADIR+"/diff_peak_calling/{c1}_vs_{c2}/{t1}_vs_{t2}.{dups}/counts/merged_counts.{metric}_based.tsv",
#     log:    run = ADIR+"/logs/{c1}_vs_{c2}/{t1}_vs_{t2}.{dups}/merge_counts_table.{metric}_based.log",
#     params: rscript = workflow.basedir+"/../wrappers/merge_counts_table/combine_counts_tabs.R",
#     conda:  "../wrappers/merge_counts_table/env.yaml"
#     script: "../wrappers/merge_counts_table/script.py"
#    
#    
# def peak_count_cov_inputs(wcs):
#     if wcs.dups == "no_dups":
#         return(ADIR+"/bams/"+wcs.sample+".no_dups.bam")
#     else:
#         return(ADIR+"/bams/"+wcs.sample+".keep_dups.bam")
#     
# rule peak_count_cov:
#     input:  bam = peak_count_cov_inputs,
#             ref = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/merged_peaks.saf",
#     output: counts = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/counts/{sample}.counts.{metric}_based.tsv",
#     log:    run = ADIR+"/logs/{comparison}/{subset}.{dups}/peak_count_cov_on_{sample}.{metric}_based.log",
#     threads:    1,
#     resources:  mem = 10,
#     params: bdg = ADIR+"/diff_peak_calling/{comparison}/{subset}.{dups}/counts/{sample}.peaks_cov.{metric}_based.bdg",
#             rscript = workflow.basedir+"/../wrappers/peak_count_cov/max_depth.R",
#     conda:  "../wrappers/peak_count_cov/env.yaml"
#     script: "../wrappers/peak_count_cov/script.py"
#   
#     
# def merge_narrowPeaks_inputs(wcs):
#     samples = list()
#     if wcs.t1 == "norep":
#         samples.append(config.loc[config.condition == wcs.c1, "full_name"].tolist())
#     elif wcs.t1 == "pooled":
#         for sa in config.loc[config.condition == wcs.c1, "full_name"].tolist():
#             samples.append(sa)
#     else:
#         samples.append(wcs.c1+"_"+wcs.t1)
#     if wcs.t2 == "norep":
#         samples.append(config.loc[config.condition == wcs.c2, "full_name"].tolist())
#     elif wcs.t2 == "pooled":
#         for sa in config.loc[config.condition == wcs.c2, "full_name"].tolist():
#             samples.append(sa)
#     else:
#         samples.append(wcs.c2+"_"+wcs.t2)
#     return(expand("results/MACS_peaks/{name}/{name}.{dups}.peaks.narrowPeak", name=set(samples), dups=wcs.dups))
#     
# rule merge_narrowPeaks:
#     input:  peaks = merge_narrowPeaks_inputs,
#     output: ref = ADIR+"/diff_peak_calling/{c1}_vs_{c2}/{t1}_vs_{t2}.{dups}/merged_peaks.saf"
#     log:    run = ADIR+"/logs/{c1}_vs_{c2}/{t1}_vs_{t2}.{dups}/merge_narrowPeaks.log"
#     params: bed = ADIR+"/diff_peak_calling/{c1}_vs_{c2}/{t1}_vs_{t2}.{dups}/merged_peaks.bed"
#     conda:  "../wrappers/merge_narrowPeaks/env.yaml"
#     script: "../wrappers/merge_narrowPeaks/script.py"


def overlap_found_peaks_input(wc):
    inputs = dict()
    c1 = wc.c1
    c2 = wc.c2
    if sample_tab.loc[sample_tab.condition==wc.c1, 'num_of_reps'].unique() == 1:
        c1 = sample_tab.loc[sample_tab.condition==wc.c1, 'peaks_name'].unique()[0]
    if sample_tab.loc[sample_tab.condition==wc.c2, 'num_of_reps'].unique() == 1:
        c2 = sample_tab.loc[sample_tab.condition==wc.c2, 'peaks_name'].unique()[0]
    inputs['c1'] = f"results/{wc.tool}_peaks/{c1}/{c1}.{wc.dups}.peaks.all.narrowPeak"
    inputs['c2'] = f"results/{wc.tool}_peaks/{c2}/{c2}.{wc.dups}.peaks.all.narrowPeak"
    return inputs

rule overlap_found_peaks:
    input:  unpack(overlap_found_peaks_input),
    output: tab= "results/overlapped_peaks/{c1}_vs_{c2}/overlapped_peaks.{dups}.by_{tool}.tsv",
            his= "results/overlapped_peaks/{c1}_vs_{c2}/overlapped_peaks.{dups}.by_{tool}.hist.tsv",
            s1 = "results/overlapped_peaks/{c1}_vs_{c2}/singletons_in_{c1}.{dups}.by_{tool}.tsv",
            s2 = "results/overlapped_peaks/{c1}_vs_{c2}/singletons_in_{c2}.{dups}.by_{tool}.tsv",
            smr= "results/overlapped_peaks/{c1}_vs_{c2}/summary_table.{dups}.by_{tool}.tsv",
    log:    run= "logs/{c1}_vs_{c2}/overlap_found_peaks_by_{tool}.{dups}.log",
    params: rscript = workflow.basedir+"/wrappers/overlap_found_peaks/overlap_peaks.R",
            odir = "results/overlapped_peaks/{c1}_vs_{c2}",
            comparison = "{c1}_vs_{c2}",
            fdr_cutof = config['diff_fdr_cutof'],
            l2fc_cutof= config['diff_l2fc_cutof'],
    conda:  "../wrappers/overlap_found_peaks/env.yaml"
    script: "../wrappers/overlap_found_peaks/script.py"
    
    
rule plot_bdgdiff_venn:
    input:  c1 = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_{c1}.{dups}.annotated.bed",
            c2 = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_{c2}.{dups}.annotated.bed",
            bt = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_both.{dups}.annotated.bed",
    output: venn = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_venn.{dups}.pdf",
            tab  = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_venn.{dups}.tsv",
    log:    run  = "logs/{c1}_vs_{c2}/plot_bdgdiff_venn.{dups}.log",
    conda:  "../wrappers/plot_bdgdiff_venn/env.yaml"
    script: "../wrappers/plot_bdgdiff_venn/script.py"


rule annotate_bdgdiff_peaks:
    input:  bed = "results/MACS_bdgdiff/{comparison}/{filename}.{dups}.bed",
            gtf = expand(reference_directory+"/annot/{ref}.gtf", ref=config['reference'])[0],
    output: bed = "results/MACS_bdgdiff/{comparison}/{filename}.{dups}.annotated.bed",
    log:    run = "logs/{comparison}/annotate_bdgdiff_peaks_in_{filename}.{dups}.log",
    resources: mem=5
    params: rscript = workflow.basedir+"/wrappers/annotate_DPC_peaks/annotate_bed_file.R",
            feat_type = config['feat_type'],
            annotate_by = config['peak_annot_by'],
    conda:  "../wrappers/annotate_DPC_peaks/env.yaml"
    script: "../wrappers/annotate_DPC_peaks/script.py"
    
    
def call_bdgdiff_input(wc):
    inputs = dict()
    c1 = wc.c1
    c2 = wc.c2
    if sample_tab.loc[sample_tab.condition==wc.c1, 'num_of_reps'].unique() == 1:
        c1 = sample_tab.loc[sample_tab.condition==wc.c1, 'peaks_name'].unique()[0]
    if sample_tab.loc[sample_tab.condition==wc.c2, 'num_of_reps'].unique() == 1:
        c2 = sample_tab.loc[sample_tab.condition==wc.c2, 'peaks_name'].unique()[0]
    inputs['trt1'] = f"results/MACS_peaks/{c1}/{c1}.{wc.dups}.bdg"
    inputs['ctl1'] = f"results/MACS_peaks/{c1}/{c1}.{wc.dups}.control.bdg"
    inputs['xls1'] = f"results/MACS_peaks/{c1}/{c1}.{wc.dups}.peaks.all.xls"
    inputs['trt2'] = f"results/MACS_peaks/{c2}/{c2}.{wc.dups}.bdg"
    inputs['ctl2'] = f"results/MACS_peaks/{c2}/{c2}.{wc.dups}.control.bdg"
    inputs['xls2'] = f"results/MACS_peaks/{c2}/{c2}.{wc.dups}.peaks.all.xls"
    return inputs
    
rule call_bdgdiff:
    input:  unpack(call_bdgdiff_input),
    output: cnd1 = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_{c1}.{dups}.bed",
            cnd2 = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_{c2}.{dups}.bed",
            comn = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_both.{dups}.bed",
    log:    run  = "logs/{c1}_vs_{c2}/call_bdgdiff.{dups}.log",
    params: logLR_cutoff = config['bdgdiff_logLR_cutoff'], # DONE: should be pipeline parameter >= 0
            min_len = config['bdgdiff_min_len'], # DONE: should be pipeline parameter(Integer)
            max_gap = config['bdgdiff_max_gap'], # DONE: should be pipeline parameter (Integer)
            cnd1 = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_{c1}.{dups}.bed.tmp",
            cnd2 = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_{c2}.{dups}.bed.tmp",
            comn = "results/MACS_bdgdiff/{c1}_vs_{c2}/enriched_peaks_in_both.{dups}.bed.tmp",
    conda:  "../wrappers/call_bdgdiff/env.yaml"
    script: "../wrappers/call_bdgdiff/script.py"

   
#########################################
# Peaks assessment
#

def peaks_summary_inputs(wcs):
    inputs = dict()
    singlerep = sample_tab.loc[sample_tab.is_control==False, "peaks_name"].unique()
    pseudorep = sample_tab.loc[(sample_tab.is_control==False)&(sample_tab.num_of_reps>1), "peaks_name"].unique()+".pseudo_reps"
    multirep  = sample_tab.loc[(sample_tab.is_control==False)&(sample_tab.num_of_reps>1), "condition"].unique()
    if any(i > 1 for i in sample_tab['num_of_reps']):
        inputs['bed'] = expand("results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv", name=singlerep, dups=wcs.dups, tool=['MACS','SEACR']) +\
             expand("results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv", name=pseudorep, dups=wcs.dups, tool=['MACS']) +\
             expand("results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv", name=multirep, dups=wcs.dups, tool=['MACS','MSPC'])
             # expand("results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv", name=multirep, dups=wcs.dups, tool=['MACS','MSPC','Fisher'])
    else:
        inputs['bed'] = expand("results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv", name=singlerep, dups=wcs.dups, tool=['MACS','SEACR'])

    if any(i > 1 for i in sample_tab['num_of_reps']):
        inputs['frip'] = expand("results/peaks_QC/FRiP/{name}.{dups}.from_{tool}.{filt}.FRiP.pdf", 
                            name=[*singlerep,*pseudorep,*multirep], 
                            dups=wcs.dups, 
                            tool=['MACS'],
                            filt=["filtered","all"]) +\
                         expand("results/peaks_QC/FRiP/{name}.{dups}.from_{tool}.{filt}.FRiP.pdf", 
                            name=multirep,  
                            dups=wcs.dups, 
                            tool=['MSPC'],
                            filt=["all"]) +\
                         expand("results/peaks_QC/FRiP/{name}.{dups}.from_{tool}.{filt}.FRiP.pdf", 
                            name=singlerep,  
                            dups=wcs.dups, 
                            tool=['SEACR'],
                            filt=["filtered","all"])
    else:
        inputs['frip'] = expand("results/peaks_QC/FRiP/{name}.{dups}.from_{tool}.{filt}.FRiP.pdf", 
                            name=singlerep, 
                            dups=wcs.dups, 
                            tool=['MACS','SEACR'],
                            filt=["filtered","all"])
    return inputs

rule peaks_summary:
    input:  unpack(peaks_summary_inputs),
    output: tab = "results/enriched_peaks_summary.{dups}.tsv",
    log:    run = "logs/peaks_summary.{dups}.log",
    params: annot = config["peak_annot_by"],
    conda:  "../wrappers/peaks_summary/env.yaml"
    script: "../wrappers/peaks_summary/script.py"
    

# def plot_profile_and_heatmap_input(wc):
#     inputs = dict()
#     if wc.name == "all_samples":
#         if wc.reps == "norep":
#             # we want results for each sample without replicates (i.e., with only one replicate)
#             samples = [x for x in set(cfg.loc[(cfg.tag == "") & (cfg.condition != "control"), "condition"].tolist())]
#         elif wc.reps == "pooled":
#             # we want pooled results per each sample with more replicates
#             samples = [x for x in set(cfg.loc[(cfg.tag != "") & (cfg.condition != "control"), "condition"].tolist())]
#         else:
#             # we want particular results of particular replicate of samples
#             samples = cfg.loc[(cfg.tag == wc.reps) & (cfg.condition != "control"), "condition"].tolist()
#         inputs['mtx'] = [ADIR+"/peaks_profile/over_"+wc.gene_set+"/"+i+"/"+wc.reps+"/coverage_matrix."+wc.dups+".mtx.gz" for i in set(samples)]
#     else: 
#         # wc.name contains condition
#         if wc.reps == "all":
#             tags = cfg.loc[cfg.condition == wc.name, "tag"].tolist()
#             if len(set(tags)) > 1:
#                 tags = tags+["pooled"]
#             files = []
#             for tag in set(tags):
#                 files.append(ADIR+"/peaks_profile/over_"+wc.gene_set+"/"+wc.name+"/"+tag+"/coverage_matrix."+wc.dups+".mtx.gz")
#             inputs['mtx'] = files
#         else:
#             # wc.reps is one of: pooled, norep, repX
#             inputs['mtx'] = [ADIR+"/peaks_profile/over_"+wc.gene_set+"/"+wc.name+"/"+wc.reps+"/coverage_matrix."+wc.dups+".mtx.gz"]
#     return inputs

rule plot_profile_and_heatmap:
    # input:  unpack(plot_profile_and_heatmap_input),
    input:  mtx = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.mtx.gz",
    output: heatmap = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.heatmap.pdf",
            profile = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.profile.pdf",
            combined= "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.combined.pdf",
    log:    run = "logs/{name}/plot_profile_and_heatmap.from_{tool}.{dups}.{filt}.over_{gs}.log",
    threads: 1
    params: dpi = 300,
            title = "{name}.{dups}.{filt}",
            profile_type = config["profile_type"],
            profile_avrg = config["profile_avrg"],
            mtx     = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.mtx.tmp.gz",
            data    = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.profile.data",
            heatmtx = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.heatmap_mtx.gz",
    conda:  "../wrappers/plot_profile_and_heatmap/env.yaml"
    script: "../wrappers/plot_profile_and_heatmap/script.py"


# rule compute_matrix_relative:
#     input:  mtx = ADIR+"/peaks_profile/over_{gene_set}/{name}/{reps}/coverage_matrix.{dups}.mtx.gz",
#     output: mtx = ADIR+"/peaks_profile/over_{gene_set}/{name}/{reps}/coverage_matrix.{dups}.rel_counts.mtx.gz",
#     log:    run = ADIR+"/logs/over_{gene_set}/compute_matrix_relative.{name}.{reps}.{dups}.log",
#     threads: 1
#     params: smooth = cfg.rel_smooth.min(),
#             rscript= workflow.basedir+"/../wrappers/compute_matrix_relative/computeMatrix_relative_output.R",
#             mtx = ADIR+"/peaks_profile/over_{gene_set}/{name}/{reps}/coverage_matrix.{dups}.rel_counts.mtx",
#     conda:  "../wrappers/compute_matrix_relative/env.yaml"
#     script: "../wrappers/compute_matrix_relative/script.py"
# 
# 
# def compute_matrix_input(wc):
#     if wc.reps == "norep":
#         name = wc.name
#         return "results/MACS_peaks/"+name+"/"+name+"."+wc.dups+".bigWig"
#     else:
#         name = wc.name+"_"+wc.reps
#         return "results/MACS_peaks/"+name+"/"+name+"."+wc.dups+".bigWig"
# 
rule compute_matrix:
    # input:  bwg = compute_matrix_input,
    input:  bwg = "results/{tool}_peaks/{name}/{name}.{dups}.bigWig",
            ref = "gene_sets/{gs}.gene_set.bed",
    output: mtx = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.mtx.gz",
    log:    run = "logs/{name}/compute_matrix.from_{tool}.{dups}.{filt}.over_{gs}.log",
    threads: 2,
    # wildcard_constraints: dups="[^\.]+",
    params: before = config['before_genes'],
            after = config['after_genes'],
            sample_name = "{name}.{dups}.{filt}",
            mtx = "results/peaks_QC/peak_profiles/over_{gs}/{name}.{dups}.from_{tool}.{filt}.over_{gs}.mtx.tmp.gz",
            rscript = workflow.basedir+"/wrappers/compute_matrix/complete_regions_set.R",
    conda:  "../wrappers/compute_matrix/env.yaml"
    script: "../wrappers/compute_matrix/script.py"


def create_gene_set_inputs(wc):
    inputs = {
      "ref": default_reference
    }
    if not ',' in gene_sets[wc.gene_set]:
        if wc.gene_set != "all_genes":
            inputs['gset'] = "gene_sets/"+gene_sets[wc.gene_set]
        else:
            inputs['gset'] = default_reference
    return inputs
  
# This rule converts GTF file into bed12 gene_set if needed or just copy bed file into project
rule create_gene_set:
    # input:  ref = default_reference
    input:  unpack(create_gene_set_inputs)
    output: ref = "gene_sets/{gene_set}.gene_set.bed"
    log:    run = "logs/create_gene_set.{gene_set}.log",
    threads: 1
    params: gene_sets = gene_sets,
            list = "gene_sets/{gene_set}.gene_set.list"
    conda:  "../wrappers/create_gene_set/env.yaml"
    script: "../wrappers/create_gene_set/script.py"


rule plot_selected_peak_profiles:
    input:  bwg = "results/MACS_peaks/{name}/{name}.{dups}.bigWig",
            bed = "results/MACS_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv",
    output: plot = "results/peaks_QC/top_{top}_peak_profiles_by_score/{name}.{dups}.top_peak_profiles.pdf",
    log:    run  = "logs/{name}/plot_{top}_selected_peak_profiles.{dups}.log",
    threads: 1,
    params: before = config['top_peaks_neighborhood'],
            after = config['top_peaks_neighborhood'],
            dpi = 300,
            profile_type = "fill",
            profile_avrg = "median",
            prefix = "results/peaks_QC/top_{top}_peak_profiles_by_score/{name}.{dups}/",
            bed = "results/peaks_QC/top_{top}_peak_profiles_by_score/{name}.{dups}/selected_peaks.bed",
            mtx = "results/peaks_QC/top_{top}_peak_profiles_by_score/{name}.{dups}/coverage.mtx.tmp.gz",
            sample_name = "{name}.{dups}",
    conda:  "../wrappers/plot_selected_peak_profiles/env.yaml"
    script: "../wrappers/plot_selected_peak_profiles/script.py"


def plot_average_peak_profile_inputs(wc):
    inputs= dict()
    inputs['bwg'] = f"results/{wc.tool}_peaks/{wc.name}/{wc.name}.{wc.dups}.bigWig"
    if wc.filt == "all":
      inputs['bed'] = f"results/{wc.tool}_peaks/{wc.name}/{wc.name}.{wc.dups}.summits.all.bed"
    else:
      inputs['bed'] = f"results/{wc.tool}_peaks/{wc.name}/{wc.name}.{wc.dups}.summits.bed"
    return inputs

rule plot_average_peak_profile:
    input:  unpack(plot_average_peak_profile_inputs),
    output: plot= "results/peaks_QC/peak_profiles/over_peaks/{name}.{dups}.from_{tool}.{filt}.average_peak_profile.pdf",
            profile="results/peaks_QC/peak_profiles/over_peaks/{name}.{dups}.from_{tool}.{filt}.average_peak_profile_only.pdf",
            profile_data = "results/peaks_QC/peak_profiles/over_peaks/{name}.{dups}.from_{tool}.{filt}.average_peak_profile.data",
    log:    run  = "logs/{name}/plot_average_peak_profile.from_{tool}.{dups}.{filt}.log",
    threads: 2,
    params: before = config['peak_profiles_radius'],
            after = config['peak_profiles_radius'],
            dpi = 300,
            title = "{name}.{dups}",
            profile_type = "std",
            profile_avrg = "median",
            sample_name = "{name}.{dups}.{filt}",
            mtx = "results/peaks_QC/peak_profiles/over_peaks/{name}.{dups}.from_{tool}.{filt}.average_peak_profile.mtx.tmp.gz",
            data = "results/peaks_QC/peak_profiles/over_peaks/{name}.{dups}.from_{tool}.{filt}.average_peak_profile.heatmap_matrix.gz",
    conda:  "../wrappers/plot_average_peak_profile/env.yaml"
    script: "../wrappers/plot_average_peak_profile/script.py"
  
  
def plot_FRiP_inputs(wc):
  inputs = dict()
  suffix = ".pseudo_reps"
  if wc.name.endswith(suffix):
    samples = sample_tab.loc[sample_tab.name == wc.name[:-len(suffix)].split('_VS_')[0], 'sample_name'].unique()
  else:
    samples = sample_tab.loc[sample_tab.name == wc.name.split('_VS_')[0], 'sample_name'].unique()
    if not '_VS_' in wc.name:
      samples = sample_tab.loc[sample_tab['condition'] == wc.name, "sample_name"].unique()
  inputs['bam'] = expand("mapped/{sample}.{dups}.bam", sample=samples, dups=wc.dups)
  if wc.filt == "all":
    inputs['bed'] = f"results/{wc.tool}_peaks/{wc.name}/{wc.name}.{wc.dups}.peaks.all.narrowPeak"
  else:
    inputs['bed'] = f"results/{wc.tool}_peaks/{wc.name}/{wc.name}.{wc.dups}.peaks.narrowPeak"
  return inputs
    
rule plot_FRiP:
    input:  unpack(plot_FRiP_inputs)
    output: plot = "results/peaks_QC/FRiP/{name}.{dups}.from_{tool}.{filt}.FRiP.pdf",
            counts = "results/peaks_QC/FRiP/{name}.{dups}.from_{tool}.{filt}.FRiP.raw_counts",
    log:    run  = "logs/{name}/plot_FRiP.from_{tool}.{dups}.{filt}.log",
    conda:  "../wrappers/plot_FRiP/env.yaml"
    script: "../wrappers/plot_FRiP/script.py"


def call_chipqc_inputs(wc):
    inputs = dict()
    inputs["peaks"] = f"results/MACS_peaks/{wc.sample}/{wc.sample}.{wc.dups}.peaks.all.narrowPeak"
    inputs["reads"] = f"mapped/{sample_tab.loc[sample_tab.peaks_name == wc.sample, 'sample_name'].unique()[0]}.{wc.dups}.bam"
    return inputs

rule call_chipqc:
    input:  unpack(call_chipqc_inputs),
    output: html = "results/ChIPQC/{sample}/{sample}.{dups}.report.html",
            Robj = "results/ChIPQC/{sample}/{sample}.{dups}.report.RData",
            Rsam = "results/ChIPQC/{sample}/{sample}.{dups}.sample.RData",
    log:    run = "logs/{sample}/call_chipqc.{dups}.log"
    params: rscript = workflow.basedir+"/wrappers/call_chipqc/chipqc_sample.R",
            odir = "results/ChIPQC/{sample}",
            prefix="{sample}.{dups}.report"
    conda:  "../wrappers/call_chipqc/env.yaml"
    script: "../wrappers/call_chipqc/script.py"
    
#########################################
# Examine Peaks reproducibility
#

def reproducible_peaks_summary_inputs(wc):
    inputs = list()
    doublerep = sample_tab.loc[(sample_tab.is_control==False) & (sample_tab.num_of_reps==2), "condition"].unique()
    inputs+= expand("results/IDR_peaks/{sample}/{sample}.{reps}.{dups}.peaks.all.bed", 
                sample=doublerep, 
                dups=wc.dups, 
                reps=['true_reps','pseudo_reps'])
    triplerep = sample_tab.loc[(sample_tab.is_control==False) & (sample_tab.num_of_reps==3), "condition"].unique()
    inputs+= expand("results/IDR_peaks/{sample}/{sample}.{reps}.{comp}.{dups}.peaks.all.bed", 
                sample=triplerep, 
                dups=wc.dups, 
                comp=["rep1_VS_rep2","rep2_VS_rep3","rep1_VS_rep3"],
                reps=['true_reps','pseudo_reps'])
    multirep = sample_tab.loc[(sample_tab.is_control==False) & (sample_tab.num_of_reps>=2), "condition"].unique()
    inputs+= expand("results/ChIP-R_peaks/{sample}/{sample}.{reps}.{dups}.peaks.all.bed", 
                sample=multirep, 
                reps=['true_reps','pseudo_reps'],
                dups=wc.dups)
    return inputs

rule reproducible_peaks_summary:
    input:  bed = reproducible_peaks_summary_inputs
    output: tab = "results/reproducible_peaks_summary.{dups}.tsv"
    log:    run = "logs/reproducible_peaks_summary.{dups}.log",
    params: cutof = config["idr_cutof"],
    conda:  "../wrappers/reproducible_peaks_summary/env.yaml"
    script: "../wrappers/reproducible_peaks_summary/script.R"

    
def call_IDR_inputs(wc):
    inputs = dict()
    # This branching decides if condition has only two replicates or three (no other options are allowed)
    if "_VS_" in wc.reps:
      # This case sets sample list with three replicates
      reps = wc.reps.split(".")[0]
      tags = wc.reps.split(".")[1].split("_VS_")
      if not tags:
        print("WARNING: something wrong happened when extracting replicates from wildcard reps in call_IDR")
      samples = sample_tab.loc[(sample_tab.condition==wc.sample) & (sample_tab.tag.isin(tags)),'peaks_name'].unique()
    else:
      # This case sets sample list with only two replicates
      reps = wc.reps
      samples = sample_tab.loc[sample_tab.condition==wc.sample,'peaks_name'].unique()
    # This branching decides about types of replicates to be compared (so far only true or pseudo are allowed)
    if wc.reps.startswith("true_reps"):
      inputs['peaks'] = expand("results/MACS_peaks/{cond}/{cond}.{dups}.peaks.all.narrowPeak", dups=wc.dups, cond=samples)
    elif wc.reps.startswith("pseudo_reps"):
      inputs['peaks'] = expand("results/MACS_peaks/{cond}.{reps}/{cond}.{reps}.{dups}.peaks.all.narrowPeak", dups=wc.dups, reps=reps, cond=samples)
    else:
      print(f"ERROR: disallowed value for wildcard reps: {wc.reps}")
    return inputs

rule call_IDR:
    input:  unpack(call_IDR_inputs)
    output: bed = "results/IDR_peaks/{sample}/{sample}.{reps}.{dups}.peaks.bed",
            all_bed = "results/IDR_peaks/{sample}/{sample}.{reps}.{dups}.peaks.all.bed",
    log:    run = "logs/{sample}/call_IDR_for_{reps}.{dups}.log",
    threads:  1
    params: tmpd = GLOBAL_TMPD_PATH,
            cutof = config['idr_cutof'],
    conda:  "../wrappers/call_IDR/env.yaml"
    script: "../wrappers/call_IDR/script.py"
    

def call_chipr_inputs(wc):
    inputs = dict()
    samples = sample_tab.loc[sample_tab.condition==wc.cond,'peaks_name'].unique()
    # This branching decides about types of replicates to be compared (so far only true or pseudo are allowed)
    if wc.reps.startswith("true_reps"):
      inputs['peaks'] = expand("results/MACS_peaks/{name}/{name}.{dups}.peaks.all.narrowPeak", dups=wc.dups, name=samples)
    elif wc.reps.startswith("pseudo_reps"):
      inputs['peaks'] = expand("results/MACS_peaks/{name}.{reps}/{name}.{reps}.{dups}.peaks.all.narrowPeak", dups=wc.dups, reps=wc.reps, name=samples)
    else:
      print(f"ERROR: disallowed value for wildcard reps: {wc.reps}")
    return inputs

rule call_chipr:
    input:  unpack(call_chipr_inputs)
    output: bed = "results/ChIP-R_peaks/{cond}/{cond}.{reps}.{dups}.peaks.all.bed"
    log:    run = "logs/{cond}/call_chipr_for_{reps}.{dups}.log"
    params: minentries = 1,
            minsize = 20,
            cutof = config['idr_cutof'],
            prefix = "results/ChIP-R_peaks/{cond}/{cond}.{reps}.{dups}.peaks",
    conda:  "../wrappers/call_chipr/env.yaml"
    script: "../wrappers/call_chipr/script.py"
    
#########################################
# Peak calling
#

# rule calculate_score_tracks:
#     input:  trt = "results/MACS_peaks/{name}/{name}.{dups}.bdg",
#             ctl = "results/MACS_peaks/{name}/{name}.{dups}.control.bdg",
#     output: ppois = "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.p-value_score.bdg",
#             qpois = "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.q-value_score.bdg",
#             maxval= "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.maximum_value_score.bdg",
#             subtr = "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.subtract_from_treatment_score.bdg",
#             FE    = "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.linear_fold_enrichment_score.bdg",
#             logFE = "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.log10_fold_enrichment_score.bdg",
#             logLR = "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.log10_likelihood_ratio_score.bdg",
#             # slogLR= "results/MACS_peaks/{name}/score_tracks/{name}.{dups}.symmetric_log10_likelihood_ratio_score.bdg", # probably not usefull, as it should be run on two ChiP enriched tracks
#     log:    run = "logs/{name}/calculate_score_tracks.{dups}.log",
#     params: prefix= "results/MACS_peaks/{name}/score_tracks/{name}.{dups}",
#             pseudo= 0.00001,
#             stats = "ppois qpois max subtract FE logFE logLR"
#     conda:  "../wrappers/calculate_score_tracks/env.yaml"
#     script: "../wrappers/calculate_score_tracks/script.py"
    

def annotate_peaks_inputs(wc):
  suffix = ".pseudo_reps"
  if wc.name.endswith(suffix):
    samples = [i for i in wc.name[:-len(suffix)].split('_VS_') if i != 'no_control']
  else:
    samples = [i for i in wc.name.split('_VS_') if i != 'no_control']
    if not '_VS_' in wc.name:
      samples = sample_tab.loc[sample_tab['condition'] == wc.name, "control"].unique()
      samples = [*samples, *sample_tab.loc[sample_tab['condition'] == wc.name, "name"].unique()]
  inputs = {
    'bed': f"results/{wc.tool}_peaks/{wc.name}/{wc.name}.{wc.dups}.peaks.all.narrowPeak",
    'gtf': expand(reference_directory+"/annot/{ref}.gtf", ref=config['reference'])[0],
    'fa' : expand(reference_directory+"/seq/{ref}.fa", ref=config['reference'])[0],
    'tagdir': expand("results/HOMER_tag_dirs/{sample}.{dups}", sample=samples, dups=wc.dups)
  }
  return inputs

rule annotate_peaks:
    input:  unpack(annotate_peaks_inputs)
    output: tsv = "results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.annotated.tsv",
            annstats = "results/{tool}_peaks/{name}/{name}.{dups}.peaks.all.ann_stats.tsv",
    log:    run = "logs/{name}/annotate_{tool}_peaks.{dups}.log",
    resources: mem = 5
    params: rscript = workflow.basedir+"/wrappers/annotate_peaks/plots_and_stats.R",
            fdr_cutof = config["macs_padj_filter"],
            best = config["top_peaks"],
            tmpd = GLOBAL_TMPD_PATH,
    conda:  "../wrappers/annotate_peaks/env.yaml"
    script: "../wrappers/annotate_peaks/script.py"
    

def create_tag_dir_by_HOMER_input(wc):
    sample_names = sample_tab.loc[sample_tab['name'] == wc.name, "sample_name"].unique()
    inputs = {
      'bam': expand("mapped/{sn}.{d}.bam", sn=sample_names, d=wc.dups),
      'fa' : expand(reference_directory+"/seq/{ref}.fa", ref=config['reference'])[0],
    }
    return inputs

rule create_tag_dir_by_HOMER:
    input:  unpack(create_tag_dir_by_HOMER_input),
    output: tagdir = directory("results/HOMER_tag_dirs/{name}.{dups}"),
            ok = "results/HOMER_tag_dirs/{name}.{dups}/all_finished",
    log:    run = "logs/{name}.{dups}/create_tag_dir_by_HOMER.log",
    resources: mem = 5,
    conda:  "../wrappers/create_tag_dir_by_HOMER/env.yaml"
    script: "../wrappers/create_tag_dir_by_HOMER/script.py"


# def merge_replicates_inputs(wc):
#     samples = sample_tab.loc[sample_tab.condition==wc.cond,'name'].unique()
#     print(samples)
#     inputs = {
#       'peaks': expand("results/MACS_peaks/{name}/{name}.{dups}.consensus_peaks.bed", name=samples, dups=wc.dups)
#     }
#     return inputs
#
# rule merge_replicates:
#     # input:  unpack(merge_replicates_inputs)
#     input:  bed = lambda wc: expand("results/MACS_peaks/{cond}/{cond}.{{dups}}.peaks.all.narrowPeak", cond=sample_tab.loc[sample_tab.condition==wc.cond,'peaks_name'].unique()),
#     output: bed = "results/Fisher_peaks/{cond}/{cond}.{dups}.peaks.all.narrowPeak"
#     log:    run = "logs/{cond}.{dups}/merge_replicates.log",
#     conda:  "../wrappers/merge_replicates/env.yaml"
#     script: "../wrappers/merge_replicates/script.py"
    
# https://github.com/ENCODE-DCC/chip-seq-pipeline2

rule call_MSPC:
    input:  bed = lambda wc: expand("results/MACS_peaks/{cond}/{cond}.{{dups}}.peaks.all.narrowPeak", cond=sample_tab.loc[sample_tab.condition==wc.cond,'peaks_name'].unique()),
    output: bed = "results/MSPC_peaks/{cond}/{cond}.{dups}.peaks.all.narrowPeak",
    log:    run = "logs/{cond}/call_MSPC.{dups}.log",
    threads: 10
    params: rscript = workflow.basedir+"/wrappers/call_MSPC/call_MSPC.R",
            rep_type = "Bio",
            strong = 1e-8,
            weak = 1e-4,
            fdr = 1,
            bed = "results/MSPC_peaks/{cond}/{cond}.{dups}.consensusPeaks.txt"
    conda:  "../wrappers/call_MSPC/env.yaml"
    script: "../wrappers/call_MSPC/script.py"


def call_SEACR_inputs(wc):
    inputs = {"ref": expand(reference_directory+"/seq/{ref}.chrom.sizes", ref=config['reference'])[0]}
    dups = wc.dups
    samples = wc.name.split('_VS_')
    # This case is for doing peak calling on true replicates
    if len(samples) == 2:
      # This case is for using specific replicate to call SEACR (with or without control)
      inputs["trt"] = expand("mapped/{sample}.{dups}.bedgraph", sample=sample_tab.loc[sample_tab.name == samples[0], "sample_name"].unique(), dups=dups)[0]
      if samples[1] != 'no_control':
        inputs['ctl'] = expand("mapped/{sample}.{dups}.bedgraph", sample=sample_tab.loc[sample_tab.name == samples[1], 'sample_name'].unique(), dups=dups)[0]
    else:
      # This case is for using all reps of one condition to call MACS (with or without control)
      print("ERROR: SEACR cannot process multi-replicate files!")
      #inputs["trt"] = expand("mapped/{sample}.{dups}.bedgraph", sample=sample_tab.loc[sample_tab.condition == samples[0], "sample_name"].unique(), dups=dups)
      #controls = sample_tab.loc[sample_tab.condition == samples[0], "control"].unique()
      #controls = sample_tab.loc[[i in controls for i in sample_tab.name], "sample_name"].unique()
      #if any(controls):
      #  inputs['ctl'] = expand("mapped/{sample}.{dups}.bedgraph", sample=controls, dups=dups)
    return(inputs)

rule call_SEACR:
    input: unpack(call_SEACR_inputs),
    output: tab = "results/SEACR_peaks/{name}/{name}.{dups}.peaks.narrowPeak",
            tab_all = "results/SEACR_peaks/{name}/{name}.{dups}.peaks.all.narrowPeak",
            bdg = "results/SEACR_peaks/{name}/{name}.{dups}.peaks.bedgraph",
            bdg_all = "results/SEACR_peaks/{name}/{name}.{dups}.peaks.all.bedgraph",
            bwg = "results/SEACR_peaks/{name}/{name}.{dups}.peaks.bigWig",
            bwg_all = "results/SEACR_peaks/{name}/{name}.{dups}.peaks.all.bigWig",
    log: "logs/{name}/call_SEACR.{dups}.log",
    threads: 1
    params: prefix = "results/SEACR_peaks/{name}/{name}.{dups}",
            norm = config["seacr_normalisation"],
            ctl_thr = config["seacr_fdr"],
            temp = GLOBAL_TMPD_PATH,
    conda:  "../wrappers/call_SEACR/env.yaml"
    script: "../wrappers/call_SEACR/script.py"
    

def call_macs2_inputs(wc):
    inputs = {"ref": expand(reference_directory+"/seq/{ref}.chrom.sizes", ref=config['reference'])[0]}
    dups = wc.dups
    samples = wc.name.split('_VS_')
    suffix = '.pseudo_reps'
    if wc.name.endswith(suffix):
      # This case is for doing peak calling on specific pseudo replicates (with or without control)
      samples[1] = samples[1][:-len(suffix)]
      name = sample_tab.loc[sample_tab.name == samples[0], ["condition","tag"]]
      inputs["trt"] = expand("mapped/pseudo/{cond}_{rep}.{dups}.bam", cond=name['condition'].unique(), rep=name['tag'].unique(), dups=dups)
      if samples[1] != 'no_control':
        name = sample_tab.loc[sample_tab.name == samples[1], ["condition","tag"]]
        if len(name["tag"].unique()) == 1:
          # This case is for having only one repllicate of DNA innput sample
          if name['tag'].unique() == "":
            # Here, DNA input sample doesn't use any tag (e.g. Replicate column is empty)
            inputs['ctl'] = expand("mapped/{cond}.{dups}.bam", cond=name['condition'].unique(), dups=dups)
          else:
            # Here, DNA input sample does use a tag (e.g., Replicate column is rep1)
            inputs['ctl'] = expand("mapped/{cond}_{rep}.{dups}.bam", cond=name['condition'].unique(), rep=name['tag'].unique(), dups=dups)
        else:
          inputs['ctl'] = expand("mapped/pseudo/{cond}_{rep}.{dups}.bam", cond=name['condition'].unique(), rep=name['tag'].unique(), dups=dups)
    else:
      # This case is for doing peak calling on true replicates
      if len(samples) == 2:
        # This case is for using specific replicate to call MACS (with or without control)
        inputs["trt"] = expand("mapped/{sample}.{dups}.bam", sample=sample_tab.loc[sample_tab.name == samples[0], "sample_name"].unique(), dups=dups)
        if samples[1] != 'no_control':
          inputs['ctl'] = expand("mapped/{sample}.{dups}.bam", sample=sample_tab.loc[sample_tab.name == samples[1], 'sample_name'].unique(), dups=dups)
      else:
        # This case is for using all reps of one condition to call MACS (with or without control)
        inputs["trt"] = expand("mapped/{sample}.{dups}.bam", sample=sample_tab.loc[sample_tab.condition == samples[0], "sample_name"].unique(), dups=dups)
        controls = sample_tab.loc[sample_tab.condition == samples[0], "control"].unique()
        controls = sample_tab.loc[[i in controls for i in sample_tab.name], "sample_name"].unique()
        if any(controls):
          inputs['ctl'] = expand("mapped/{sample}.{dups}.bam", sample=controls, dups=dups)
    return inputs
        
#TODO: zakomponovat pouziti faCounts skriptu na vypocitani effective genome size (cize pocet baz (ACGT) minus pocet N)
#TODO: zakomponovat pouziti broad_peaks nasavenia
rule call_macs2:
    input:  unpack(call_macs2_inputs),
    output: trt_bdg = "results/MACS_peaks/{name}/{name}.{dups}.bdg",
            trt_bwg = "results/MACS_peaks/{name}/{name}.{dups}.bigWig",
            ctl_bdg = "results/MACS_peaks/{name}/{name}.{dups}.control.bdg",
            ctl_bwg = "results/MACS_peaks/{name}/{name}.{dups}.control.bigWig",
            sum_tab = "results/MACS_peaks/{name}/{name}.{dups}.summits.bed",
            nar_tab = "results/MACS_peaks/{name}/{name}.{dups}.peaks.narrowPeak",
            xls_tab_all = "results/MACS_peaks/{name}/{name}.{dups}.peaks.all.xls",
            sum_tab_all = "results/MACS_peaks/{name}/{name}.{dups}.summits.all.bed",
            nar_tab_all = "results/MACS_peaks/{name}/{name}.{dups}.peaks.all.narrowPeak",
    log:    run = "logs/{name}/call_macs2.{dups}.log",
    threads: 1
    params: trt_bdg = "results/MACS_peaks/{name}/{name}.{dups}_treat_pileup.bdg",
            trt_bwg = "results/MACS_peaks/{name}/{name}.{dups}_treat_pileup.bigWig",
            ctl_bdg = "results/MACS_peaks/{name}/{name}.{dups}_control_lambda.bdg",
            ctl_bwg = "results/MACS_peaks/{name}/{name}.{dups}_control_lambda.bigWig",
            xls_tab = "results/MACS_peaks/{name}/{name}.{dups}_peaks.xls",
            sum_tab = "results/MACS_peaks/{name}/{name}.{dups}_summits.bed",
            nar_tab = "results/MACS_peaks/{name}/{name}.{dups}_peaks.narrowPeak",
            effective_GS = config["effective_genome_size"],
            frag_len = config["fragment_length"],
            qval_cutof = config["macs_padj_filter"],
            broad = config["macs_broad_peaks"], # TODO: doriesit
            brcut = config["macs_broad_cutof"],
            dir = "results/MACS_peaks/{name}/",
            name= "{name}.{dups}",
            temp= GLOBAL_TMPD_PATH,
    conda:  "../wrappers/call_macs2/env.yaml"
    script: "../wrappers/call_macs2/script.py"
    
##########################################
# PREPARE INPUT
#

rule convert_bam_to_bedgraph:
    input:  bam = "{file_name}.bam",
            ref = expand(reference_directory+"/seq/{ref}.chrom.sizes", ref=config['reference'])[0],
    output: bdg = "{file_name}.bedgraph",
    log:    "logs/{file_name}.bam2bedgraph.log",
    threads: 2
    params: max_frag_len = "1000",
            temp= GLOBAL_TMPD_PATH,
    conda:  "../wrappers/convert_bam_to_bedgraph/env.yaml"
    script: "../wrappers/convert_bam_to_bedgraph/script.py"


rule prepare_pseudo_reps:
    input:  bam = lambda wc: expand("mapped/{sample}.{{dups}}.bam", sample=sample_tab.loc[sample_tab.condition==wc.cond, 'sample_name'].unique())
    output: rep1 = "mapped/pseudo/{cond}_rep1.{dups}.bam",
            rep2 = "mapped/pseudo/{cond}_rep2.{dups}.bam",
            rep3 = "mapped/pseudo/{cond}_rep3.{dups}.bam",
            rep4 = "mapped/pseudo/{cond}_rep4.{dups}.bam",
    log:    run = "logs/{cond}/prepare_pseudo_reps.{dups}.log",
    threads: 1
    resources:  mem = 25,
    params: num_of_reps = lambda wc: sample_tab.loc[sample_tab.condition==wc.cond, 'num_of_reps'].unique()[0],
            tags = lambda wc: sample_tab.loc[sample_tab.condition==wc.cond, 'tag'].unique(),
            merged = "mapped/pseudo/{cond}_merged.{dups}.bam",
            header = "mapped/pseudo/{cond}_merged.{dups}.header",
            prefix = "mapped/pseudo/{cond}_rep",
    conda:  "../wrappers/prepare_pseudo_reps/env.yaml"
    script: "../wrappers/prepare_pseudo_reps/script.py"

    
