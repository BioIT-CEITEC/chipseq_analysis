#########################################
# wrapper for rule: call_macs2
#########################################
import os
import sys
import math
import subprocess
import re
from snakemake.shell import shell

f = open(snakemake.log.run, 'a+')
f.write("\n##\n## RULE: call_macs2 \n##\n")
f.close()

shell.executable("/bin/bash")

version = str(subprocess.Popen("conda list 2>&1 ", shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
f = open(snakemake.log.run, 'at')
f.write("## CONDA:\n"+version+"\n")
f.close()

input_line = "-t "+" ".join(snakemake.input.trt)
inputs = list(snakemake.input.trt)
if hasattr(snakemake.input, 'ctl'):
  input_line += " -c "+" ".join(snakemake.input.ctl)
  inputs += list(snakemake.input.ctl)

# First we need to check if input files are paired-end (by default) or single-end
paired = True
for inp in inputs:
    command = 'samtools view '+inp+' 2>> '+snakemake.log.run+' | head -1 2>> '+snakemake.log.run+' | cut -f 2 2>> '+snakemake.log.run
    # command = 'bc <<< "$(samtools view '+inp+' | head -1 | cut -f 2) % 2"'
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    flag = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    f.write("## FLAG:"+flag)
    f.write("## INFO: file "+inp+" is "+("paired-end" if int(flag)%2==1 else "single-end")+"\n\n")
    f.close()
    if int(flag)%2==0:
        paired = False
        break

if snakemake.params.spikein:
  # In case of spike-in normalisation the input files need to be normalised and converted into bedgraph explicitly
  spike_inputs = list(snakemake.input.trt_spike)
  if hasattr(snakemake.input, 'ctl_spike'):
    spike_inputs += list(snakemake.input.ctl_spike)

  for i in range(len(spike_inputs)):
    bam = inputs[i]
    bdg = os.path.join(snakemake.params.dir, os.path.basename(bam).replace('.bam','.bedgraph'))
    sbam = spike_inputs[i]
    # Compute spike-in normalisation factor
    if paired:
      ## Flag 2816 is combination of 1) not primary alignment, 2) read fails platform/vendor quality checks and 3) supplementary alignment; which we don't want
      command = "$(which time) --verbose samtools view -c -f 66 -F 2816 -@ "+str(snakemake.threads)+" "+sbam+" 2>> "+snakemake.log.run
    else:
      command = "$(which time) --verbose samtools view -c -F 2816 -@ "+str(snakemake.threads)+" "+sbam+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    spike_frags = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
    f.write("## INFO: Spike-in fragments: "+str(spike_frags))
    scaling_spikein = round(float(snakemake.params.scalefac)/int(spike_frags), 8)
    f.write("## INFO: Spike-in scale factor: "+str(scaling_spikein)+"\n\n")
    f.close()

    dlen = snakemake.params.frag_len
    if str(dlen) == 'unk':
      command = "$(which time) --verbose macs2 predictd"+\
                " -i "+bam+\
                " -g "+snakemake.params.effective_GS+\
                " 2>&1 | tee -a "+snakemake.log.run
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.flush()
      predictd_out = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
      for line in predictd_out.split('\n'):
        if '# predicted fragment length is' in line:
          dlen = re.findall('^.*# predicted fragment length is (-?[0-9]+).*$',line)[0]
      f.write("## INFO: predicted fragment length (d-length) is:"+dlen+"\n\n")
      f.close()
    else:
      dlen = int(float(dlen))
      f = open(snakemake.log.run, 'at')
      f.write("## INFO: explicit fragment length (d-length) is: "+str(dlen)+"\n\n")
      f.close()

    # Converting BAM file into BED file containing only reads properly aligned as primary (and paired, if possible)
    if paired:
      # TODO: consider running `macs2 preditd` or using config['fragment_length'] to extend the insert size if shorter
      command = "$(which time) --verbose samtools view -uh -f 2 -F 2816 -@ "+str(snakemake.threads)+" "+bam+" 2>> "+snakemake.log.run+\
                " | $(which time) --verbose samtools sort -n -@ "+str(snakemake.threads)+" 2>> "+snakemake.log.run+\
                " | $(which time) --verbose bedtools bamtobed -bedpe -i stdin 2>> "+snakemake.log.run+\
                " | awk '$1==$4' 2>> "+snakemake.log.run+\
                " | cut -f 1,2,6 | $(which time) --verbose sort -k1,1 -k2,2n -k3,3n 2>> "+snakemake.log.run+" > "+snakemake.params.bed
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)
    else:
      # TODO: Here should be an extraction of single-end bed file from BAM following by extension using macs2 pileup (https://github.com/macs3-project/MACS/wiki/Advanced:-Call-peaks-using-MACS2-subcommands#step-3-extend-chip-sample-to-get-chip-coverage-track)
      #       use bedtools bamtobed to extract filtered reads in BED format, then if it's IP sample, extend the SE read to the d (fragment) length on 3'end considering the strand or half-d length on both sides if it's control sample.
      #       The d length is taken either from config['fragment_length'] if it's a number or from `macs2 preditd` command if it's 'unk'.
      command = "$(which time) --verbose bedtools bamtobed -i "+bam+\
                " 2>> "+snakemake.log.run+\
                " | awk '{{ if($3-$2 < "+str(dlen)+") $3=$2+"+str(dlen)+"; print $0 }}' OFS='\t' 2>> "+snakemake.log.run+\
                " | $(which time) --verbose sort -k1,1 -k2,2n -k3,3n 2>> "+snakemake.log.run+" > "+snakemake.params.bed
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)

      exit("Not finished yet!")

    # Converting BED file into bedgraph track file for peak calling
    if bam in snakemake.input.trt:
      # This case is for converting treatment (ChIP) samples into bedgraph track file for peak calling
      command = "$(which time) --verbose bedtools genomecov -bga -scale "+str(scaling_spikein)+\
                " -i "+snakemake.params.bed+\
                " -g "+snakemake.input.ref+\
                " > "+bdg+" 2>> "+snakemake.log.run
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)

    elif bam in snakemake.input.ctl:
      # This case is for converting control samples into lambda bedgraph track file for peak calling (more complicated because of all lambda tracks)
      # Extraction of fragment-length background
      bdg_dlen = os.path.join(snakemake.params.dir, os.path.basename(bam).replace('.bam','.dlen_bg.bedgraph'))
      command = "$(which time) --verbose bedtools genomecov -bga -i "+snakemake.params.bed+\
                " -g "+snakemake.input.ref+\
                " > "+bdg_dlen+" 2>> "+snakemake.log.run
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)
      # Extraction of slocal background
      scaling_local = 0.2 # This value is important, because we extend every fragment in bed file five times of its original length, therefore, the pile-up must be scaled down to one fifth
      bdg_slocal = os.path.join(snakemake.params.dir, os.path.basename(bam).replace('.bam','.slocal_bg.bedgraph'))
      command = "$(which time) --verbose awk '{{ext=2*($3-$2); print $1, ($2-ext<0)?0:$2-ext, $3+ext}}' OFS='\t' "+snakemake.params.bed+" 2>> "+snakemake.log.run+\
                " | $(which time) --verbose sort -k1,1 -k2,2n -k3,3n 2>> "+snakemake.log.run+\
                " | $(which time) --verbose bedtools genomecov -bga -scale "+str(scaling_local)+" -i stdin -g "+snakemake.input.ref+" 2>> "+snakemake.log.run+\
                " > "+bdg_slocal
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)
      # Extraction of llocal background
      scaling_local = 0.04 # This value is important, because we extend every fragment in bed file five times of its original length, therefore, the pile-up must be scaled down to one fifth
      bdg_llocal = os.path.join(snakemake.params.dir, os.path.basename(bam).replace('.bam','.llocal_bg.bedgraph'))
      command = "$(which time) --verbose awk '{{ext=12*($3-$2); print $1, ($2-ext<0)?0:$2-ext, $3+ext}}' OFS='\t' "+snakemake.params.bed+" 2>> "+snakemake.log.run+\
                " | $(which time) --verbose sort -k1,1 -k2,2n -k3,3n 2>> "+snakemake.log.run+\
                " | $(which time) --verbose bedtools genomecov -bga -scale "+str(scaling_local)+" -i stdin -g "+snakemake.input.ref+" 2>> "+snakemake.log.run+\
                " > "+bdg_llocal
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)
      # Extraction of genome background (computing total fragment length)
      command = "awk '{{sum+=($3-$2)}}END{{print sum}}' "+snakemake.params.bed+" 2>> "+snakemake.log.run
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      total_len = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
      f.write("## INFO: total fragment length: "+total_len+"\n")
      # Extraction of genome background (computing total reference genome length)
      command = "awk '{{sum+=$2}}END{{print sum}}' "+snakemake.input.ref+" 2>> "+snakemake.log.run
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      genome_len = str(subprocess.Popen(command, shell=True, stdout=subprocess.PIPE).communicate()[0], 'utf-8')
      f.write("## INFO: total genome length: "+str(genome_len)+"\n")
      scaling_bg = round(float(total_len)/float(genome_len), 8)
      f.write("## INFO: genome background factor: "+str(total_len)+"/"+str(genome_len)+"="+str(scaling_bg)+"\n")
      f.close()
      # Combination of all backgrounds and normalisation
      command = "$(which time) --verbose macs2 bdgcmp -m max -t "+bdg_dlen+" -c "+bdg_slocal+" -o /dev/stdout 2>> "+snakemake.log.run+\
                " | $(which time) --verbose macs2 bdgcmp -m max -t /dev/stdin -c "+bdg_llocal+" -o /dev/stdout 2>> "+snakemake.log.run+\
                " | $(which time) --verbose macs2 bdgopt -m max -i /dev/stdin -p "+str(scaling_bg)+" -o /dev/stdout 2>> "+snakemake.log.run+\
                " | $(which time) --verbose macs2 bdgopt -i /dev/stdin -m multiply -p "+str(scaling_spikein)+" -o "+bdg+" 2>> "+snakemake.log.run
      f = open(snakemake.log.run, 'at')
      f.write("## COMMAND: "+command+"\n")
      f.close()
      shell(command)
    else:
      print("ERROR: file "+bam+" is neither in input trt nor in input ctl files!")
      print(exit)
      exit()
      
  if len(snakemake.input.trt) > 1:
    command = "$(which time) --verbose macs2 cmbreps -i "+\
              " ".join([os.path.join(snakemake.params.dir, os.path.basename(i).replace('.bam','.bedgraph')) for i in snakemake.input.trt])+\
              " -m max -o "+snakemake.output.trt_bdg+" 2>> "+snakemake.log.run
  else:
    command = "mv "+os.path.join(snakemake.params.dir, os.path.basename(snakemake.input.trt[0]).replace('.bam','.bedgraph'))+\
              " "+snakemake.output.trt_bdg+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  if len(snakemake.input.ctl) > 1:
    command = "$(which time) --verbose macs2 cmbreps -i "+\
              " ".join([os.path.join(snakemake.params.dir, os.path.basename(i).replace('.bam','.bedgraph')) for i in snakemake.input.ctl])+\
              " -m max -o "+snakemake.output.ctl_bdg+" 2>> "+snakemake.log.run
  else:
    command = "mv "+os.path.join(snakemake.params.dir, os.path.basename(snakemake.input.ctl[0]).replace('.bam','.bedgraph'))+\
              " "+snakemake.output.ctl_bdg+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  # TODO: add the branching for broadpeaks
  final_track = os.path.join(snakemake.params.dir, snakemake.params.name) + '.qpois.bdg'
  command = "$(which time) --verbose macs2 bdgcmp -t "+snakemake.output.trt_bdg+" -c "+snakemake.output.ctl_bdg+" -m qpois -o "+final_track+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  try:
    cutoff_stats = final_track.replace('.bdg', '.cutoff_stats')
    command = "$(which time) --verbose macs2 bdgpeakcall -i "+final_track+" --cutoff-analysis -c 1 -l 115 -g 75 -o "+cutoff_stats+" 2>> "+snakemake.log.run
    f = open(snakemake.log.run, 'at')
    f.write("## COMMAND: "+command+"\n")
    f.close()
    shell(command)
  except:
    f = open(snakemake.log.run, 'at')
    f.write("## INFO: Cutoff-analysis failed!\n")
    f.close()

  command = "$(which time) --verbose macs2 bdgpeakcall -i "+final_track+" -c 1 -l 115 -g 75 -o "+snakemake.params.nar_tab+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

  command = "cat "+snakemake.params.nar_tab+" | awk '{{if(NR==1){{print $0}}else{{$9=$5/10; print $0}}}}' OFS='\t' > "+snakemake.output.nar_tab_all+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

  # create empty output files
  command = "touch "+snakemake.output.sum_tab_all+" "+snakemake.output.sum_tab+" "+snakemake.output.xls_tab_all+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

else:
  # Without spike-in normalisation, macs2 could be used in single-command way with implicit normalisation and conversion of input files

  # Here should be a command using the faCount tool to count the real effective genome size if needed.
  if str(snakemake.params.effective_GS) == "unk":
    f = open(snakemake.log.run, 'at')
    f.write("## INFO: Effective genome size parameter (-g) can no longer be 'unk'! It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7) and 'dm' for fruitfly (1.2e8). Please, change the value. \n")
    f.close()
    exit()

  # Set the proper type of input data (single- vs. paired-end)
  if paired:
    input_line += " -f BAMPE"
  else:
    input_line += " -f BAM"

  keep_dups = "all"

  # Set the proper MACS2 parameters to estimate fragment length if unknown
  if snakemake.params.frag_len == "unk":
    nomodel = "--fix-bimodal"
  else:
    nomodel = "--nomodel --extsize "+str(snakemake.params.frag_len)

  command = "$(which time) macs2 callpeak "+input_line+\
            " --keep-dup "+keep_dups+\
            " -g "+str(snakemake.params.effective_GS)+\
            " --outdir "+snakemake.params.dir+\
            " --name "+snakemake.params.name+\
            " "+nomodel+\
            " --bdg"+\
            " --tempdir "+snakemake.params.temp+\
            " -q 0.1 >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

  # Final post-processing of resulting peaks (i.e., renaming, filtering and conversion)
  command = "mv "+snakemake.params.trt_bdg+" "+snakemake.output.trt_bdg+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  command = "mv "+snakemake.params.ctl_bdg+" "+snakemake.output.ctl_bdg+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  # rename original output file
  command = "mv "+snakemake.params.xls_tab+" "+snakemake.output.xls_tab_all+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)
  
  # rename original output file with no cutof
  command = "mv "+snakemake.params.sum_tab+" "+snakemake.output.sum_tab_all+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

  # use Qvalue cutof
  command = "$(which time) awk -v OFS='\t' -F '\t' '$5 > -log("+str(snakemake.params.qval_cutof)+")/log(10)' "+snakemake.output.sum_tab_all+" > "+snakemake.output.sum_tab+" 2>> "+snakemake.log.run
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

  # rename original output file with no cutof
  command = "mv "+snakemake.params.nar_tab+" "+snakemake.output.nar_tab_all+" >> "+snakemake.log.run+" 2>&1"
  f = open(snakemake.log.run, 'at')
  f.write("## COMMAND: "+command+"\n")
  f.close()
  shell(command)

# use Qvalue cutof
command = "$(which time) awk -v OFS='\t' -F '\t' '$9 > -log("+str(snakemake.params.qval_cutof)+")/log(10)' "+snakemake.output.nar_tab_all+" > "+snakemake.output.nar_tab+" 2>> "+snakemake.log.run
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

# bedGraphToBigWig ${i} /mnt/ssd/ssd_3/references/saccharomyces_cerevisiae/R64-1-1.100/seq/chrom.sizes ${i%.bdg}.bigWig
command = "$(which time) bedGraphToBigWig "+snakemake.output.trt_bdg+" "+snakemake.input.ref+" "+snakemake.output.trt_bwg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)

command = "$(which time) bedGraphToBigWig "+snakemake.output.ctl_bdg+" "+snakemake.input.ref+" "+snakemake.output.ctl_bwg+" >> "+snakemake.log.run+" 2>&1"
f = open(snakemake.log.run, 'at')
f.write("## COMMAND: "+command+"\n")
f.close()
shell(command)
