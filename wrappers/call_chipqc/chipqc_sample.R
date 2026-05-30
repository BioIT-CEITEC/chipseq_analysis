library(ChIPQC)
library(Rsamtools)
library(rtracklayer)
options("browser"="false")

args = commandArgs(trailingOnly=TRUE)
print(args)
outdir = args[1]
html_prefix = args[2]
sample_data = args[3]
peaks = args[4]
reads = args[5]

trim_bam_header <- function(bam_file, peak_file, out_bam) {

  # 1. Get chromosomes with actual reads (from BAM index)
  bam_index <- idxstatsBam(bam_file)
  chroms_with_reads <- bam_index$seqname[bam_index$mapped > 0]

  # 2. Get chromosomes present in peaks
  peaks <- rtracklayer::import(peak_file)
  chroms_with_peaks = as.character(unique(seqnames(peaks)))

  # 3. Union: keep contigs that have reads OR peaks
  chroms_to_keep = union(chroms_with_reads, chroms_with_peaks)

  # 4. Filter the BAM to kept chromosomes and reheader
  # filterBam writes only reads on kept chroms, and we rebuild seqinfo
  bam_header  = scanBamHeader(bam_file)
  old_targets = bam_header$mapped$targets                        # named int vector
  new_targets = old_targets[names(old_targets) %in% chroms_to_keep]

  # Replace targets in header
  bam_header$mapped$targets = new_targets

  # 5. Write filtered BAM via filterBam with a ScanBamParam restricted to
  #    the kept chromosomes, then fix the header with Rsamtools reheader
  param = ScanBamParam(
    which = as(Seqinfo(names(new_targets), new_targets), "GRanges")
  )

  tmp_bam = paste0(out_bam, ".tmp.bam")
  filterBam(bam_file, tmp_bam, param = param)

  # 6. Reheader: write a SAM header text file, call samtools reheader
  header_lines <- c(
    bam_header$text[["@HD"]],   # keep original @HD line
    paste0("@SQ\tSN:", names(new_targets), "\tLN:", new_targets),
    # keep @RG and @PG lines if present
    unlist(bam_header$text[names(bam_header$text) != "@HD" &
                            names(bam_header$text) != "@SQ"])
  )

  header_file <- tempfile(fileext = ".sam")
  writeLines(header_lines, header_file)

  system2("samtools", c("reheader", header_file, tmp_bam),
          stdout = out_bam)
  file.remove(tmp_bam)

  # 7. Index the output
  indexBam(out_bam)

  message("Done. Kept ", length(chroms_to_keep), " of ",
          length(old_targets), " chromosomes/contigs.")
  message("Output: ", out_bam)
}

new_reads = paste0(reads, ".new_header.bam")
trim_bam_header(reads, peaks, new_reads)

print("# Loading data into ChIPQCsample object")
sample = ChIPQCsample(new_reads, peaks)
print("# Saving ChIPQCsample object")
save(sample, file=sample_data)
print("# Producing ChIPQC report")
ChIPQCreport(sample, reportName=html_prefix, reportFolder=outdir)
