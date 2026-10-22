#!/usr/bin/env bash
# Step 1: mobile-element insertion calling with xTea v0.1.9 (https://github.com/parklab/xTea).
# xTea is third-party software and is not included here; install it following its own instructions
# and cite: Chu C, et al. Nat Commun 2021;12:3836.
#
# Exact command used in the study (Methods 2.3), run once per sample on GRCh38/hg38 BAM files:
#   sample_id.txt          one sample ID per line
#   illumina_bam_list.txt  sample ID and BAM path per line
#   rep_lib_annotation/    xTea repeat library for hg38
#   gencode.v33.annotation.gff3

xtea -i sample_id.txt -b illumina_bam_list.txt -x null -p ./path_work_folder/ \
  -o submit_jobs.sh -l rep_lib_annotation/ -r hg38.fa -g gencode.v33.annotation.gff3 \
  --xtea <xtea_lib_path> -y 1 --user --nclip 1 --cr 1 --nd 1 -f 5907 -n 16 -m 64 -x null --lsf

# Then submit the generated jobs (LSF in the study):  sh submit_jobs.sh
# Output: one VCF per sample and element family (Alu, LINE-1).
