#!/usr/bin/env bash
# Extract mobile-element insertions from the gnomAD-SV v4.1 sites VCF (GRCh38).
# Download the sites VCF from https://gnomad.broadinstitute.org/downloads (gnomAD-SV v4.1).
# Output columns: chrom, start, end, type (e.g. INS:ME:ALU), AF  -> data/gnomad_v4.1_INS_ME_subset.tsv.gz
# Requires bcftools.
set -euo pipefail
VCF=${1:?usage: make_gnomad_subset.sh gnomad.v4.1.sv.sites.vcf.gz}
bcftools query -i 'INFO/SVTYPE="INS"' -f '%CHROM\t%POS\t%INFO/END\t%ALT\t%INFO/AF\n' "$VCF" \
  | sed -e 's/<//' -e 's/>//' \
  | awk -F'\t' '$4 ~ /^INS(:ME:(ALU|LINE1|SVA))?$/' \
  | gzip > data/gnomad_v4.1_INS_ME_subset.tsv.gz
