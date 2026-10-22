# Input data (not included)

Patient-level and sequencing data are not part of this repository. Place the following
tab-separated files here, or point `MEI_DATA_DIR` to their folder.

| File | Columns | Notes |
|---|---|---|
| `calls.tsv.gz` | patient_id, compartment (T/NT/WB), element (ALU/LINE1), chrom, pos, AF, LCLIP, RCLIP, LDISC, RDISC, SUBTYPE, REF_REP, GENE_INFO | Calls after blood subtraction and T–NT comparison, AF ≤ 0.3 (see `02_filtering/UPSTREAM_STEPS.md`). The read-support and annotation fields are taken directly from the xTea VCF INFO column. |
| `samples.tsv` | patient_id, compartment, coverage, mean_insert_size | One row per sequenced sample |
| `patients.tsv` | patient_id, site, age, histology, pack_years | `pack_years` may be empty |
| `wb_specific_counts.tsv` | patient_id, WB_ALU, WB_LINE1 | WB-specific AF ≤ 0.3 counts per patient (three-compartment comparison) |
| `exclusions.tsv` | patient_id, element, compartment, reason | Patient values excluded from specific analyses |
| `gnomad_v4.1_INS_ME_subset.tsv.gz` | chrom, start, end, type, AF (no header) | Made with `02_filtering/make_gnomad_subset.sh` |

The primary call set itself (93,161 candidate tissue-restricted insertions) is provided as
Supplementary Table S10 of the article.
