# Upstream steps (before `filter_calls.py`)

These steps produce `data/calls.tsv.gz` from the per-sample xTea VCFs. They are described in
Methods 2.4 of the paper:

1. **Coordinate harmonization.** Candidate insertions are matched across the three samples of a
   patient (T, NT, WB) using ±50-bp padded genomic intervals (interval tree).
2. **Blood subtraction.** Events detected in the matched whole blood are treated as constitutional
   and removed from that patient's T–NT comparison.
3. **T–NT comparison.** Remaining calls present in only one lung tissue are kept as T-specific or
   NT-specific calls.
4. **Allele-frequency filter.** Calls with AF ≤ 0.3 are kept.

<!-- TODO (authors): add the scripts used for steps 1–4 to this folder. -->

The output table must contain the columns listed in `data/README.md`. WB-specific calls are
included in the same table (compartment = WB) because the recurrent-locus filter considers calls
in any compartment of other patients.
