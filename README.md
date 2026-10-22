# Somatic LINE-1 and Alu retrotransposition in lung tumors and adjacent normal lung

Analysis code for:

> Jacob-Hirsch J, Eyal E. *Whole-Genome Sequencing Reveals Somatic LINE-1 and Alu Retrotransposition
> Extending Beyond Lung Tumors into Adjacent Normal Tissue.* Genes 2026. DOI: [to be added]

The code starts from mobile-element insertion calls made with xTea and reproduces the call filtering,
detection-power model, statistics and figures reported in the article.

## Pipeline

| Step | Folder | What it does |
|---|---|---|
| 1 | `01_calling/` | xTea v0.1.9 command used for Alu and LINE-1 calling (xTea itself: https://github.com/parklab/xTea) |
| 2a | `02_filtering/UPSTREAM_STEPS.md` | Coordinate harmonization (±50 bp), blood subtraction, T–NT comparison, AF ≤ 0.3 |
| 2b | `02_filtering/filter_calls.py` | Recurrent-locus filter (±50 bp, other patients), gnomAD-SV v4.1 same-family filter (±100 bp), genomic-context annotation and impact tiers |
| 3 | `03_detection_power/detection_power.py` | Per-call detection power in the matched tissue and in blood; well-powered and two-sided TPRT subsets |
| 4 | `04_statistics/burden_statistics.py` | Paired Wilcoxon and Spearman tests, outlier sensitivity, age associations with coverage/insert-size adjustment, dominance classes (Table 3), three-compartment Friedman tests, technical covariates, pack-years, patient-level genomic context with GEE (Table 4) |
| 5 | `05_figures/make_figures.py` | Figures 2–5 and Supplementary Figure S5 |

## Running

```bash
pip install -r requirements.txt
# put the input tables in ./data (see data/README.md), then:
./run_all.sh
```

Results are written to `./results` (`statistics.xlsx`, `detection_power_summary.tsv`, `figures/`).
Paths can be changed with the environment variables `MEI_DATA_DIR` and `MEI_RESULTS_DIR`.

## Data availability

No sequencing or patient-level data are included here. See the Data Availability Statement of the
article for how to access the underlying data. The primary call set is available as Supplementary
Table S10 of the article.

## License

MIT (see `LICENSE`). xTea and gnomAD are third-party resources and are subject to their own terms.
