"""Step 3: per-call detection power in the matched tissue and in blood.

For each T- or NT-specific call, expected read support in the matched sample is the
observed support (split + discordant reads) scaled by the coverage ratio; detection
probability is P(Poisson(lambda) >= k_min), where k_min is the minimum support among
PASS calls. A conservative variant deflates lambda by 25%.
WB power is computed for a constitutional heterozygous insertion (expected AF 0.5).

Input : results/calls_annotated.tsv.gz, data/samples.tsv, data/patients.tsv
Output: results/calls_annotated.tsv.gz (columns added), results/detection_power_summary.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from common import RESULTS, load_samples  # noqa: E402

OTHER = {"T": "NT", "NT": "T"}


def p_detect(lam, k):
    return 1 - stats.poisson.cdf(k - 1, lam)


def main():
    d = pd.read_csv(RESULTS / "calls_annotated.tsv.gz", sep="\t", dtype={"chrom": str})
    cov, _, _ = load_samples()
    k_min = int(d.support.min())
    lung = d.compartment.isin(["T", "NT"])
    L = d[lung]
    c_self = np.array([cov.loc[p, s] for p, s in zip(L.patient_id, L.compartment)], float)
    c_oth = np.array([cov.loc[p, OTHER[s]] for p, s in zip(L.patient_id, L.compartment)], float)
    c_wb = np.array([cov.loc[p, "WB"] if "WB" in cov.columns else np.nan for p in L.patient_id], float)
    lam = L.support.values * c_oth / c_self
    d.loc[lung, "power_matched"] = p_detect(lam, k_min)
    d.loc[lung, "power_matched_conservative"] = p_detect(lam * 0.75, k_min)
    af = L.AF.clip(lower=0.02).values
    lam_wb = L.support.values * c_wb / c_self * np.minimum(0.5 / af, 10)
    d.loc[lung, "power_wb_heterozygous"] = p_detect(lam_wb, k_min)
    d["well_powered"] = d.power_matched >= 0.9
    d.to_csv(RESULTS / "calls_annotated.tsv.gz", sep="\t", index=False)

    P = d[lung & ~d.recurrent_50bp & ~d.gnomad_same_family]
    P.to_csv(RESULTS / "primary_call_set.tsv.gz", sep="\t", index=False)
    P = d[lung]  # summary over all T/NT-specific calls, as reported in the manuscript
    wb = P.power_wb_heterozygous.dropna()
    summ = pd.Series({
        "k_min (minimum read support)": k_min,
        "median power in matched tissue": P.power_matched.median(),
        "% calls with >=90% power": 100 * (P.power_matched >= 0.9).mean(),
        "% calls with >=90% power (conservative)": 100 * (P.power_matched_conservative >= 0.9).mean(),
        "% calls with >=99% WB power (heterozygous)": 100 * (wb >= 0.99).mean(),
        "% calls with two-sided TPRT": 100 * P.two_sided_tprt.mean(),
    })
    summ.to_csv(RESULTS / "detection_power_summary.tsv", sep="\t", header=["value"])
    print(summ.to_string())


if __name__ == "__main__":
    main()
