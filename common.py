"""Shared input loading and helper functions.

All paths are relative to the data directory, which defaults to ./data and can be
changed with the environment variable MEI_DATA_DIR. Outputs go to ./results
(or MEI_RESULTS_DIR). No patient data are stored in this repository.
"""
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("MEI_DATA_DIR", ROOT / "data"))
RESULTS = Path(os.environ.get("MEI_RESULTS_DIR", ROOT / "results"))
RESULTS.mkdir(parents=True, exist_ok=True)

ELEMENTS = ["ALU", "LINE1"]
LABEL = {"ALU": "Alu", "LINE1": "LINE-1"}

# Patients excluded from specific analyses (see manuscript Methods); edit data/exclusions.tsv
def load_exclusions():
    f = DATA / "exclusions.tsv"
    if not f.exists():
        return set()
    x = pd.read_csv(f, sep="\t")
    return {(r.patient_id, r.element, r.compartment) for r in x.itertuples()}


def load_calls():
    """Event-level calls after blood subtraction and T-vs-NT comparison (AF <= 0.3).

    Required columns: patient_id, compartment (T, NT or WB), element (ALU/LINE1),
    chrom, pos, AF, LCLIP, RCLIP, LDISC, RDISC, SUBTYPE, REF_REP, GENE_INFO.
    """
    d = pd.read_csv(DATA / "calls.tsv.gz", sep="\t", dtype={"chrom": str})
    d = d[d.element.isin(ELEMENTS)].copy()
    d = d.drop_duplicates(["patient_id", "compartment", "element", "chrom", "pos"]).reset_index(drop=True)
    d["support"] = d[["LCLIP", "RCLIP", "LDISC", "RDISC"]].sum(axis=1)
    return d


def load_samples():
    """Per-sample metrics and per-patient clinical variables.

    samples.tsv: patient_id, compartment, coverage, mean_insert_size
    patients.tsv: patient_id, site, age, histology, pack_years (may be empty)
    """
    s = pd.read_csv(DATA / "samples.tsv", sep="\t")
    p = pd.read_csv(DATA / "patients.tsv", sep="\t").set_index("patient_id").sort_index()
    cov = s.pivot_table(index="patient_id", columns="compartment", values="coverage", aggfunc="first")
    ins = s.pivot_table(index="patient_id", columns="compartment", values="mean_insert_size", aggfunc="first")
    return cov.reindex(p.index), ins.reindex(p.index), p


def load_wb_counts():
    """Patient-level WB-specific AF <= 0.3 counts: patient_id, WB_ALU, WB_LINE1."""
    return pd.read_csv(DATA / "wb_specific_counts.tsv", sep="\t").set_index("patient_id")


def counts(df, pids):
    """Patient x (element, compartment) count matrix."""
    return (df.groupby(["patient_id", "element", "compartment"]).size()
              .unstack(["element", "compartment"]).reindex(pids).fillna(0))


def partial_spearman(x, y, covariates):
    """Spearman correlation of x and y after rank-residualizing on covariates."""
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    Z = np.column_stack([stats.rankdata(z) for z in covariates] + [np.ones(len(x))])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    r = np.corrcoef(ex, ey)[0, 1]
    n, k = len(x), len(covariates)
    t = r * np.sqrt((n - 2 - k) / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), n - 2 - k)


def tukey_mask(v, k=3.0):
    q1, q3 = np.percentile(v, [25, 75])
    return (v >= q1 - k * (q3 - q1)) & (v <= q3 + k * (q3 - q1))
