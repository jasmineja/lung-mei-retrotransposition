"""Step 2b: recurrent-locus and gnomAD-SV filters, genomic-context annotation.

Input : data/calls.tsv.gz  (calls after blood subtraction and T-vs-NT comparison, AF <= 0.3)
        data/gnomad_v4.1_INS_ME_subset.tsv.gz  (see make_gnomad_subset.sh)
Output: results/calls_annotated.tsv.gz   (all calls with filter flags)
        results/primary_call_set.tsv.gz  (T- and NT-specific calls passing all filters)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from common import DATA, RESULTS, load_calls  # noqa: E402

RECUR_WIN = 50    # bp; same-family call in any compartment of a different patient
GNOMAD_WIN = 100  # bp; same-family MEI catalogued in gnomAD-SV v4.1


def recurrent_flags(df, win=RECUR_WIN):
    flag = np.zeros(len(df), bool)
    for _, g in df.groupby(["chrom", "element"]):
        g = g.sort_values("pos")
        pos, pid, idx = g.pos.values, g.patient_id.values, g.index.values
        for i in range(len(g)):
            lo = np.searchsorted(pos, pos[i] - win, "left")
            hi = np.searchsorted(pos, pos[i] + win, "right")
            if (pid[lo:hi] != pid[i]).any():
                flag[idx[i]] = True
    return flag


def gnomad_flags(df, gsub, win=GNOMAD_WIN):
    out = pd.Series(False, index=df.index)
    for c, sub in df.groupby("chrom"):
        gs = gsub[gsub.chrom == c]
        if gs.empty:
            continue
        lo, hi = gs.start.values - win, gs.end.values + win
        o = np.argsort(lo)
        lo, hi = lo[o], np.maximum.accumulate(hi[o])
        k = np.searchsorted(lo, sub.pos.values, "right")
        hit = (k > 0) & (hi[np.maximum(k - 1, 0)] >= sub.pos.values)
        out.loc[sub.index] = hit
    return out.values


CONTEXT_PRIORITY = ["exon", "CDS", "stop_codon", "UTR5", "UTR3", "up_stream", "down_stream", "intron", "not_gene_region"]
CONTEXT_MAP = {"exon": "Exon", "CDS": "Exon", "stop_codon": "Exon", "UTR5": "UTR", "UTR3": "UTR",
               "up_stream": "Promoter", "down_stream": "Downstream", "intron": "Intron", "not_gene_region": "Intergenic"}
TIER = {"Nested": "Low", "Intergenic": "Low", "Intron": "Intermediate", "Downstream": "Intermediate",
        "Promoter": "High", "Exon": "High", "UTR": "High"}


def context(ref_rep, gene_info):
    if str(ref_rep).startswith("Fall_in_"):
        return "Nested"
    toks = {e.split(":")[0] for e in str(gene_info).split("/")}
    for p in CONTEXT_PRIORITY:
        if p in toks:
            return CONTEXT_MAP[p]
    return "Intergenic"


def main():
    d = load_calls()
    d["recurrent_50bp"] = recurrent_flags(d)

    g = pd.read_csv(DATA / "gnomad_v4.1_INS_ME_subset.tsv.gz", sep="\t", header=None,
                    names=["chrom", "start", "end", "type", "AF"])
    g["chrom"] = g.chrom.astype(str).str.replace("chr", "").replace({"X": "23", "Y": "24"})
    fam = {"ALU": "INS:ME:ALU", "LINE1": "INS:ME:LINE1"}
    d["gnomad_same_family"] = False
    for el, t in fam.items():
        m = d.element == el
        d.loc[m, "gnomad_same_family"] = gnomad_flags(d[m], g[g.type == t])

    d["context"] = [context(a, b) for a, b in zip(d.REF_REP, d.GENE_INFO)]
    d["tier"] = d.context.map(TIER)
    d["two_sided_tprt"] = d.SUBTYPE.astype(str).str.startswith("two_side")

    d.to_csv(RESULTS / "calls_annotated.tsv.gz", sep="\t", index=False)
    lung = d.compartment.isin(["T", "NT"])
    prim = d[lung & ~d.recurrent_50bp & ~d.gnomad_same_family]
    prim.to_csv(RESULTS / "primary_call_set.tsv.gz", sep="\t", index=False)
    print(f"T/NT calls: {int(lung.sum())}; recurrent: {int((lung & d.recurrent_50bp).sum())}; "
          f"gnomAD-SV same family: {int((lung & d.gnomad_same_family).sum())}; primary call set: {len(prim)}")
    print(prim.groupby(["element", "compartment"]).size().to_string())


if __name__ == "__main__":
    main()
