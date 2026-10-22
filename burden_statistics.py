"""Step 4: patient-level statistics reported in the manuscript.

Input : results/primary_call_set.tsv.gz, results/calls_annotated.tsv.gz,
        data/samples.tsv, data/patients.tsv, data/wb_specific_counts.tsv, data/exclusions.tsv
Output: results/statistics.xlsx with one sheet per analysis
"""
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from common import (ELEMENTS, LABEL, RESULTS, counts, load_exclusions, load_samples,  # noqa: E402
                    load_wb_counts, partial_spearman, tukey_mask)

warnings.filterwarnings("ignore")
cov, ins, pat = load_samples()
PIDS = list(pat.index)
EXCL = load_exclusions()


def ids_for(element, compartment=None):
    return [p for p in PIDS if (p, element, compartment or "T") not in EXCL]


def burden(df, label):
    """Totals, medians, paired Wilcoxon, T-NT Spearman, outlier sensitivity, age associations."""
    c = counts(df, PIDS)
    rows = []
    for el in ELEMENTS:
        ids = ids_for(el)
        T, N = c[(el, "T")].reindex(ids), c[(el, "NT")].reindex(ids)
        w, s = stats.wilcoxon(T, N), stats.spearmanr(T, N)
        keep = tukey_mask((T - N).values)
        so = stats.spearmanr(T[keep], N[keep])
        r = dict(call_set=label, element=LABEL[el], n_pairs=len(ids), T_total=int(c[(el, "T")].sum()),
                 NT_total=int(c[(el, "NT")].sum()), T_median_all=c[(el, "T")].median(), NT_median_all=c[(el, "NT")].median(),
                 T_IQR=f'{c[(el, "T")].quantile(.25):.0f}-{c[(el, "T")].quantile(.75):.0f}',
                 NT_IQR=f'{c[(el, "NT")].quantile(.25):.0f}-{c[(el, "NT")].quantile(.75):.0f}',
                 wilcoxon_p=w.pvalue, spearman_rho=s.statistic, spearman_p=s.pvalue,
                 spearman_rho_no_outliers=so.statistic, spearman_p_no_outliers=so.pvalue, n_no_outliers=int(keep.sum()))
        for comp in ["T", "NT"]:
            idc = ids_for(el, comp)
            vec = c[(el, comp)].reindex(idc)
            age = pat.age.reindex(idc).values
            rr = stats.spearmanr(age, vec)
            pr = partial_spearman(age, vec.values, [cov[comp].reindex(idc).values, ins[comp].reindex(idc).values])
            r.update({f"age_{comp}_rho": rr.statistic, f"age_{comp}_p": rr.pvalue,
                      f"age_{comp}_partial_rho": pr[0], f"age_{comp}_partial_p": pr[1]})
        rows.append(r)
    return rows


def dominance(df, thr=0.2):
    """Table 3: (T - NT)/(T + NT) on combined Alu + LINE-1."""
    c = counts(df, PIDS)
    t = sum(c[(e, "T")] for e in ELEMENTS)
    n = sum(c[(e, "NT")] for e in ELEMENTS)
    rel = (t - n) / (t + n).replace(0, np.nan)
    cls = np.where(rel >= thr, "Tumor-dominant", np.where(rel <= -thr, "NT-dominant", "Balanced/mixed"))
    out = pd.Series(cls, index=PIDS).value_counts().rename("n").to_frame()
    out["%"] = 100 * out.n / out.n.sum()
    return out


def three_compartment(df):
    """Friedman + paired Wilcoxon (BH) for T, NT and WB-specific burden; primary set excludes one site's WB."""
    c, wb = counts(df, PIDS), load_wb_counts()
    rows = []
    for scope, excl_site in [("primary (site VU WB excluded)", "VU"), ("all sites", None)]:
        for el in ELEMENTS:
            x = pd.DataFrame({"T": c[(el, "T")], "NT": c[(el, "NT")], "WB": wb[f"WB_{el}"].reindex(PIDS)})
            x = x[x.WB.notna()]
            if excl_site:
                x = x[pat.site.reindex(x.index) != excl_site]
            f = stats.friedmanchisquare(x["T"], x["NT"], x["WB"])
            q = multipletests([stats.wilcoxon(x["T"], x["WB"]).pvalue, stats.wilcoxon(x["NT"], x["WB"]).pvalue,
                               stats.wilcoxon(x["T"], x["NT"]).pvalue], method="fdr_bh")[1]
            rows.append(dict(scope=scope, element=LABEL[el], n=len(x), friedman_p=f.pvalue,
                             q_T_vs_WB=q[0], q_NT_vs_WB=q[1], q_T_vs_NT=q[2]))
    return pd.DataFrame(rows)


def technical(df):
    """Delta(T-NT) burden vs Delta coverage / insert size; WB coverage, insert size and pack-years vs burden."""
    c = counts(df, PIDS)
    rows = []
    for el in ELEMENTS:
        ids = ids_for(el)
        dB = c[(el, "T")].reindex(ids) - c[(el, "NT")].reindex(ids)
        dC = (cov["T"] - cov["NT"]).reindex(ids)
        dI = (ins["T"] - ins["NT"]).reindex(ids)
        for name, x in [("delta coverage", dC), ("delta insert size", dI)]:
            r = stats.spearmanr(x, dB)
            rows.append(dict(test=f"Delta burden vs {name}", element=LABEL[el], n=len(ids), rho=r.statistic, p=r.pvalue))
        pr = partial_spearman(dC.values, dB.values, [dI.values])
        rows.append(dict(test="Delta burden vs delta coverage (partial, adj. insert size)", element=LABEL[el],
                         n=len(ids), rho=pr[0], p=pr[1]))
        for comp in ["T", "NT"]:
            idc = ids_for(el, comp)
            y = c[(el, comp)].reindex(idc)
            r = stats.spearmanr(ins[comp].reindex(idc), y)
            rows.append(dict(test=f"insert size vs burden ({comp})", element=LABEL[el], n=len(idc), rho=r.statistic, p=r.pvalue))
            wbi = [p for p in idc if pd.notna(cov.loc[p].get("WB", np.nan))]
            r = stats.spearmanr(cov.loc[wbi, "WB"], y.reindex(wbi))
            rows.append(dict(test=f"WB coverage vs burden ({comp})", element=LABEL[el], n=len(wbi), rho=r.statistic, p=r.pvalue))
            pyi = [p for p in idc if pd.notna(pat.loc[p, "pack_years"])]
            if pyi:
                r = stats.spearmanr(pat.loc[pyi, "pack_years"].astype(float), y.reindex(pyi))
                rows.append(dict(test=f"pack-years vs burden ({comp})", element=LABEL[el], n=len(pyi), rho=r.statistic, p=r.pvalue))
    out = pd.DataFrame(rows)
    out["q_within_test_family"] = out.groupby(out.test.str.replace(r" \((T|NT)\)", "", regex=True)).p.transform(
        lambda p: multipletests(p, method="fdr_bh")[1])
    return out


def genomic_context(df, min_events=10):
    """Patient-level paired Wilcoxon on category proportions and GEE (binomial, exchangeable, clustered by patient)."""
    cats = ["Nested", "Intergenic", "Intron", "Downstream", "Promoter", "Exon", "UTR"]
    n_ev = df.groupby(["patient_id", "compartment"]).size().unstack()
    elig = n_ev[(n_ev.get("T", 0) >= min_events) & (n_ev.get("NT", 0) >= min_events)].index
    s = df[df.patient_id.isin(elig)]
    prop = s.groupby(["patient_id", "compartment"]).context.value_counts(normalize=True).unstack().fillna(0)
    tier = s.groupby(["patient_id", "compartment"]).tier.value_counts(normalize=True).unstack().fillna(0)
    rows = []
    for cat in cats + ["High tier"]:
        src = tier["High"] if cat == "High tier" else prop.get(cat, pd.Series(0, index=prop.index))
        pt = src.xs("T", level=1).reindex(elig).fillna(0)
        pn = src.xs("NT", level=1).reindex(elig).fillna(0)
        wp = stats.wilcoxon(pt, pn).pvalue if (pt - pn).abs().sum() > 0 else 1.0
        y = (s.tier == "High") if cat == "High tier" else (s.context == cat)
        g = pd.DataFrame({"y": y.astype(int).values, "T": (s.compartment == "T").astype(int).values,
                          "pid": s.patient_id.values}).sort_values("pid")
        m = smf.gee("y ~ T", "pid", g, family=sm.families.Binomial(), cov_struct=sm.cov_struct.Exchangeable()).fit()
        lo, hi = np.exp(m.conf_int().loc["T"])
        rows.append(dict(category=cat, n_patients=len(elig), median_pct_T=100 * pt.median(), median_pct_NT=100 * pn.median(),
                         wilcoxon_p=wp, GEE_OR=np.exp(m.params["T"]), GEE_CI_low=lo, GEE_CI_high=hi, GEE_p=m.pvalues["T"]))
    out = pd.DataFrame(rows)
    out["wilcoxon_q"] = multipletests(out.wilcoxon_p, method="fdr_bh")[1]
    out["GEE_q"] = multipletests(out.GEE_p, method="fdr_bh")[1]
    return out


def main():
    P = pd.read_csv(RESULTS / "primary_call_set.tsv.gz", sep="\t", dtype={"chrom": str})
    sets = {
        "Primary call set": P,
        "Well-powered (>=90%)": P[P.power_matched >= 0.9],
        "Well-powered, conservative": P[P.power_matched_conservative >= 0.9],
        "Two-sided TPRT": P[P.two_sided_tprt],
        "Two-sided TPRT + well-powered": P[P.two_sided_tprt & (P.power_matched >= 0.9)],
    }
    for thr in [0.2, 0.1]:
        sets[f"AF <= {thr}"] = P[P.AF <= thr]
    B = pd.DataFrame(sum([burden(v, k) for k, v in sets.items()], []))
    with pd.ExcelWriter(RESULTS / "statistics.xlsx") as xw:
        B.to_excel(xw, sheet_name="burden_and_sensitivity", index=False)
        dominance(P).to_excel(xw, sheet_name="dominance_Table3")
        three_compartment(P).to_excel(xw, sheet_name="three_compartment", index=False)
        technical(P).to_excel(xw, sheet_name="technical_covariates", index=False)
        genomic_context(P).to_excel(xw, sheet_name="genomic_context_Table4", index=False)
    pd.set_option("display.width", 200)
    print(B[B.call_set == "Primary call set"].T.to_string())
    print(dominance(P).to_string())
    print(three_compartment(P).to_string())


if __name__ == "__main__":
    main()
