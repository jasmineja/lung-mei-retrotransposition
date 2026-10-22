"""Step 5: Figures 2-5 and Supplementary Figure S5 of the manuscript.

Input : results/primary_call_set.tsv.gz, data/patients.tsv, data/wb_specific_counts.tsv, data/exclusions.tsv
Output: results/figures/*.png, *.tiff (300 dpi) and *.pdf
Requires the Liberation Sans font (or edit FONT below), and adjustText for label placement.
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from adjustText import adjust_text  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from common import RESULTS, load_exclusions, load_samples, load_wb_counts  # noqa: E402

FONT = "Liberation Sans"
plt.rcParams.update({"font.family": FONT, "font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 1.0, "savefig.dpi": 300})
OUT = str(RESULTS / "figures") + "/"
Path(OUT).mkdir(parents=True, exist_ok=True)

P = pd.read_csv(RESULTS / "primary_call_set.tsv.gz", sep="\t", dtype={"chrom": str})
_, _, pat = load_samples()
info = pat.rename(columns={"site": "loc", "histology": "subtype"})
PIDS = list(info.index)
c = (P.groupby(["patient_id", "element", "compartment"]).size()
       .unstack(["element", "compartment"]).reindex(PIDS).fillna(0))
wb = load_wb_counts()
tb = pd.DataFrame({"Somatic_AF<=0.3_Count_WB_ALU": wb.WB_ALU, "Somatic_AF<=0.3_Count_WB_L1": wb.WB_LINE1})
EXCL = load_exclusions()
L1_OK = [p for p in PIDS if (p, "LINE1", "T") not in EXCL]
def hist(s):
    return 'Adenocarcinoma' if s == 'Adenocarcinoma' else ('Squamous cell carcinoma' if s == 'Squamous cell carcinoma' else 'Other histologies')
info['grp'] = info.subtype.map(hist)
def fmt_p(p):
    if p < 1e-3:
        m, e = f"{p:.1e}".split('e'); return f"{m} × 10$^{{{int(e)}}}$"
    return f"{p:.2g}" if p < 0.01 else f"{p:.2f}" if p >= 0.1 else f"{p:.3f}".rstrip('0')
stats_out = {}

# ------------------------------------------------------------------ Figure 2
keep = [p for p in PIDS if info.loc[p, 'loc'] != 'VU' and p in tb.index and pd.notna(tb.loc[p, 'Somatic_AF<=0.3_Count_WB_ALU'])]
fig, axes = plt.subplots(1, 2, figsize=(11, 5.6))
for ax, (mei, lab, wbcol, letter) in zip(axes, [('LINE1', 'LINE-1', 'Somatic_AF<=0.3_Count_WB_L1', 'A'), ('ALU', 'Alu', 'Somatic_AF<=0.3_Count_WB_ALU', 'B')]):
    df = pd.DataFrame({'T': c.loc[keep, (mei, 'T')], 'NT': c.loc[keep, (mei, 'NT')], 'WB': tb.loc[keep, wbcol]})
    f = stats.friedmanchisquare(df['T'], df['NT'], df['WB']).pvalue
    q = multipletests([stats.wilcoxon(df['T'], df['WB']).pvalue, stats.wilcoxon(df['NT'], df['WB']).pvalue, stats.wilcoxon(df['T'], df['NT']).pvalue], method='fdr_bh')[1]
    stats_out['fig2_' + mei] = dict(n=len(df), friedman=f, q=list(q))
    Y = np.log10(1 + df.values)
    for row in Y:
        ax.plot([0, 1, 2], row, color='#c8c8c8', lw=0.7, zorder=1)
    ax.boxplot([Y[:, 0], Y[:, 1], Y[:, 2]], positions=[0, 1, 2], widths=0.45, showfliers=False,
               medianprops=dict(color='black', lw=2), boxprops=dict(lw=1.2), whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2), zorder=2)
    for k in range(3):
        ax.scatter(np.full(len(Y), k), Y[:, k], s=28, color='black', zorder=3)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['Tumor\n(T)', 'Adjacent normal\n(NT)', 'Whole blood\n(WB)'])
    ax.set_ylabel(f'Compartment-specific {"Alu" if mei=="ALU" else "LINE-1"} insertions\nlog$_{{10}}$(1 + count)')
    top = Y.max() * 1.08 + 0.1
    ax.set_ylim(0, top * 1.42)
    ax.set_xlim(-0.5, 2.5)
    yb = top * 1.03
    hh = top * 0.022
    for (a, b, txt, qq) in [(0, 0.95, 'T vs WB', q[0]), (1.05, 2, 'NT vs WB', q[1])]:
        pass
    # pairwise brackets
    for (x0, x1, txt, qq, h) in [(-0.2, 0.2, '', None, 0)]:
        pass
    pairs = [((0, 2), 'T vs WB', q[0], 0), ((1, 2), 'NT vs WB', q[1], 1), ((0, 1), 'T vs NT', q[2], 2)]
    xs = {0: (-0.25, 0.75), 1: (0.85, 1.65), 2: (1.75, 2.45)}
    for (pr, txt, qq, slot) in pairs:
        x0, x1 = xs[slot]
        ax.plot([x0, x0, x1, x1], [yb, yb + hh, yb + hh, yb], color='black', lw=1)
        ax.text((x0 + x1) / 2, yb + 2*hh, f'{txt}\n$q$ = {fmt_p(qq)}', ha='center', va='bottom', fontsize=10)
    yf = top * 1.22
    ax.plot([-0.25, -0.25, 2.45, 2.45], [yf, yf + hh, yf + hh, yf], color='black', lw=1)
    ax.text(1.1, yf + 2*hh, f'n = {len(df)} patients (complete T/NT/WB)\nFriedman test: P = {fmt_p(f)}', ha='center', va='bottom', fontsize=10.5)
    ax.set_yticks([t for t in ax.get_yticks() if t <= top])
    ax.spines['left'].set_bounds(0, top)
    ax.set_title(f'{letter}   {"Alu" if mei=="ALU" else lab}', loc='left', fontweight='bold', fontsize=16, pad=12)
fig.tight_layout()
fig.savefig(OUT + 'Figure2.png'); fig.savefig(OUT + 'Figure2.tiff', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + 'Figure2.pdf'); plt.close(fig)

# ------------------------------------------------------------------ Figure 3
fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.6))
for ax, (mei, lab, letter, ids) in zip(axes, [('LINE1', 'LINE-1', 'A', L1_OK), ('ALU', 'Alu', 'B', PIDS)]):
    T = c.loc[ids, (mei, 'T')].values; N = c.loc[ids, (mei, 'NT')].values
    w = stats.wilcoxon(T, N); dlt = T - N; dlt = dlt[dlt != 0]
    rk = stats.rankdata(np.abs(dlt)); rrb = abs(rk[dlt > 0].sum() - rk[dlt < 0].sum()) / rk.sum()
    stats_out['fig3_' + mei] = dict(n=len(ids), p=w.pvalue, rrb=rrb, detT=int((T > 0).sum()), detN=int((N > 0).sum()))
    yT, yN = np.log10(1 + T), np.log10(1 + N)
    for a, b in zip(yT, yN):
        ax.plot([0, 1], [a, b], color='#d0d0d0', lw=0.8, zorder=1)
    ax.boxplot([yT, yN], positions=[0, 1], widths=0.38, whis=(0, 100), medianprops=dict(color='black', lw=2), boxprops=dict(lw=1.4), whiskerprops=dict(lw=1.2), capprops=dict(lw=1.2), zorder=2)
    ax.scatter(np.zeros(len(yT)), yT, s=40, color='black', zorder=3)
    ax.scatter(np.ones(len(yN)), yN, s=40, facecolor='white', edgecolor='black', lw=1.1, zorder=3)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['Tumor', 'Adjacent normal'], fontsize=12)
    ax.set_xlim(-0.45, 1.45)
    ax.yaxis.grid(True, color='#e6e6e6', lw=0.8); ax.set_axisbelow(True)
    ax.set_ylabel('Tissue-restricted insertions, log$_{10}$(1 + count)', fontsize=12)
    ax.set_title(f'{letter}   {"Alu" if mei=="ALU" else lab}', loc='left', fontweight='bold', fontsize=15, pad=46)
    ax.text(0.5, 1.10, f'n = {len(ids)} paired patients   |   Wilcoxon P = {w.pvalue:.3f}   |   $r_{{rb}}$ = {rrb:.2f}', transform=ax.transAxes, ha='center', fontsize=11)
    ax.text(0.5, 1.04, f'Detection: Tumor {int((T>0).sum())}/{len(ids)}   |   Adjacent normal {int((N>0).sum())}/{len(ids)}', transform=ax.transAxes, ha='center', fontsize=10, color='#444444')
fig.tight_layout()
fig.savefig(OUT + 'Figure3.png'); fig.savefig(OUT + 'Figure3.tiff', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + 'Figure3.pdf'); plt.close(fig)

# ------------------------------------------------------------------ Figure 4 (two images)
for mei, lab, letter, ids in [('LINE1', 'LINE-1', 'A', L1_OK), ('ALU', 'Alu', 'B', PIDS)]:
    dl = (c.loc[ids, (mei, 'T')] - c.loc[ids, (mei, 'NT')]).sort_values()
    fig, ax = plt.subplots(figsize=(10, 5.4))
    x = np.arange(len(dl))
    cols = ['#c0c0c0' if v < 0 else 'black' for v in dl.values]
    edges = ['black' if (i < 5 and dl.values[i] < 0) else 'none' for i in range(len(dl))]
    lws = [1.8 if e == 'black' else 0 for e in edges]
    ax.bar(x, dl.values, color=cols, edgecolor=edges, linewidth=lws, width=0.85)
    for i in range(len(dl)):
        if dl.values[i] < 0 and edges[i] == 'none':
            ax.bar(x[i], dl.values[i], color='#c0c0c0', edgecolor='#555555', linewidth=0.5, width=0.85)
    ax.axhline(0, color='black', lw=1.2)
    rng = dl.abs().max()
    for i in list(range(5)) + list(range(len(dl) - 5, len(dl))):
        v = dl.values[i]
        if v == 0: continue
        off = rng * 0.02
        ax.text(x[i], v - off if v < 0 else v + off, dl.index[i], rotation=90, ha='center', va='top' if v < 0 else 'bottom', fontsize=8.5)
    ax.set_xticks([]); ax.spines['bottom'].set_visible(False)
    ax.set_xlabel('Patients ranked by Δ(T−NT)', fontsize=12)
    ax.set_ylabel(f'Δ(Tumor − NT) tissue-restricted {"Alu" if mei=="ALU" else "LINE-1"} insertions', fontsize=12)
    ax.set_ylim(-rng * 1.30, rng * 1.42)
    ax.text(0.03, 0.95, 'Adjacent normal enriched', transform=ax.transAxes, fontsize=11, color='#555555', va='top')
    ax.text(0.86, 0.95, 'Tumor enriched', transform=ax.transAxes, fontsize=11, ha='right', va='top')
    ax.text(0.03, 0.86, f'n = {len(dl)}', transform=ax.transAxes, fontsize=11, va='top')
    ax.set_title(f'{letter}   {"Alu" if mei=="ALU" else lab}', loc='left', fontweight='bold', fontsize=15)
    fig.tight_layout()
    nm = f'Figure4{letter}'
    fig.savefig(OUT + nm + '.png'); fig.savefig(OUT + nm + '.tiff', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + nm + '.pdf'); plt.close(fig)

# ------------------------------------------------------------------ Figure 5 (two images)
GRP = ['Adenocarcinoma', 'Squamous cell carcinoma', 'Other histologies']
STY = {'Adenocarcinoma': dict(marker='o', facecolor='white', edgecolor='black'),
       'Squamous cell carcinoma': dict(marker='s', facecolor='#808080', edgecolor='black'),
       'Other histologies': dict(marker='^', facecolor='black', edgecolor='black')}
def fence_mask(v):
    q1, q3 = np.percentile(v, [25, 75]); return (v >= q1 - 3 * (q3 - q1)) & (v <= q3 + 3 * (q3 - q1))
rng = np.random.default_rng(7)
jit = {p: rng.uniform(-0.33, 0.33) for p in PIDS}
nout = {}
for variant, letter in [('all', 'A'), ('excl', 'B')]:
    fig, axes = plt.subplots(1, 2, figsize=(14, 7.2))
    for ax, (mei, ylabel, ids) in zip(axes, [('ALU', 'Alu Δ(T − NT) (tumor-specific − NT-specific insertions)', PIDS), ('LINE1', 'L1 Δ(T − NT) (tumor-specific − NT-specific insertions)', L1_OK)]):
        dl = (c.loc[ids, (mei, 'T')] - c.loc[ids, (mei, 'NT')])
        ok = pd.Series(fence_mask(dl.values), index=dl.index)
        nout[mei] = int((~ok).sum())
        texts = []
        for gi, g in enumerate(GRP):
            gid = [p for p in ids if info.loc[p, 'grp'] == g]
            show = gid if variant == 'all' else [p for p in gid if ok[p]]
            xs = np.array([gi + jit[p] for p in show]); ys = dl[show].values
            ax.scatter(xs, ys, s=55, linewidths=1.0, zorder=3, **STY[g])
            med = np.median(dl[show].values) if len(show) else 0
            ax.plot([gi - 0.35, gi + 0.35], [med, med], color='black', lw=2.2, zorder=2)
            if variant == 'all':
                R_ = dl.abs().max()
                for p, xx in zip(show, xs):
                    if not ok[p]:
                        v = dl[p]; yy = v
                        if abs(v) < 0.2 * R_:
                            yy = v + (0.09 if v >= 0 else -0.09) * R_
                        texts.append(ax.text(xx + 0.12, yy, p, fontsize=8.5))
        ax.set_xticks(range(3)); ax.set_xticklabels([f'{g} (n={sum(info.loc[ids,"grp"]==g)})' for g in GRP], rotation=20, ha='right')
        ax.set_xlim(-0.6, 2.6); ax.yaxis.grid(True, color='#e0e0e0'); ax.set_axisbelow(True)
        ax.set_ylabel(ylabel, fontsize=11)
        if texts:
            adjust_text(texts, ax=ax, x=list(np.concatenate([o.get_offsets()[:,0] for o in ax.collections])), y=list(np.concatenate([o.get_offsets()[:,1] for o in ax.collections])), arrowprops=dict(arrowstyle='-', color='#555555', lw=0.7), expand=(1.6, 2.2), force_text=(0.6, 0.9))
    sub = 'all samples; extreme outliers labeled, Tukey 3×IQR' if variant == 'all' else 'extreme outliers excluded, Tukey 3×IQR'
    fig.suptitle(f'Tumor-minus-NT difference in retrotransposition burden by histologic subtype\n({sub})', fontsize=14)
    handles = [plt.Line2D([], [], linestyle='none', markersize=9, markeredgewidth=1.0, marker=STY[g]['marker'], markerfacecolor=STY[g]['facecolor'], markeredgecolor='black', label=f'{g} (n={sum(info.grp==g)})') for g in GRP]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0.06, 1, 0.95))
    nm = f'Figure5{letter}'
    fig.savefig(OUT + nm + '.png'); fig.savefig(OUT + nm + '.tiff', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + nm + '.pdf'); plt.close(fig)
stats_out['fig5_outliers'] = nout

# ------------------------------------------------------------------ Supplementary Figure S5
sm = pat[['pack_years']]
fig, axes = plt.subplots(2, 2, figsize=(10, 8.6))
for ax, (mei, comp, title, letter) in zip(axes.flat, [('ALU', 'T', 'Alu — Tumor', 'A'), ('ALU', 'NT', 'Alu — Adjacent normal', 'B'), ('LINE1', 'T', 'LINE-1 — Tumor', 'C'), ('LINE1', 'NT', 'LINE-1 — Adjacent normal', 'D')]):
    ids = [p for p in PIDS if p in sm.index and pd.notna(sm.loc[p, 'pack_years']) and not (mei == 'LINE1' and comp == 'T' and p not in L1_OK)]
    x = sm.loc[ids, 'pack_years'].astype(float).values; yraw = c.loc[ids, (mei, comp)].values; y = np.log10(1 + yraw)
    r = stats.spearmanr(x, yraw)
    stats_out[f'S5_{mei}_{comp}'] = dict(n=len(ids), rho=r.statistic, p=r.pvalue)
    ax.scatter(x, y, s=45, color='#8c8c8c', edgecolor='black', lw=0.7, zorder=3)
    b = np.polyfit(x, y, 1); xx = np.linspace(x.min(), x.max(), 50); ax.plot(xx, np.polyval(b, xx), color='black', lw=1.3)
    ax.grid(True, color='#e6e6e6', lw=0.8); ax.set_axisbelow(True)
    ax.set_xlabel('Pack-years', fontsize=12); ax.set_ylabel('log$_{10}$(1 + tissue-restricted insertion count)', fontsize=11)
    ax.set_title(f'{title}\n', fontweight='bold', fontsize=12)
    ax.text(0.5, 1.03, f'n = {len(ids)}   Spearman $\\rho$ = {f'{r.statistic:.2f}'.replace('-','−')}   $P$ = {r.pvalue:.3f}', transform=ax.transAxes, ha='center', fontsize=10)
    ax.text(-0.12, 1.12, letter, transform=ax.transAxes, fontweight='bold', fontsize=15)
fig.tight_layout(h_pad=2.5)
fig.savefig(OUT + 'FigureS5.png'); fig.savefig(OUT + 'FigureS5.tiff', pil_kwargs={'compression': 'tiff_lzw'}); fig.savefig(OUT + 'FigureS5.pdf'); plt.close(fig)
json.dump(stats_out, open(OUT + 'fig_stats.json', 'w'), indent=1, default=float)
print(json.dumps(stats_out, indent=1, default=float))
