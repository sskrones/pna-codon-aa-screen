"""Plot the independent OpenMM re-runs (validation/results/*.json) against the original table.

Output: figures/fig7_validation_dissociation.png, figures/fig8_validation_reproducibility.png,
        validation/validation_summary.csv
Usage:  python validation/plot_validation.py
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BLUE, ORANGE, INK, INK2, MUTED, GRID = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "savefig.dpi": 220, "savefig.bbox": "tight", "axes.titleweight": "bold", "axes.titlesize": 11,
    "axes.titlecolor": INK, "axes.titlelocation": "left",
})

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
runs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(here, "results", "*.json")))]
orig = pd.read_csv(os.path.join(root, "results", "mmgbsa_results_md.csv"))

rows = []
for r in runs:
    o = orig[(orig.codon == r["codon"]) & (orig.amino_acid == r["aa"])].iloc[0]
    rows.append({"codon": r["codon"], "amino_acid": r["aa"], "seed": r["seed"],
                 "original_dG_mean": round(o.mmgbsa_dg_bind_mean_kcal_mol, 3),
                 "original_dG_sd": round(o.mmgbsa_dg_bind_std_kcal_mol, 3),
                 "rerun_dG_mean": r["dg_md_mean"], "rerun_dG_sd": r["dg_md_std"],
                 "frac_frames_in_contact": r["frac_frames_contact_lt4A"],
                 "min_dist_last_frame_A": r["mindist_last"], "dG_docked_pose": r["dg_docked_pose"]})
tab = pd.DataFrame(rows).sort_values(["codon", "amino_acid", "seed"])
tab.to_csv(os.path.join(here, "validation_summary.csv"), index=False)
print(tab.to_string(index=False))

# fig7: minimum ligand-receptor heavy-atom distance over production MD (one run per pair)
first = {}
for r in runs:
    first.setdefault((r["codon"], r["aa"]), r)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
t = np.arange(1, 51) * 3
for (c, aa), r in sorted(first.items()):
    ax.plot(t, r["mindist_series"], color=BLUE, alpha=0.55, linewidth=1.4)
ax.axhline(4.0, color=ORANGE, linewidth=1.5, linestyle="--")
ax.text(150, 4.6, "contact (< 4 Å)", color=INK2, fontsize=8.5, ha="right")
ax.set_xlim(0, 150); ax.set_ylim(0, None)
ax.set_xlabel("Production MD time (ps)"); ax.set_ylabel("Min. heavy-atom distance (Å)")
ax.set_title(f"Most ligands drift off the trimer in unrestrained MD ({len(first)} pairs, OpenMM re-run)")
fig.savefig(os.path.join(root, "figures", "fig7_validation_dissociation.png")); plt.close(fig)

# fig8: original table value vs repeated runs for the two top-ranked pairs (+ others)
pairs = tab.groupby(["codon", "amino_acid"]).size().reset_index()[["codon", "amino_acid"]]
pairs["orig"] = [tab[(tab.codon == c) & (tab.amino_acid == a)].original_dG_mean.iloc[0] for c, a in zip(pairs.codon, pairs.amino_acid)]
pairs = pairs.sort_values("orig").reset_index(drop=True)
fig, ax = plt.subplots(figsize=(7.2, 3.6))
for i, p in pairs.iterrows():
    sub = tab[(tab.codon == p.codon) & (tab.amino_acid == p.amino_acid)]
    ax.scatter([i] * len(sub), sub.rerun_dG_mean, s=30, color=BLUE, edgecolor="white", linewidth=1, zorder=3,
               label="OpenMM re-run (one per seed)" if i == 0 else None)
    ax.scatter(i, p.orig, s=70, marker="D", facecolor="none", edgecolor=ORANGE, linewidth=2, zorder=2,
               label="Original table (sander)" if i == 0 else None)
ax.axhline(0, color=INK2, linewidth=1)
ax.set_xticks(range(len(pairs)))
ax.set_xticklabels([f"{c}\n{a}" for c, a in zip(pairs.codon, pairs.amino_acid)], fontsize=8)
ax.set_ylabel("MM-GBSA ΔG (kcal/mol)")
ax.set_title("Top-ranked pairs do not reproduce: their ΔG depends on how long the ligand lingers")
ax.legend(frameon=False, fontsize=8.5, loc="lower right")
fig.savefig(os.path.join(root, "figures", "fig8_validation_reproducibility.png")); plt.close(fig)
