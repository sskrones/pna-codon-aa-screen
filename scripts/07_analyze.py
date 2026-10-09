"""Step 7 - statistics and figures from results/mmgbsa_results_md.csv.

Outputs
  results/tables/vina_matrix.csv, mmgbsa_matrix.csv          64 x 15 matrices
  results/tables/amino_acid_summary.csv, codon_summary.csv
  results/tables/variance_decomposition.csv                  two-way additive decomposition
  results/tables/chirality_tests.csv                         paired Wilcoxon, L vs D
  results/tables/summary_stats.json                          numbers quoted in the report
  figures/fig1..fig6 (*.png)

Usage: python scripts/07_analyze.py [--results results/mmgbsa_results_md.csv] [--ligands data/ligands]
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import stats

# ---- palette (validated reference palette, light mode) ----
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1"
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.dpi": 100, "savefig.dpi": 220, "savefig.bbox": "tight", "axes.titleweight": "bold",
    "axes.titlesize": 11, "axes.titlecolor": INK, "axes.titlelocation": "left",
})

AA_NAMES = {
    "AIB": "Aib", "ALA": "L-Ala", "ASP": "L-Asp", "BAL": "β-Ala", "DAL": "D-Ala", "DAS": "D-Asp",
    "DGN": "D-Gln", "DGU": "D-Glu", "DVA": "D-Val", "GAB": "GABA", "GLU": "L-Glu", "GLY": "Gly",
    "ISV": "L-Iva", "MET": "L-Met", "VAL": "L-Val",
}
DL_PAIRS = {"ALA": "DAL", "ASP": "DAS", "GLU": "DGU", "VAL": "DVA"}


def heavy_atoms(lig_dir, code):
    return sum(1 for l in open(os.path.join(lig_dir, f"{code.lower()}.pdbqt"))
               if l.startswith(("ATOM", "HETATM")) and l[77:79].strip() not in ("H", "HD"))


def decompose(m):
    """Share of total variance explained by row (codon), column (amino acid) and residual."""
    g = m.values.mean()
    tot = ((m.values - g) ** 2).sum()
    aa = m.shape[0] * ((m.mean(axis=0) - g) ** 2).sum()
    co = m.shape[1] * ((m.mean(axis=1) - g) ** 2).sum()
    resid = m.sub(m.mean(axis=1), axis=0).sub(m.mean(axis=0), axis=1) + g
    return {"amino_acid": aa / tot, "codon": co / tot, "interaction_residual": 1 - (aa + co) / tot,
            "residual_sd_kcal_mol": float(resid.values.std())}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="results/mmgbsa_results_md.csv")
    ap.add_argument("--ligands", default="data/ligands")
    ap.add_argument("--tables", default="results/tables")
    ap.add_argument("--figures", default="figures")
    a = ap.parse_args()
    os.makedirs(a.tables, exist_ok=True); os.makedirs(a.figures, exist_ok=True)

    d = pd.read_csv(a.results)
    d = d[d.status == "OK"].copy()
    d["heavy_atoms"] = d.amino_acid.map({c: heavy_atoms(a.ligands, c) for c in d.amino_acid.unique()})
    d["n_DAP"] = d.codon.str.count("D")
    V = d.pivot(index="codon", columns="amino_acid", values="vina_affinity_kcal_mol")
    G = d.pivot(index="codon", columns="amino_acid", values="mmgbsa_dg_bind_mean_kcal_mol")
    V.to_csv(f"{a.tables}/vina_matrix.csv"); G.to_csv(f"{a.tables}/mmgbsa_matrix.csv")

    aa_sum = d.groupby("amino_acid").agg(
        vina_mean=("vina_affinity_kcal_mol", "mean"), vina_sd=("vina_affinity_kcal_mol", "std"),
        mmgbsa_mean=("mmgbsa_dg_bind_mean_kcal_mol", "mean"), mmgbsa_sd=("mmgbsa_dg_bind_mean_kcal_mol", "std"),
        heavy_atoms=("heavy_atoms", "first")).sort_values("vina_mean")
    aa_sum.round(4).to_csv(f"{a.tables}/amino_acid_summary.csv")
    co_sum = d.groupby("codon").agg(
        vina_mean=("vina_affinity_kcal_mol", "mean"), mmgbsa_mean=("mmgbsa_dg_bind_mean_kcal_mol", "mean"),
        n_DAP=("n_DAP", "first")).sort_values("vina_mean")
    co_sum.round(4).to_csv(f"{a.tables}/codon_summary.csv")

    dv, dg = decompose(V), decompose(G)
    pd.DataFrame({"vina": dv, "mmgbsa": dg}).round(4).to_csv(f"{a.tables}/variance_decomposition.csv")

    chir = []
    for l, dd in DL_PAIRS.items():
        for name, M in (("vina", V), ("mmgbsa", G)):
            diff = M[dd] - M[l]
            p = stats.wilcoxon(M[l], M[dd]).pvalue
            chir.append({"pair": f"{l}/{dd}", "metric": name, "median_D_minus_L": diff.median(),
                         "mean_D_minus_L": diff.mean(), "wilcoxon_p": p, "bonferroni_p": min(1.0, p * 8)})
    chir = pd.DataFrame(chir); chir.round(5).to_csv(f"{a.tables}/chirality_tests.csv", index=False)

    sem = d.mmgbsa_dg_bind_std_kcal_mol / np.sqrt(d.n_snapshots)
    rho_size = stats.spearmanr(aa_sum.vina_mean, aa_sum.heavy_atoms)
    rho_dap = stats.spearmanr(co_sum.n_DAP, co_sum.vina_mean)
    rho_vg = stats.spearmanr(d.vina_affinity_kcal_mol, d.mmgbsa_dg_bind_mean_kcal_mol)
    rho_ms = stats.spearmanr(d.mmgbsa_dg_bind_mean_kcal_mol, d.mmgbsa_dg_bind_std_kcal_mol)
    summ = {
        "n_pairs": int(len(d)), "n_codons": int(d.codon.nunique()), "n_amino_acids": int(d.amino_acid.nunique()),
        "vina_min": d.vina_affinity_kcal_mol.min(), "vina_max": d.vina_affinity_kcal_mol.max(),
        "vina_mean": d.vina_affinity_kcal_mol.mean(),
        "mmgbsa_min": d.mmgbsa_dg_bind_mean_kcal_mol.min(), "mmgbsa_max": d.mmgbsa_dg_bind_mean_kcal_mol.max(),
        "mmgbsa_frac_positive": float((d.mmgbsa_dg_bind_mean_kcal_mol > 0).mean()),
        "mmgbsa_frac_abs_mean_lt_2sd": float((d.mmgbsa_dg_bind_mean_kcal_mol.abs() < 2 * d.mmgbsa_dg_bind_std_kcal_mol).mean()),
        "mmgbsa_frac_abs_mean_gt_2sem": float((d.mmgbsa_dg_bind_mean_kcal_mol.abs() > 2 * sem).mean()),
        "mmgbsa_n_below_minus1": int((d.mmgbsa_dg_bind_mean_kcal_mol < -1).sum()),
        "variance_vina": dv, "variance_mmgbsa": dg,
        "spearman_aa_vina_vs_heavy_atoms": [rho_size.statistic, rho_size.pvalue],
        "spearman_codon_vina_vs_nDAP": [rho_dap.statistic, rho_dap.pvalue],
        "spearman_vina_vs_mmgbsa_all_pairs": [rho_vg.statistic, rho_vg.pvalue],
        "spearman_mmgbsa_mean_vs_sd": [rho_ms.statistic, rho_ms.pvalue],
        "top_codons_vina": co_sum.head(5).vina_mean.round(3).to_dict(),
        "bottom_codons_vina": co_sum.tail(5).vina_mean.round(3).to_dict(),
    }
    json.dump(summ, open(f"{a.tables}/summary_stats.json", "w"), indent=2, default=float)

    # ---------------- figures ----------------
    cmap = LinearSegmentedColormap.from_list("seq", SEQ[::-1])  # darker = more negative = stronger
    aa_order = aa_sum.index.tolist()
    co_order = co_sum.index.tolist()

    # fig1 heatmap (amino acids x codons)
    fig, ax = plt.subplots(figsize=(13, 4.6))
    M = V.loc[co_order, aa_order].T
    im = ax.imshow(M.values, aspect="auto", cmap=cmap)
    ax.set_xticks(range(len(co_order))); ax.set_xticklabels(co_order, rotation=90, fontsize=7, family="DejaVu Sans Mono")
    ax.set_yticks(range(len(aa_order))); ax.set_yticklabels([AA_NAMES[x] for x in aa_order], fontsize=9)
    ax.grid(False)
    for s in ax.spines.values(): s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("Vina affinity (kcal/mol)", color=INK2)
    cb.outline.set_visible(False)
    ax.set_title("Vina affinity, 64 codons × 15 amino acids (both axes sorted by mean)")
    ax.set_xlabel("Codon (strongest → weakest mean)"); ax.set_ylabel("")
    fig.savefig(f"{a.figures}/fig1_vina_heatmap.png"); plt.close(fig)

    # fig2 variance decomposition
    fig, ax = plt.subplots(figsize=(6.4, 2.6))
    cats = ["amino_acid", "codon", "interaction_residual"]
    labels = ["Amino acid", "Codon", "Codon × amino acid\n(specificity + noise)"]
    cols = [BLUE, AQUA, ORANGE]
    for i, (name, dd_) in enumerate((("Vina", dv), ("MM-GBSA", dg))):
        left = 0
        for c, lab, col in zip(cats, labels, cols):
            w = dd_[c]
            ax.barh(i, w, left=left, color=col, height=0.55, edgecolor="white", linewidth=2,
                    label=lab if i == 0 else None)
            if w > 0.06:
                ax.text(left + w / 2, i, f"{w:.0%}", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
            left += w
    ax.set_yticks([0, 1]); ax.set_yticklabels(["Vina", "MM-GBSA"]); ax.invert_yaxis()
    ax.set_xlim(0, 1); ax.xaxis.set_major_formatter(matplotlib.ticker.PercentFormatter(1.0)); ax.grid(False)
    ax.legend(ncol=3, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.18), fontsize=8.5)
    ax.set_title("Where the variance comes from (two-way additive decomposition)")
    fig.savefig(f"{a.figures}/fig2_variance_decomposition.png"); plt.close(fig)

    # fig3 amino-acid size effect
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    ax.scatter(aa_sum.heavy_atoms, aa_sum.vina_mean, s=46, color=BLUE, edgecolor="white", linewidth=1.5, zorder=3)
    # hand-placed label offsets (points) so that close/identical points stay readable
    offs = {"GLY": (8, 0), "BAL": (8, 0), "ALA": (8, -7), "DAL": (8, 5), "GAB": (8, 2), "AIB": (8, 0),
            "VAL": (8, -9), "DVA": (8, 0), "ISV": (8, 9), "MET": (8, 0), "ASP": (8, -7), "DAS": (8, 5),
            "DGU": (8, 4), "GLU": (8, -2), "DGN": (8, -6)}
    for code, r in aa_sum.iterrows():
        dx, dy = offs.get(code, (8, 0))
        ax.annotate(AA_NAMES[code], (r.heavy_atoms, r.vina_mean), xytext=(dx, dy),
                    textcoords="offset points", fontsize=7.5, color=INK2, va="center")
    ax.set_xlabel("Heavy atoms in amino acid"); ax.set_ylabel("Mean Vina affinity over 64 codons (kcal/mol)")
    ax.set_title(f"Docking score tracks ligand size (Spearman ρ = {rho_size.statistic:.2f})")
    ax.set_xlim(4.3, 11.8)
    fig.savefig(f"{a.figures}/fig3_aa_size_effect.png"); plt.close(fig)

    # fig4 codon DAP effect
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    rng = np.random.default_rng(0)
    for n in range(4):
        y = co_sum[co_sum.n_DAP == n].vina_mean
        ax.scatter(n + rng.uniform(-0.12, 0.12, len(y)), y, s=22, color=BLUE, alpha=0.75, edgecolor="white", linewidth=0.8, zorder=3)
        ax.plot([n - 0.25, n + 0.25], [y.mean()] * 2, color=INK, linewidth=2, zorder=4)
    for c in co_sum.head(3).index:
        r = co_sum.loc[c]
        ax.annotate(c, (r.n_DAP, r.vina_mean), xytext=(8, -10), textcoords="offset points", fontsize=7.5,
                    color=INK2, family="DejaVu Sans Mono", va="center")
    ax.set_xticks(range(4)); ax.set_xlabel("Number of D (2,6-diaminopurine) bases in codon")
    ax.set_ylabel("Mean Vina affinity over 15 amino acids (kcal/mol)")
    ax.set_title(f"Codon ranking follows DAP content (ρ = {rho_dap.statistic:.2f})")
    fig.savefig(f"{a.figures}/fig4_codon_dap_effect.png"); plt.close(fig)

    # fig5 MM-GBSA ranked with SD band
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    s = d.sort_values("mmgbsa_dg_bind_mean_kcal_mol").reset_index(drop=True)
    x = np.arange(len(s))
    ax.fill_between(x, s.mmgbsa_dg_bind_mean_kcal_mol - s.mmgbsa_dg_bind_std_kcal_mol,
                    s.mmgbsa_dg_bind_mean_kcal_mol + s.mmgbsa_dg_bind_std_kcal_mol, color=SEQ[0], linewidth=0,
                    label="± 1 SD over 50 frames")
    ax.plot(x, s.mmgbsa_dg_bind_mean_kcal_mol, color=BLUE, linewidth=2, label="Mean ΔG")
    ax.axhline(0, color=INK2, linewidth=1)
    ax.set_xlim(0, len(s) - 1); ax.set_xlabel("960 codon–amino-acid pairs, ranked"); ax.set_ylabel("MM-GBSA ΔG (kcal/mol)")
    ax.set_title("Every MM-GBSA estimate overlaps zero")
    ax.legend(frameon=False, loc="lower right", fontsize=8.5)
    fig.savefig(f"{a.figures}/fig5_mmgbsa_ranked.png"); plt.close(fig)

    # fig6 chirality pairs (Vina)
    fig, axes = plt.subplots(1, 4, figsize=(11, 3.0), sharex=False)
    for ax, (l, dd_) in zip(axes, DL_PAIRS.items()):
        lo = min(V[l].min(), V[dd_].min()) - 0.05; hi = max(V[l].max(), V[dd_].max()) + 0.05
        ax.plot([lo, hi], [lo, hi], color=MUTED, linewidth=1, linestyle="--")
        ax.scatter(V[l], V[dd_], s=16, color=BLUE, alpha=0.8, edgecolor="white", linewidth=0.6, zorder=3)
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        r = chir[(chir.pair == f"{l}/{dd_}") & (chir.metric == "vina")].iloc[0]
        ax.set_title(f"{AA_NAMES[l]} vs {AA_NAMES[dd_]}", fontsize=10)
        ax.text(0.04, 0.96, f"Δmedian {r.median_D_minus_L:+.3f}\np = {r.wilcoxon_p:.3f}", transform=ax.transAxes,
                va="top", fontsize=7.5, color=INK2)
        ax.set_xlabel("L (kcal/mol)", fontsize=8.5)
        ax.tick_params(labelsize=7.5)
    axes[0].set_ylabel("D (kcal/mol)", fontsize=8.5)
    fig.suptitle("Enantiomers dock alike on every codon (Vina, 64 codons per panel)", x=0.01, ha="left",
                 fontweight="bold", fontsize=11, color=INK)
    fig.tight_layout()
    fig.savefig(f"{a.figures}/fig6_chirality_pairs.png"); plt.close(fig)

    print(json.dumps(summ, indent=2, default=float))


if __name__ == "__main__":
    main()
