# Stereochemical recognition between DMHU-PNA codons and Bennu amino acids: a 960-pair computational screen

Pelinsu Eda Şentürk · Bioengineering, Istanbul Medeniyet University 

*[Türkçe sürüm → rapor_tr.md](rapor_tr.md)*

## Abstract

In this model, free amino acids show no codon-specific recognition by single-stranded DMHU-PNA codons.

**Question.** Could a peptide nucleic acid (PNA) genetic system built from four non-canonical bases assign amino acids to codons through direct stereochemical interaction? The bases are 2,6-diaminopurine (D), 6-methyluracil (M), hypoxanthine (H) and uracil (U).

**Methods.** All 64 three-base PNA codons were screened against 15 amino acids reported in Bennu samples and carbonaceous meteorites, including D-enantiomers, Aib, isovaline, β-alanine and GABA. Each of the 960 pairs went through AutoDock Vina docking, 100 ps heating and 150 ps unrestrained MD in implicit solvent (GB, igb=5), and single-trajectory MM-GBSA over 50 frames. Selected pairs were re-simulated independently in OpenMM.

**Results.**

1. 92 % of Vina variance is additive: 49 % amino acid, 43 % codon. The amino-acid term tracks ligand size (ρ = −0.94) and the codon term tracks 2,6-diaminopurine count (ρ = −0.80). Codon × amino-acid specificity is 8 % (SD 0.09 kcal/mol).
2. Every MM-GBSA mean lies within 2 SD of zero. Ligands leave the trimer during MD, and top-ranked pairs do not reproduce.
3. No L/D difference exceeds method accuracy.

**Conclusion.** Under these conditions, direct triplet–amino-acid recognition is not supported as a coding mechanism. Had such a system existed, assignment through adaptor-like intermediates would be more plausible.

## Introduction

The stereochemical hypothesis holds that the first codon assignments arose from direct chemical affinity between amino acids and nucleotide triplets. RNA aptamer studies reported cognate codons and anticodons in amino-acid binding sites more often than chance (Yarus et al., 2009). The hypothesis remains debated, and alternatives such as coevolution or a "frozen accident" need no direct affinity (Koonin and Novozhilov, 2009).

PNA is a nucleic-acid analogue with a neutral N-(2-aminoethyl)glycine backbone instead of sugar-phosphate (Nielsen et al., 1991). It is discussed as a candidate prebiotic genetic polymer. Samples returned from asteroid Bennu by OSIRIS-REx contain 14 of the 20 protein amino acids and all five DNA/RNA nucleobases. Chiral non-protein amino acids in them were racemic or nearly so (Glavin et al., 2025), so D- and L-forms may have co-existed prebiotically.

This study asked: in a PNA genetic system built from D, M, H and U, do Bennu-relevant amino acids bind their codons directly and specifically? If they do, does binding distinguish enantiomers?

## Methods

All calculations used AmberTools, Open Babel and AutoDock Vina on Google Colab. The code for each step is in `scripts/`.

| Step | Tool | Details |
| --- | --- | --- |
| Monomer parametrisation | antechamber, parmchk2 | Capped PNA monomers; GAFF2 atom types, AM1-BCC charges; caps removed to make tleap residues (DAP, MTU, HYX, URA) |
| Codon trimers | tleap, sander | 4³ = 64 trimers; 2000-step GB minimisation (igb=5) |
| Amino acids | antechamber, tleap | 15 amino acids; GAFF2/AM1-BCC; net charge 0 (neutral COOH/NH₂ form) |
| Docking | Vina 1.2 (Python) | Minimised trimer as rigid receptor; box = receptor bounding box + 6 Å; exhaustiveness 8; top pose kept |
| Complex | ParmEd | Pose mapped onto the ligand topology by atom name; complex = receptor + ligand |
| MD | sander, igb=5, no cutoff | 500-step minimisation (heavy atoms restrained, 5 kcal/mol/Å²) → 100 ps heating 0→300 K (2 kcal/mol/Å²) → 150 ps unrestrained MD; Langevin γ = 5 ps⁻¹, 2 fs, SHAKE |
| MM-GBSA | sander single point | One frame every 3 ps (50 frames): ΔG = E(complex) − E(receptor) − E(ligand); no entropy term |

The 64 × 15 Vina and MM-GBSA matrices were decomposed with a two-way additive model: grand mean + codon effect + amino-acid effect + residual. Codon-specific recognition would have to appear in the residual. The size effect was tested by Spearman correlation between amino-acid means and heavy-atom count. Chirality was tested for four L/D pairs (Ala, Asp, Glu, Val) with paired Wilcoxon tests over the 64 codons, Bonferroni-corrected.

## Results

All 960 pairs completed without error. Vina affinities range from −3.51 to −1.76 kcal/mol (mean −2.50), which corresponds to millimolar-level binding on Vina's nominal scale. The strongest scores are DMH–D-Gln (−3.51), DMH–L-Glu (−3.49) and DMH–L-Asp (−3.48).

![Figure 1. Vina affinity for 64 codons × 15 amino acids, both axes sorted by mean.](../figures/fig1_vina_heatmap.png)

*Figure 1. Vina affinity for 64 codons × 15 amino acids, both axes sorted by mean. Row and column patterns are independent: each amino acid keeps roughly the same rank on every codon.*

### Docking scores are additive and size-driven

The amino-acid term explains 49 % of Vina variance and the codon term 43 %. The residual, which would carry codon-specific recognition, is only 8 % (SD 0.09 kcal/mol). That is far below Vina's typical error of a few kcal/mol.

![Figure 2. Two-way additive decomposition of Vina and MM-GBSA variance.](../figures/fig2_variance_decomposition.png)

*Figure 2. Two-way additive decomposition of Vina and MM-GBSA variance.*

The amino-acid term comes almost entirely from ligand size: Spearman ρ = −0.94 (p < 10⁻⁶) between mean Vina score and heavy-atom count. The codon term tracks the number of 2,6-diaminopurine (D) bases (ρ = −0.80). D is the largest of the four bases and has the most H-bond donors.

| Amino-acid size | D bases per codon |
| --- | --- |
| ![Figure 3](../figures/fig3_aa_size_effect.png) | ![Figure 4](../figures/fig4_codon_dap_effect.png) |

*Figures 3 and 4. Mean Vina score per amino acid vs heavy-atom count (left); mean Vina score per codon vs number of D bases (right).*

### MM-GBSA shows no binding signal

MM-GBSA means range from −2.48 to +0.45 kcal/mol, and ΔG is positive for 61 % of pairs. For all 960 pairs, |mean| < 2 SD. Mean and SD are strongly anti-correlated (ρ = −0.69), so negative ΔG values come only from highly fluctuating runs. MM-GBSA and Vina rankings are weakly and inversely related (ρ = −0.23). The most negative pairs are UMH–L-Iva (−2.48 ± 2.60), MDU–L-Met (−2.14 ± 1.51) and HHH–L-Iva (−1.71 ± 1.53).

![Figure 5. MM-GBSA ΔG for all 960 pairs, ranked; light band ± 1 SD.](../figures/fig5_mmgbsa_ranked.png)

*Figure 5. MM-GBSA ΔG for all 960 pairs, ranked; light band ± 1 SD.*

### No chirality preference

L/D pairs score almost identically on all 64 codons.

| Pair | Vina, median D − L (kcal/mol) | Vina p (Bonferroni) | MM-GBSA, median D − L (kcal/mol) | MM-GBSA p (Bonferroni) |
| --- | --- | --- | --- | --- |
| L-Ala / D-Ala | +0.016 | 0.59 | +0.013 | 1.00 |
| L-Asp / D-Asp | +0.027 | 1.00 | −0.006 | 1.00 |
| L-Glu / D-Glu | +0.049 | 0.039 | −0.008 | 1.00 |
| L-Val / D-Val | +0.014 | 1.00 | −0.008 | 1.00 |

Only the Glu/D-Glu Vina difference survives correction. It is 0.05 kcal/mol, far below method resolution, and absent from MM-GBSA.

![Figure 6. Vina scores of L and D enantiomers on each codon; dashed line y = x.](../figures/fig6_chirality_pairs.png)

*Figure 6. Vina scores of L and D enantiomers on each codon; dashed line y = x.*

## Independent validation

In OpenMM re-runs, the ligand leaves the trimer within a few tens of picoseconds in most pairs. MM-GBSA ΔG therefore measures how long the ligand stays near the trimer, not affinity.

The screen protocol was rebuilt in OpenMM 8.6: OBC2 (the counterpart of igb=5), mbondi radii, no cutoff, and the same restraints, heating and production lengths. The starting pose came from a random rigid-body search on the GB interaction energy instead of Vina. Nine pairs were chosen: the four with the most negative MM-GBSA, the strongest Vina score (DMH–D-Gln) and four pairs near zero (DDD–L-Asp, DDD–L-Met, HHH–Gly, UUU–L-Glu). The top two pairs were repeated with four seeds each. Code and results are in `validation/`.

![Figure 7. Minimum ligand–trimer heavy-atom distance during production MD (9 pairs).](../figures/fig7_validation_dissociation.png)

*Figure 7. Minimum ligand–trimer heavy-atom distance during production MD (9 pairs). Dashed line: 4 Å contact threshold.*

The fraction of frames in contact ranges from 2 % to 56 %. At the start of MD the ligand touches the trimer with an interaction energy of −1.2 to −3.9 kcal/mol; that contact is not kept once restraints are lifted. For near-zero pairs the re-runs match the original table, e.g. DDD–L-Asp 0.12 ± 0.27 originally and 0.06 ± 0.13 kcal/mol in the re-run.

![Figure 8. Original table values vs OpenMM re-runs.](../figures/fig8_validation_reproducibility.png)

*Figure 8. Original table values (open diamonds) vs OpenMM re-runs (blue dots, one per seed).*

| Pair | Original ΔG (kcal/mol) | Re-runs ΔG (kcal/mol) | Frames in contact |
| --- | --- | --- | --- |
| UMH–L-Iva | −2.48 ± 2.60 | −0.64; −0.10; −0.15; −0.25 | 8–30 % |
| MDU–L-Met | −2.14 ± 1.51 | +0.09; +0.11; −0.28; −0.08 | 4–38 % |
| HHH–L-Iva | −1.71 ± 1.53 | −0.16 | 36 % |
| MMH–D-Val | −1.71 ± 1.72 | −0.31 | 20 % |
| DDD–L-Met | −0.34 ± 0.86 | −0.58 | 56 % |

The top of the table does not reproduce. That is what regression to the mean predicts when the extremes are picked from 960 noisy estimates. These pairs rank high because the ligand happened to stay in contact longer. Hydrophobic ligands (Iva, Met, Val) appear often among them: van der Waals contact makes them slightly stickier, but not specifically bound.

## Discussion and limitations

The interaction between a free amino acid and a single-stranded PNA triplet is both weak and non-specific. What drives the docking score is not "this codon recognises this amino acid" but "a bigger ligand touches more surface", consistent with Vina's known size bias. In RNA aptamer studies, amino-acid selectivity appears in pockets of folded, much larger structures (Yarus et al., 2009), so the absence of a signal on bare triplets does not contradict that work.

The negative result holds within the limits below. "Not detected" is not the same as "impossible".

- **Receptor geometry.** Docking used the minimised trimer, but the complex was built on unminimised tleap coordinates (RMSD 0.4–1.1 Å).
- **Open trimer termini.** No terminal groups were added after cap removal; the first backbone N and last carbonyl C have two bonds each.
- **Neutral amino acids.** Ligands were modelled in the neutral form, not as the zwitterions and anions present at pH 7. This can change electrostatics.
- **GB radii.** mbondi was used instead of the mbondi2 set recommended for igb=5.
- **Short, unrestrained MD.** The weak complex falls apart within 150 ps. Single-trajectory MM-GBSA without entropy then does not measure binding free energy.
- **Model scope.** Single-stranded trimer, implicit solvent, one pose per pair, no ions or mineral surfaces, a general force field (GAFF2) for non-standard PNA residues.

The validation re-runs avoid the first two issues (consistent receptor coordinates, pose searched on the same geometry) and behave the same way, so these issues are not expected to change the main result. The effect of zwitterionic forms was not tested.

## Conclusion and next steps

In this model, no evidence was found that Bennu-relevant amino acids bind DMHU-PNA codons directly, specifically or chirally. Had such a genetic system existed, codon–amino-acid assignment would more plausibly have relied on adaptor-like intermediates than on direct triplet binding. The study also makes two pitfalls of high-throughput docking + MM-GBSA screens concrete: ligand-size bias, and ΔG measuring residence time when complexes are unstable.

1. **Technical fixes.** Consistent receptor coordinates, capped trimer termini, zwitterionic forms, mbondi2.
2. **A suitable free-energy method.** Restrained MD with umbrella sampling or alchemical ΔG for weak binders, or many short explicit-solvent replicas per pair.
3. **Larger structures.** PNA duplexes, hairpins and loops where binding pockets can form.
4. **Size control.** Use docking only as a filter and add size-normalised measures such as ligand efficiency.

## References

1. Glavin D. P. et al. Abundant ammonia and nitrogen-rich soluble organic matter in samples from asteroid (101955) Bennu. *Nature Astronomy* (2025). [doi:10.1038/s41550-024-02472-9](https://doi.org/10.1038/s41550-024-02472-9)
2. Yarus M., Widmann J. J., Knight R. RNA–amino acid binding: a stereochemical era for the genetic code. *J. Mol. Evol.* 69, 406–429 (2009). [doi:10.1007/s00239-009-9270-1](https://doi.org/10.1007/s00239-009-9270-1)
3. Koonin E. V., Novozhilov A. S. Origin and evolution of the genetic code: the universal enigma. *IUBMB Life* 61, 99–111 (2009). [arXiv:0807.4749](https://arxiv.org/abs/0807.4749)
4. Nielsen P. E., Egholm M., Berg R. H., Buchardt O. Sequence-selective recognition of DNA by strand displacement with a thymine-substituted polyamide. *Science* 254, 1497–1500 (1991).
5. Eberhardt J., Santos-Martins D., Tillack A. F., Forli S. AutoDock Vina 1.2.0. *J. Chem. Inf. Model.* 61, 3891–3898 (2021).
6. Onufriev A., Bashford D., Case D. A. Exploring protein native states and large-scale conformational changes with a modified generalized Born model. *Proteins* 55, 383–394 (2004).
7. Genheden S., Ryde U. The MM/PBSA and MM/GBSA methods to estimate ligand-binding affinities. *Expert Opin. Drug Discov.* 10, 449–461 (2015).
