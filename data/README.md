# Data

| Folder | Files | Notes |
|---|---|---|
| `receptors/` | `{CODON}.prmtop`, `{CODON}.inpcrd` (64 each) | tleap output, **unminimised** coordinates, mbondi radii |
| | `{codon}_min.pdbqt` (64) | GB-minimised trimers, used as rigid Vina receptors |
| `ligands/` | `{aa}.prmtop`, `{aa}.inpcrd`, `{aa}.pdbqt` (15 each) | GAFF2 + AM1-BCC, net charge 0 (neutral, non-zwitterionic) |
| `mdin/` | `min_complex.in`, `md_equil.in`, `md_prod.in`, `sp_energy.in` | sander inputs used by `scripts/06_dock_md_mmgbsa.py` |
| `raw/` | *(to add)* `{D,H,M,U}-capped_monomer.cjson`, amino-acid `.sdf` files | inputs of steps 1 and 4 |

Codon letters: D = 2,6-diaminopurine (DAP), M = 6-methyluracil (MTU), H = hypoxanthine (HYX), U = uracil (URA).
Amino-acid codes: AIB α-aminoisobutyric acid, ALA L-Ala, ASP L-Asp, BAL β-alanine, DAL D-Ala, DAS D-Asp,
DGN D-Gln, DGU D-Glu, DVA D-Val, GAB γ-aminobutyric acid, GLU L-Glu, GLY Gly, ISV L-isovaline, MET L-Met, VAL L-Val.
