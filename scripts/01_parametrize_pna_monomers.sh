#!/usr/bin/env bash
# Step 1 - GAFF2 / AM1-BCC parametrisation of the four PNA monomers.
#
# Inputs  (data/raw/, not yet in the repo): {D,H,M,U}-capped_monomer.cjson
#         capped PNA monomers built in Avogadro (backbone capped on both ends).
# Outputs (build/monomers/): *-capped_monomer.gaff2.mol2, *-capped_monomer.frcmod
#
# Base letters / residue names:
#   D = 2,6-diaminopurine (DAP)   M = 6-methyluracil (MTU)
#   H = hypoxanthine      (HYX)   U = uracil         (URA)
set -euo pipefail
RAW=${RAW:-data/raw}
OUT=${OUT:-build/monomers}
mkdir -p "$OUT"

declare -A RESNAME=( [D]=DAP [M]=MTU [H]=HYX [U]=URA )

for b in D M H U; do
  obabel "$RAW/${b}-capped_monomer.cjson" -O "$OUT/${b}-capped_monomer.mol2"
  antechamber -i "$OUT/${b}-capped_monomer.mol2" -fi mol2 \
              -o "$OUT/${b}-capped_monomer.gaff2.mol2" -fo mol2 \
              -c bcc -at gaff2 -nc 0 -m 1 -rn "${RESNAME[$b]}" -pf y
  parmchk2 -i "$OUT/${b}-capped_monomer.gaff2.mol2" -f mol2 \
           -o "$OUT/${b}-capped_monomer.frcmod" -s gaff2
done
echo "Monomer parametrisation done -> $OUT"
