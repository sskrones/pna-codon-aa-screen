"""Step 4 - GAFF2 / AM1-BCC parametrisation of the 15 Bennu-relevant amino acids.

Inputs  (data/raw/amino_acids/, not yet in the repo): one SDF per amino acid.
Outputs (build/amino_acids/): {code}.gaff2.mol2, {code}.frcmod, {code}.prmtop, {code}.inpcrd

All amino acids are built with net charge 0 (-nc 0), i.e. in the neutral
COOH / NH2 form, not the zwitterion found in water at pH 7 (see README, limitations).
"""
import os
import subprocess

RAW = os.environ.get("RAW", "data/raw/amino_acids")
OUT = os.environ.get("OUT", "build/amino_acids")
os.makedirs(OUT, exist_ok=True)

AA_LIST = [  # (SDF file, residue code)
    ("alpha aib.sdf", "AIB"),          # 2-aminoisobutyric acid
    ("beta alanine.sdf", "BAL"),       # beta-alanine
    ("d-alanine.sdf", "DAL"),
    ("d-aspartic acid.sdf", "DAS"),
    ("d-glutamic acid.sdf", "DGU"),
    ("d-glutamine.sdf", "DGN"),
    ("d-valine.sdf", "DVA"),
    ("gamma aba.sdf", "GAB"),          # gamma-aminobutyric acid
    ("glycine.sdf", "GLY"),
    ("l-alanine.sdf", "ALA"),
    ("l-aspartic acid.sdf", "ASP"),
    ("l-glutamic acid.sdf", "GLU"),
    ("l-isovaline.sdf", "ISV"),
    ("l-methionine.sdf", "MET"),
    ("l-valine.sdf", "VAL"),
]

for sdf, code in AA_LIST:
    c = code.lower()
    subprocess.run(["antechamber", "-i", os.path.join(RAW, sdf), "-fi", "mdl",
                    "-o", f"{OUT}/{c}.gaff2.mol2", "-fo", "mol2",
                    "-c", "bcc", "-at", "gaff2", "-nc", "0", "-m", "1",
                    "-rn", code, "-pf", "y"], check=True)
    subprocess.run(["parmchk2", "-i", f"{OUT}/{c}.gaff2.mol2", "-f", "mol2",
                    "-o", f"{OUT}/{c}.frcmod", "-s", "gaff2"], check=True)
    leap = f"""source leaprc.gaff2
loadamberparams {OUT}/{c}.frcmod
mol = loadmol2 {OUT}/{c}.gaff2.mol2
check mol
saveamberparm mol {OUT}/{c}.prmtop {OUT}/{c}.inpcrd
quit
"""
    with open(f"{OUT}/{c}.leap", "w") as f:
        f.write(leap)
    subprocess.run(["tleap", "-f", f"{OUT}/{c}.leap"], check=True)
    print(code, "done")
