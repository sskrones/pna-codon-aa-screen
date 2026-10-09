"""Step 5 - convert amino acids (ligands) and minimised trimers (receptors) to PDBQT.

GAFF2 atom types are mapped to SYBYL types so Open Babel perceives elements and
hybridisation correctly; AM1-BCC charges are written back into the PDBQT charge
column; hydrogens bonded to N/O/S are typed HD (donor), the rest H; ring carbons
of the bases are typed A (aromatic) from their original GAFF2 type.
All hydrogens are kept so that every PDBQT atom name matches the prmtop atom name
(needed later to put the docked pose back into the Amber topology).

Usage:  python 05_make_pdbqt.py --aa-dir build/amino_acids --trimer-dir build/trimers
"""
import argparse
import itertools
import subprocess

GAFF2_TO_SYBYL = {
    "oh": "O.3", "o": "O.2", "n8": "N.3", "nt": "N.am",
    "c3": "C.3", "c": "C.2",
    "hc": "H", "h1": "H", "hn": "H", "ho": "H",
    "ss": "S.3",
    "n": "N.3", "n2": "N.2", "n3": "N.3", "na": "N.ar", "nb": "N.ar",
    "nc": "N.ar", "nd": "N.ar", "ns": "N.ar", "nv": "N.ar",
    "ca": "C.ar", "cc": "C.2", "cd": "C.2", "cq": "C.ar", "cp": "C.ar",
    "cf": "C.2", "ce": "C.2",
    "os": "O.3", "op": "O.2", "oq": "O.2",
    "h4": "H", "h5": "H", "ha": "H",
}
AROMATIC_C_TYPES = {"ca", "cc", "cd", "cp", "cq"}   # ce/cf are conjugated chain carbons
AA_CODES = ["aib", "ala", "asp", "bal", "dal", "das", "dgn", "dgu",
            "dva", "gab", "glu", "gly", "isv", "met", "val"]


def retype_mol2(mol2_in, mol2_out):
    """Rewrite GAFF2 atom types as SYBYL; return per-atom info for later fixes."""
    names, resids, charges, gtypes, bonds = {}, {}, {}, {}, []
    out, section, unmapped = [], None, set()
    for line in open(mol2_in):
        if line.startswith("@<TRIPOS>ATOM"):
            section = "ATOM"; out.append(line); continue
        if line.startswith("@<TRIPOS>BOND"):
            section = "BOND"; out.append(line); continue
        if line.startswith("@<TRIPOS>"):
            section = None; out.append(line); continue
        if section == "ATOM" and line.strip():
            p = line.split()
            aid, name, gt, resid, q = int(p[0]), p[1], p[5], int(p[6]), float(p[8])
            names[aid], resids[aid] = name, resid
            charges[(resid, name)], gtypes[(resid, name)] = q, gt
            if gt not in GAFF2_TO_SYBYL:
                unmapped.add(gt)
            out.append(f"{aid:>7} {name:<10}{float(p[2]):>10.4f}{float(p[3]):>10.4f}"
                       f"{float(p[4]):>10.4f} {GAFF2_TO_SYBYL.get(gt, gt):<10}{resid:>5}"
                       f" {p[7]:<8}{q:>12.6f}\n")
        else:
            if section == "BOND" and line.strip():
                p = line.split(); bonds.append((int(p[1]), int(p[2])))
            out.append(line)
    if unmapped:
        raise ValueError(f"{mol2_in}: unmapped GAFF2 types {unmapped}")
    open(mol2_out, "w").writelines(out)
    h_parent = {}
    for a1, a2 in bonds:
        n1, n2 = names[a1], names[a2]
        if n1.startswith("H") and not n2.startswith("H"):
            h_parent[(resids[a1], n1)] = n2[0]
        elif n2.startswith("H") and not n1.startswith("H"):
            h_parent[(resids[a2], n2)] = n1[0]
    aromatic = {k for k, t in gtypes.items() if t in AROMATIC_C_TYPES}
    return charges, h_parent, aromatic


def fix_pdbqt(pdbqt, charges, h_parent, aromatic):
    out = []
    for line in open(pdbqt):
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            try:
                key = (int(line[22:26]), name)
            except ValueError:
                key = (None, name)
            if key in charges:
                line = line[:70] + f"{charges[key]:+.3f}" + line[76:]
            if name.startswith("H") and key in h_parent:
                rest, _ = line.rstrip("\n").rsplit(None, 1)
                line = f"{rest} {'HD' if h_parent[key] in ('N', 'O', 'S') else 'H'}\n"
            elif key in aromatic:
                rest, _ = line.rstrip("\n").rsplit(None, 1)
                line = f"{rest} A\n"
        out.append(line)
    open(pdbqt, "w").writelines(out)


def ligands(aa_dir):
    for c in AA_CODES:
        info = retype_mol2(f"{aa_dir}/{c}.gaff2.mol2", f"{aa_dir}/{c}.fixed.mol2")
        subprocess.run(["obabel", f"{aa_dir}/{c}.fixed.mol2", "-O", f"{aa_dir}/{c}.pdbqt", "-xh"], check=True)
        fix_pdbqt(f"{aa_dir}/{c}.pdbqt", *info)
        print("ligand", c, "done")


def receptors(trimer_dir):
    for combo in itertools.product("dmhu", repeat=3):
        codon = "".join(combo); C = codon.upper()
        mol2 = f"{trimer_dir}/{codon}_min.mol2"
        with open(f"{trimer_dir}/{codon}_cpptraj.in", "w") as f:
            f.write(f"parm {trimer_dir}/{C}.prmtop\ntrajin {trimer_dir}/{C}_min.rst\ntrajout {mol2} mol2\nrun\n")
        subprocess.run(["cpptraj", "-i", f"{trimer_dir}/{codon}_cpptraj.in"], check=True)
        info = retype_mol2(mol2, f"{trimer_dir}/{codon}_fixed.mol2")
        # -xr: rigid receptor (no torsion tree)
        subprocess.run(["obabel", f"{trimer_dir}/{codon}_fixed.mol2", "-O",
                        f"{trimer_dir}/{codon}_min.pdbqt", "-xr", "-xh"], check=True)
        fix_pdbqt(f"{trimer_dir}/{codon}_min.pdbqt", *info)
        print("receptor", C, "done")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--aa-dir", default="build/amino_acids")
    ap.add_argument("--trimer-dir", default="build/trimers")
    a = ap.parse_args()
    ligands(a.aa_dir)
    receptors(a.trimer_dir)
