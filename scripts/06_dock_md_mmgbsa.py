"""Step 6 - main screen: Vina docking -> complex assembly -> GB MD -> single-trajectory MM-GBSA
for all 64 codons x 15 amino acids (960 pairs).

This is the final pipeline cell of notebooks/parametrizasyon_colab_export.py,
with Colab/Drive code removed and paths made configurable. The protocol is
unchanged, so it reproduces results/mmgbsa_results_md.csv (up to MD noise).

Per pair:
  1. Vina docking of the ligand onto {codon}_min.pdbqt
     (box = receptor bounding box + 6 A padding, exhaustiveness 8, top pose kept)
  2. Pose copied onto the ligand Amber topology by atom name; complex = receptor + ligand (ParmEd)
  3. sander, igb=5, no cutoff:
        min_complex.in  500 steps, heavy atoms restrained (5 kcal/mol/A^2)
        md_equil.in     100 ps, 0 -> 300 K, heavy atoms restrained (2 kcal/mol/A^2)
        md_prod.in      150 ps unrestrained, frame every 3 ps (50 frames)
  4. For each frame: dG = E(complex) - E(receptor) - E(ligand)  (sp_energy.in, igb=5)
     -> mean and SD over frames. No entropy term.

KNOWN ISSUES (kept on purpose so the published numbers are reproducible; see README):
  * receptor coordinates come from {CODON}.inpcrd (unminimised tleap geometry)
    while docking used {codon}_min.pdbqt (minimised geometry);
  * production MD is unrestrained and the weakly bound ligand usually leaves
    the trimer within tens of ps, so dG mostly reflects residence time;
  * prmtops carry mbondi radii (mbondi2 is recommended for igb=5).

Usage:
  python scripts/06_dock_md_mmgbsa.py --receptors data/receptors --ligands data/ligands \
         --mdin data/mdin --out results/mmgbsa_results_md.csv
"""
import argparse
import csv
import glob
import os
import shutil
import subprocess
import tempfile

import numpy as np
import parmed as pmd
from vina import Vina

AA_CODES = ["aib", "ala", "asp", "bal", "dal", "das", "dgn", "dgu",
            "dva", "gab", "glu", "gly", "isv", "met", "val"]
FIELDS = ["codon", "amino_acid", "vina_affinity_kcal_mol",
          "mmgbsa_dg_bind_mean_kcal_mol", "mmgbsa_dg_bind_std_kcal_mol",
          "n_snapshots", "status"]


def receptor_box(pdbqt_path, padding=6.0):
    xyz = np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])]
                    for l in open(pdbqt_path) if l.startswith(("ATOM", "HETATM"))])
    lo, hi = xyz.min(0), xyz.max(0)
    return ((lo + hi) / 2).tolist(), (hi - lo + 2 * padding).tolist()


def pose_coords(pdbqt_path):
    return {l[12:16].strip(): (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            for l in open(pdbqt_path) if l.startswith(("ATOM", "HETATM"))}


def sander(mdin, prmtop, inpcrd, out, rst, ref=None, traj=None):
    cmd = ["sander", "-O", "-i", mdin, "-p", prmtop, "-c", inpcrd, "-r", rst, "-o", out]
    if traj:
        cmd += ["-x", traj]
    if ref:
        cmd += ["-ref", ref]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"sander failed ({mdin}): {r.stdout[-400:]}{r.stderr[-400:]}")


def last_energy(log):
    """Total energy from the last 'NSTEP ENERGY' block of a sander minimisation log."""
    lines = open(log).read().splitlines()
    idx = [i for i, l in enumerate(lines) if "NSTEP" in l and "ENERGY" in l]
    if not idx:
        return None
    for l in lines[idx[-1] + 1:]:
        if l.strip():
            return float(l.split()[1])
    return None


def run_pair(codon, aa, a, v, receptor, n_rec, tmp):
    j = lambda name: os.path.join(tmp, name)
    v.set_ligand_from_file(os.path.join(a.ligands, f"{aa}.pdbqt"))
    v.dock(exhaustiveness=a.exhaustiveness, n_poses=9)
    affinity = v.energies(n_poses=1)[0][0]
    v.write_poses(j("pose.pdbqt"), n_poses=1, overwrite=True)

    pose = pose_coords(j("pose.pdbqt"))
    lig = pmd.load_file(os.path.join(a.ligands, f"{aa}.prmtop"))
    missing = [at.name for at in lig.atoms if at.name not in pose]
    if missing:
        raise ValueError(f"atom names not in pose: {missing}")
    lig.coordinates = [pose[at.name] for at in lig.atoms]
    cpx = receptor + lig
    cpx.save(j("cpx.prmtop"), overwrite=True)
    cpx.save(j("cpx.inpcrd"), overwrite=True)

    mdin = lambda n: os.path.join(a.mdin, n)
    sander(mdin("min_complex.in"), j("cpx.prmtop"), j("cpx.inpcrd"), j("min.out"), j("min.rst7"), ref=j("cpx.inpcrd"))
    sander(mdin("md_equil.in"), j("cpx.prmtop"), j("min.rst7"), j("eq.out"), j("eq.rst7"), ref=j("min.rst7"))
    sander(mdin("md_prod.in"), j("cpx.prmtop"), j("eq.rst7"), j("prod.out"), j("prod.rst7"), traj=j("prod.nc"))

    with open(j("extract.cpptraj"), "w") as f:
        f.write(f"parm {j('cpx.prmtop')}\ntrajin {j('prod.nc')}\ntrajout {j('snap.rst7')} restart\n")
    subprocess.run(["cpptraj", "-i", j("extract.cpptraj")], check=True, capture_output=True)
    snaps = sorted(glob.glob(j("snap.rst7*")), key=lambda p: int(p.rsplit(".", 1)[-1]) if p[-1].isdigit() else 0)

    rec_prm = os.path.join(a.receptors, f"{codon}.prmtop")
    lig_prm = os.path.join(a.ligands, f"{aa}.prmtop")
    dg = []
    for s in snaps:
        xyz = pmd.load_file(j("cpx.prmtop"), s).coordinates
        r = pmd.load_file(rec_prm); r.coordinates = xyz[:n_rec]; r.save(j("r.rst7"), overwrite=True)
        l = pmd.load_file(lig_prm); l.coordinates = xyz[n_rec:]; l.save(j("l.rst7"), overwrite=True)
        sp = mdin("sp_energy.in")
        sander(sp, j("cpx.prmtop"), s, j("c_sp.out"), j("x.rst7"))
        sander(sp, rec_prm, j("r.rst7"), j("r_sp.out"), j("x.rst7"))
        sander(sp, lig_prm, j("l.rst7"), j("l_sp.out"), j("x.rst7"))
        e = [last_energy(j(n)) for n in ("c_sp.out", "r_sp.out", "l_sp.out")]
        if None not in e:
            dg.append(e[0] - e[1] - e[2])
    if not dg:
        raise RuntimeError("no snapshot energies could be read")
    dg = np.array(dg)
    return affinity, dg.mean(), dg.std(), len(dg)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--receptors", default="data/receptors")
    ap.add_argument("--ligands", default="data/ligands")
    ap.add_argument("--mdin", default="data/mdin")
    ap.add_argument("--out", default="results/mmgbsa_results_md.csv")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--codons", nargs="*", help="subset of codons, e.g. DDD DMH")
    ap.add_argument("--aas", nargs="*", help="subset of amino acids, e.g. asp gly")
    a = ap.parse_args()

    for tool in ("sander", "cpptraj"):
        if shutil.which(tool) is None:
            raise SystemExit(f"{tool} not found - install AmberTools (see environment.yml)")

    receptors = sorted(glob.glob(os.path.join(a.receptors, "*_min.pdbqt")))
    if a.codons:
        receptors = [p for p in receptors if os.path.basename(p)[:3].upper() in {c.upper() for c in a.codons}]
    aas = [x.lower() for x in a.aas] if a.aas else AA_CODES

    done = set()
    if os.path.exists(a.out):
        done = {(r["codon"], r["amino_acid"]) for r in csv.DictReader(open(a.out))}
        print(f"{len(done)} pairs already in {a.out}; resuming")
    new = not os.path.exists(a.out)
    fh = open(a.out, "a", newline="")
    w = csv.DictWriter(fh, fieldnames=FIELDS)
    if new:
        w.writeheader()

    for rp in receptors:
        codon = os.path.basename(rp)[:3].upper()
        todo = [x for x in aas if (codon, x.upper()) not in done]
        if not todo:
            continue
        center, size = receptor_box(rp)
        v = Vina(sf_name="vina", verbosity=0)
        v.set_receptor(rp)
        v.compute_vina_maps(center=center, box_size=size)
        receptor = pmd.load_file(os.path.join(a.receptors, f"{codon}.prmtop"),
                                 os.path.join(a.receptors, f"{codon}.inpcrd"))
        for aa in todo:
            row = {"codon": codon, "amino_acid": aa.upper()}
            tmp = tempfile.mkdtemp(prefix=f"{codon}_{aa}_")
            try:
                aff, m, s, n = run_pair(codon, aa, a, v, receptor, len(receptor.atoms), tmp)
                row.update(vina_affinity_kcal_mol=aff, mmgbsa_dg_bind_mean_kcal_mol=m,
                           mmgbsa_dg_bind_std_kcal_mol=s, n_snapshots=n, status="OK")
                print(f"{codon}-{aa.upper()}: vina={aff:.3f}  dG={m:.3f} +/- {s:.3f}  n={n}")
            except Exception as e:  # keep going, record the failure
                row.update(n_snapshots=0, status=f"ERROR: {e}")
                print(f"{codon}-{aa.upper()}: ERROR {e}")
            finally:
                shutil.rmtree(tmp, ignore_errors=True)
            w.writerow(row)
            fh.flush()
    fh.close()


if __name__ == "__main__":
    main()
