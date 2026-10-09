"""Step 3 - build all 64 PNA codon trimers with tleap and minimise them in GB solvent.

Run inside build/monomers (needs D/M/H/U.lib and *-capped_monomer.frcmod).
Outputs: build/trimers/{CODON}.prmtop / .inpcrd   (tleap geometry)
         build/trimers/{CODON}_min.rst / _min.pdb  (sander-minimised geometry)

NOTE: data/receptors/{CODON}.inpcrd in this repo is the *unminimised* tleap
geometry, while {codon}_min.pdbqt (used for docking) is the minimised one.
"""
import itertools
import os
import subprocess

BASES = ["d", "m", "h", "u"]
OUT = os.environ.get("OUT", "../trimers")
os.makedirs(OUT, exist_ok=True)

MIN_IN = """Trimer minimisation (GB)
 &cntrl
  imin=1, maxcyc=2000, ncyc=1000,
  cut=999.0, ntb=0, igb=5,
 /
"""

for combo in itertools.product(BASES, repeat=3):
    codon = "".join(combo).upper()          # e.g. "DMH"
    leap = f"""source leaprc.gaff2
loadamberparams D-capped_monomer.frcmod
loadamberparams M-capped_monomer.frcmod
loadamberparams H-capped_monomer.frcmod
loadamberparams U-capped_monomer.frcmod
loadoff D.lib
loadoff M.lib
loadoff H.lib
loadoff U.lib
trimer = sequence {{ {combo[0]} {combo[1]} {combo[2]} }}
saveamberparm trimer {OUT}/{codon}.prmtop {OUT}/{codon}.inpcrd
quit
"""
    leap_path = f"{OUT}/{codon}.leap"
    with open(leap_path, "w") as f:
        f.write(leap)
    subprocess.run(["tleap", "-f", leap_path], check=True)

    min_in = f"{OUT}/{codon}_min.in"
    with open(min_in, "w") as f:
        f.write(MIN_IN)
    subprocess.run(["sander", "-O", "-i", min_in, "-o", f"{OUT}/{codon}_min.out",
                    "-p", f"{OUT}/{codon}.prmtop", "-c", f"{OUT}/{codon}.inpcrd",
                    "-r", f"{OUT}/{codon}_min.rst"], check=True)
    with open(f"{OUT}/{codon}_min.pdb", "w") as f:
        subprocess.run(["ambpdb", "-p", f"{OUT}/{codon}.prmtop", "-c", f"{OUT}/{codon}_min.rst"],
                       stdout=f, check=True)
    print(f"{codon} done")
