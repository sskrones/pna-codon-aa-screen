"""Independent re-run of one codon–amino-acid pair with OpenMM, mirroring the
sander protocol in the notebook (igb=5/OBC2, mbondi radii from prmtop, no SA term,
no cutoff, SHAKE on H bonds, Langevin gamma=5, 2 fs):
  restrained min (k=5) -> 100 ps restrained heating 0->300 K (k=2) -> 150 ps free MD,
  frames every 3 ps (50 frames) -> single-trajectory dG = E_cpx - E_rec - E_lig.
Docking is replaced by a rigid-body random search on the GB interaction energy
(so the check needs only OpenMM + ParmEd); the question is only whether the ligand stays bound.
This is a check of the protocol, not a replacement for the screen: the starting pose
differs from Vina's, and OpenMM's OBC2 is the counterpart of sander igb=5.

Usage: python validation/rerun_pair_openmm.py CODON aa [seed] [--receptors data/receptors] [--ligands data/ligands]
       > validation/results/seed_CODON_aa_SEED.json
"""
import argparse, numpy as np, parmed as pmd, tempfile, os, json
import openmm as mm, openmm.app as app, openmm.unit as u

ap = argparse.ArgumentParser()
ap.add_argument("codon"); ap.add_argument("aa"); ap.add_argument("seed", nargs="?", type=int, default=1)
ap.add_argument("--receptors", default="data/receptors"); ap.add_argument("--ligands", default="data/ligands")
args = ap.parse_args()
codon, aa, seed = args.codon.upper(), args.aa.lower(), args.seed
rng = np.random.default_rng(seed)
tmp = tempfile.mkdtemp()

rec = pmd.load_file(f"{args.receptors}/{codon}.prmtop", f"{args.receptors}/{codon}.inpcrd")
lig = pmd.load_file(f"{args.ligands}/{aa}.prmtop", f"{args.ligands}/{aa}.inpcrd")
cpx = rec + lig
cpx.save(f"{tmp}/c.prmtop");
nrec = len(rec.atoms)

def system(prm):
    p = app.AmberPrmtopFile(prm)
    s = p.createSystem(nonbondedMethod=app.NoCutoff, implicitSolvent=app.OBC2,
                       gbsaModel=None, constraints=app.HBonds)
    return p, s

pc, sc = system(f"{tmp}/c.prmtop")
rec.save(f"{tmp}/r.prmtop"); lig.save(f"{tmp}/l.prmtop")
_, sr = system(f"{tmp}/r.prmtop"); _, sl = system(f"{tmp}/l.prmtop")
plat = mm.Platform.getPlatformByName("CPU")
def ctx(s): return mm.Context(s, mm.VerletIntegrator(1e-3), plat)
cc, cr, cl = ctx(sc), ctx(sr), ctx(sl)
KCAL = u.kilocalorie_per_mole
def E(c, x): c.setPositions(x * u.angstrom); return c.getState(getEnergy=True).getPotentialEnergy().value_in_unit(KCAL)
def dG(x): return E(cc, x) - E(cr, x[:nrec]) - E(cl, x[nrec:])

xr = np.array(rec.coordinates); xl0 = np.array(lig.coordinates); xl0 -= xl0.mean(0)
# ---- rigid-body search for a bound starting pose ----
def randrot():
    q = rng.normal(size=4); q /= np.linalg.norm(q); a, b, c_, d = q
    return np.array([[a*a+b*b-c_*c_-d*d, 2*(b*c_-a*d), 2*(b*d+a*c_)],
                     [2*(b*c_+a*d), a*a-b*b+c_*c_-d*d, 2*(c_*d-a*b)],
                     [2*(b*d-a*c_), 2*(c_*d+a*b), a*a-b*b-c_*c_+d*d]])
lo, hi = xr.min(0) - 2, xr.max(0) + 2
best = (1e9, None)
for _ in range(4000):
    xl = xl0 @ randrot().T + rng.uniform(lo, hi)
    dmin = np.min(np.linalg.norm(xr[:, None] - xl[None], axis=2))
    if dmin < 1.6 or dmin > 4.0: continue
    e = dG(np.vstack([xr, xl]))
    if e < best[0]: best = (e, xl)
x0 = np.vstack([xr, best[1]])
dg_pose = dG(x0)

# ---- MD protocol ----
heavy = [a.index for a in pc.topology.atoms() if a.element.symbol != "H"]
def add_restraint(s, ref, k):
    f = mm.CustomExternalForce("k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    f.addGlobalParameter("k", k * KCAL / u.angstrom**2)
    for p in ("x0", "y0", "z0"): f.addPerParticleParameter(p)
    for i in heavy: f.addParticle(i, (ref[i] / 10.0).tolist())
    return s.addForce(f)
import copy
s_min = copy.deepcopy(sc); add_restraint(s_min, x0, 5.0)
integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 5 / u.picosecond, 2 * u.femtosecond)
integ.setRandomNumberSeed(seed)
sim = app.Simulation(pc.topology, s_min, integ, plat)
sim.context.setPositions(x0 * u.angstrom); sim.minimizeEnergy(maxIterations=500)
xm = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.angstrom)
dg_min = dG(xm)
# heating with k=2 restraint to minimized structure
s_eq = copy.deepcopy(sc); add_restraint(s_eq, xm, 2.0)
integ = mm.LangevinMiddleIntegrator(0 * u.kelvin, 5 / u.picosecond, 2 * u.femtosecond); integ.setRandomNumberSeed(seed)
sim = app.Simulation(pc.topology, s_eq, integ, plat); sim.context.setPositions(xm * u.angstrom)
for i in range(50):   # 50 x 1000 steps = 100 ps, temperature ramp 0 -> 300 K
    integ.setTemperature(300 * (i + 1) / 50 * u.kelvin); sim.step(1000)
st = sim.context.getState(getPositions=True, getVelocities=True)
integ = mm.LangevinMiddleIntegrator(300 * u.kelvin, 5 / u.picosecond, 2 * u.femtosecond); integ.setRandomNumberSeed(seed + 7)
sim = app.Simulation(pc.topology, sc, integ, plat)
sim.context.setPositions(st.getPositions()); sim.context.setVelocities(st.getVelocities())
rows = []
for f in range(50):   # 150 ps, frame every 1500 steps (3 ps)
    sim.step(1500)
    x = sim.context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(u.angstrom)
    hr = [i for i in heavy if i < nrec]; hl = [i for i in heavy if i >= nrec]
    dmin = float(np.min(np.linalg.norm(x[hr][:, None] - x[hl][None], axis=2)))
    com = float(np.linalg.norm(x[:nrec].mean(0) - x[nrec:].mean(0)))
    rows.append((dG(x), dmin, com))
r = np.array(rows)
out = dict(codon=codon, aa=aa.upper(), seed=seed, dg_docked_pose=round(dg_pose, 2), dg_after_min=round(dg_min, 2),
           dg_md_mean=round(r[:, 0].mean(), 3), dg_md_std=round(r[:, 0].std(), 3),
           frac_frames_contact_lt4A=round(float((r[:, 1] < 4.0).mean()), 2),
           mindist_first=round(r[0, 1], 1), mindist_last=round(r[-1, 1], 1), com_last=round(r[-1, 2], 1),
           mindist_series=[round(v, 1) for v in r[:, 1]], dg_series=[round(v, 2) for v in r[:, 0]])
print(json.dumps(out))
