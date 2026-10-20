"""run1b.py -- e3-bpexact: the exact test on ALL 30 neighbour slots of F (not only the 16 at distance 4): every slot that
no other generator can reach must carry no junk.  Reachable (free) slots, as in the relaxation: the neighbours of O1's cell
lam1 (the top reading O1 members) and the two G-fixed slots (t1 itself, read directly by the readout, and the third point
of the line t1-F, reached by the top reading O2).  F2 aside.  sigma = trivial, V, conj(V).  argv: seed psi nBm"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
t0 = time.time()
seed = int(sys.argv[1]); psi = int(sys.argv[2]); nBm = int(sys.argv[3])
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 1); nrng = np.random.default_rng(seed + 2)
ys = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'ys_seed{seed}_psi{psi}.npy'))
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < nBm:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
mi = minv(L.m); ks = [mm(L.m, b, mi) for b in Bm]
D = D4(L); P = L.P
Rt = [D.R(k, sig_triv, 1).real for k in ks]
fixed = [d for d in range(30) if all(R[d, d] == 1 for R in Rt)]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(fixed))
dist = [int(P.dZ[D.slot[d]]) for d in range(30)]
print('[run1b] slot input distances', dist, ' free', Fr, ' d0', D.d0, flush=True)
req = [1, 2, 3, 4, 5, 6]
saveSx = D.Sx
D.Sx = [d for d in range(30) if d != D.d0 and d not in Fr]
print('[run1b] constrained slots', len(D.Sx), flush=True)
for name, sig, ds in (('trivial', sig_triv, 1), ('natural V', sig_nat, 16), ('conj V', sig_conj, 16)):
    res = D.test(sig, ds, ys, ks, req, None, nrng, full_rows=True)
    print(f'[all30] sigma = {name}: null dim {res["dim_null"]} (junk check {res["junk"]:.1e}); kernel-block rank (need {ds}) ' +
          ', '.join(f'a={a}:{res[a]}' for a in range(7)), f' {time.time() - t0:.0f}s', flush=True)
# a = 0 (F1's own character) alone, all 30 slots: F2 must serve it at the frame level
for name, sig, ds in (('trivial', sig_triv, 1), ('natural V', sig_nat, 16)):
    res = D.test(sig, ds, ys, ks, [0], None, nrng, full_rows=True)
    print(f'[all30] sigma = {name}, a = 0 alone: null dim {res["dim_null"]}, kernel-block rank {res[0]} (need {ds})', flush=True)
# control: a generic kappa (random y on the same letters, NOT designed): must fail in sigma = trivial
yr = nrng.standard_normal(len(ks)) + 1j * nrng.standard_normal(len(ks))
for name, sig, ds in (('trivial', sig_triv, 1), ('natural V', sig_nat, 16)):
    res = D.test(sig, ds, yr, ks, req, None, nrng, full_rows=True)
    print(f'[all30] CONTROL random y: sigma = {name}: null dim {res["dim_null"]}; kernel ranks ' +
          ', '.join(f'a={a}:{res[a]}' for a in req), flush=True)
    D.Sx = saveSx
    res = D.test(sig, ds, yr, ks, req, None, nrng, full_rows=True)
    print(f'[d4]    CONTROL random y: sigma = {name}: null dim {res["dim_null"]}; kernel ranks ' +
          ', '.join(f'a={a}:{res[a]}' for a in req), flush=True)
    D.Sx = [d for d in range(30) if d != D.d0 and d not in Fr]
print('done %.0fs' % (time.time() - t0))
