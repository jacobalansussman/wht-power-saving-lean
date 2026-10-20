"""fullrow.py -- e3-bpexact: UPPER bound version of the per-irrep distance-4 test: the top's row block w is allowed to be
ANY ds x 30ds matrix (not only Wsp = {E_d1 R_sigma(A)}).  r_a(full) >= r_a(Wsp); the F2 criterion r_a + m_a(pi) >= ds is
necessary with r_a(full) too, so a failure here is a failure for every top content.  Generic (designed) lift of run1.
argv: seed psi nBm sigma_names"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
t0 = time.time()
seed, psi, nBm, names = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4].split(',')
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
cons = {'all30': [d for d in range(30) if d != D.d0 and d not in Fr], 'd4': list(D.Sx)}
for name in names:
    sig, ds = (sig_nat, 16) if name == 'V' else make_sig(name)
    Rk = sum(y * D.R(k, sig, ds) for y, k in zip(ys, ks))
    Rp = [D.R(p, sig, ds) for p in L.p]
    col0 = np.arange(D.d0 * ds, (D.d0 + 1) * ds)
    for mode, cs in cons.items():
        colx = np.concatenate([np.arange(d * ds, (d + 1) * ds) for d in cs])
        out = []
        for a in range(7):
            Ma = (sum(L.chars[a][j] * Rp[j] for j in range(7)) / 7) @ Rk
            C = Ma[:, colx]; U, s, Vh = np.linalg.svd(C, full_matrices=True); r = int((s > 1e-9 * max(1, s[0])).sum())
            Nul = U[:, r:].conj().T; B0 = Nul @ Ma[:, col0]
            sv = np.linalg.svd(B0, compute_uv=False) if B0.size else np.zeros(1)
            ra = int((sv > 1e-8 * np.abs(Ma).max()).sum())
            ma = int(round(np.trace(sum(L.chars[a][j] * sig(L.p[j]) for j in range(7)) / 7).real))
            out.append((a, ra, ma))
        need = sum(ds - ra for a, ra, ma in out)
        print(f'[fullrow] {name} (dim {ds}) {mode}: (a: r_a(full) / m_a(p)) ' + ', '.join(f'{a}: {ra}/{ma}' for a, ra, ma in out) +
              f'; total deficit {need} vs ds {ds} -> {"EXCLUDED for every F2 twist" if need > ds else "not excluded by this count"}'
              f'  {time.time() - t0:.0f}s', flush=True)
print('done %.0fs' % (time.time() - t0))
