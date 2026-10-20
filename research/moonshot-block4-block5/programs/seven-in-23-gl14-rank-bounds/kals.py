"""kals.py -- e3-bpexact: alternating least squares over (top row block w, kappa) in one irrep sigma and one character a,
looking for SPECIAL (w, kappa) that reach a kernel block of rank r_target with no junk (kmax only tried a generic w).
Unknowns: w (r_target x 30 d) arbitrary rows (upper-bound version), kappa = sum_k y_k k over N letters of B'_m.
Equations:  w R(e_a) R(kappa) Pi_x = 0,   w R(e_a) R(kappa) Pi_d0 = Z  (Z: r_target x d, full row rank; for r_target < d its row space is updated to the current kernel block (polar factor)).
w-step and y-step are linear least squares.  Several random starts.  Residual -> 0 means such (w, kappa) exist.
argv: seed psi N sigma a r_target starts iters mode(d4|all30)"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
t0 = time.time()
seed, psi, N, sname, a, rt, starts, iters, mode = (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4],
    int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), sys.argv[9])
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 7); nrng = np.random.default_rng(seed + 10)
sig, ds = (sig_nat, 16) if sname == 'V' else make_sig(sname)
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm); seen_s = set()
mi = minv(L.m); D = D4(L); P = L.P
ks = []
while len(ks) < N:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) in seen: continue
    seen.add(key(w)); Bm.append(w); k = mm(L.m, w, mi)
    if sname in ('gl14', 'gl15'):
        sk = sigma(k)
        if sk in seen_s: continue
        seen_s.add(sk)
    ks.append(k)
# R_sigma(k) restricted: we need M_k = R(e_a) R(k) [:, cols]; store R(k)[:, cols] blocks via act
acts = [D.act(k) for k in ks]
Rp = [D.R(p, sig, ds) for p in L.p]
Ea = sum(L.chars[a][j] * Rp[j] for j in range(7)) / 7
Fr_fixed = [d for d in range(30) if all(act[d][0] == d for act in acts[:300])]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(Fr_fixed))
cons = [d for d in range(30) if d != D.d0 and d not in Fr] if mode == 'all30' else D.Sx
blocks = cons + [D.d0]; pos = {d: i for i, d in enumerate(blocks)}
sigb = [[sig(beta) for (dd, beta) in act] for act in acts]
print(f'[kals] sigma {sname} (dim {ds}), a = {a}, target rank {rt}, {len(ks)} letters, mode {mode}; setup {time.time() - t0:.0f}s', flush=True)
def apply_all(wp):
    """wp: (r, 30 ds) -> array (N, r * len(blocks) * ds): vec(wp R(k)|cols) for every letter"""
    r = wp.shape[0]; W3 = wp.reshape(r, 30, ds)
    out = np.zeros((len(ks), r, len(blocks), ds), dtype=complex)
    for j, act in enumerate(acts):
        for d, (dd, beta) in enumerate(act):
            if dd in pos: out[j, :, pos[dd], :] += W3[:, d, :] @ sigb[j][d]
    return out.reshape(len(ks), -1)
def RK(y):
    M = np.zeros((30 * ds, len(blocks) * ds), dtype=complex)
    for j, act in enumerate(acts):
        if abs(y[j]) < 1e-15: continue
        for d, (dd, beta) in enumerate(act):
            if dd in pos: M[d * ds:(d + 1) * ds, pos[dd] * ds:(pos[dd] + 1) * ds] += y[j] * sigb[j][d]
    return M
best = 9
for st in range(starts):
    Z = nrng.standard_normal((rt, ds)) + 1j * nrng.standard_normal((rt, ds))
    T = np.zeros((rt, len(blocks), ds), dtype=complex); T[:, -1, :] = Z; T = T.reshape(rt, -1)
    y = nrng.standard_normal(len(ks)) + 1j * nrng.standard_normal(len(ks))
    for it in range(iters):
        M = Ea @ RK(y)                                   # (30ds, cols)
        w, *_ = np.linalg.lstsq(M.T, T.T, rcond=None); w = w.T     # w M = T
        wp = w @ Ea
        A = apply_all(wp)                                # (N, rt*cols)
        y, *_ = np.linalg.lstsq(A.T, T.reshape(-1), rcond=None)
        res = np.linalg.norm(A.T @ y - T.reshape(-1)) / np.linalg.norm(T)
        if rt < ds:      # the kernel block's ROW SPACE is free: move Z to the polar factor of the current kernel block
            C0 = (A.T @ y).reshape(rt, len(blocks), ds)[:, -1, :]
            U_, s_, Vh_ = np.linalg.svd(C0, full_matrices=False)
            Z = (U_ @ Vh_) * np.sqrt(ds)
            T = np.zeros((rt, len(blocks), ds), dtype=complex); T[:, -1, :] = Z; T = T.reshape(rt, -1)
        if it % 5 == 0 or it == iters - 1:
            print(f'   start {st} iter {it}: relative residual {res:.3e}  {time.time() - t0:.0f}s', flush=True)
    best = min(best, res)
print(f'[kals] {sname} a={a} rank-{rt} target: best relative residual {best:.3e} (0 = reachable with no junk)')
print('done %.0fs' % (time.time() - t0))
