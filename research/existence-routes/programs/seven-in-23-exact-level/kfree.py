"""kfree.py -- e3-bpexact: is the sigma = V failure of the lifted survivor a property of the LIFT or of the FAMILY?
kappa is made free: kappa = sum_k y_k k over N letters k of B'_m = B cap m B m^-1 (products of the sampled ones).
For a fixed generic row block w' = W0 R_sigma(e_a) (rows in the a-isotypic part, W0 random), the condition
      w' R_sigma(kappa) Pi_(constrained slots) = 0,   w' R_sigma(kappa) Pi_d0 = I
is LINEAR in y.  Solvable for a generic w' => sigma does not exclude the family (for that character).  Also solved
jointly with the trivial-irrep (frame-sum) conditions for all six characters on the SAME y (generic gamma per character).
Slots: 'all30' = every slot except d0 and the 4 reachable ones; 'd4' = the 15 other distance-4 slots.
argv: seed psi N mode(all30|d4) a"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
t0 = time.time()
seed, psi, N, mode, a_ch = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5])
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 7); nrng = np.random.default_rng(seed + 8)
base = [b for (b, s) in L.Bm]
Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < N:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
mi = minv(L.m); ks = [mm(L.m, b, mi) for b in Bm]
D = D4(L); P = L.P
acts = [D.act(k) for k in ks]
print(f'[kfree] {len(ks)} letters of B\'_m, act computed  {time.time() - t0:.0f}s', flush=True)
Rt0 = [D.R(k, sig_triv, 1).real for k in ks[:400]]
fixed = [d for d in range(30) if all(R[d, d] == 1 for R in Rt0)]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(fixed))
cons = [d for d in range(30) if d != D.d0 and d not in Fr] if mode == 'all30' else D.Sx
def rows_for(sig, ds, wprime, target):
    """columns (one per letter) of vec(w' R_sigma(k) restricted to cons + d0), and the target vector"""
    blocks = cons + [D.d0]
    pos = {d: i for i, d in enumerate(blocks)}
    nr = wprime.shape[0]
    A = np.zeros((nr * ds * len(blocks), len(ks)), dtype=complex)
    wp = wprime.reshape(nr, 30, ds)
    for j, act in enumerate(acts):
        out = np.zeros((nr, len(blocks), ds), dtype=complex)
        for d, (dd, beta) in enumerate(act):
            if dd in pos: out[:, pos[dd], :] += wp[:, d, :] @ sig(beta)
        A[:, j] = out.reshape(-1)
    return A, target.reshape(-1)
Rp_cache = {}
def Re(sig, ds, a):
    Rp = [D.R(p, sig, ds) for p in L.p]
    return sum(L.chars[a][j] * Rp[j] for j in range(7)) / 7
blocks = []; rhs = []; labels = []
# sigma = V, character a_ch: generic w'
ds = 16
W0 = nrng.standard_normal((16, 480)) + 1j * nrng.standard_normal((16, 480))
wV = W0 @ Re(sig_nat, 16, a_ch)
tgt = np.zeros((16, len(cons) + 1, 16), dtype=complex); tgt[:, -1, :] = np.eye(16)
A, b = rows_for(sig_nat, 16, wV, tgt); blocks.append(A); rhs.append(b); labels.append(('V', a_ch, A.shape[0]))
nV = A.shape[0]
xV, *_ = np.linalg.lstsq(A, b, rcond=None)
rV = np.linalg.norm(A @ xV - b) / np.linalg.norm(b)
print(f'[kfree] sigma = V alone, a = {a_ch}, mode {mode}: {A.shape[0]} equations, {A.shape[1]} unknowns; relative residual {rV:.2e}  {time.time() - t0:.0f}s', flush=True)
# trivial irrep, all six characters, generic gamma per character (the relaxation's W_a), free slots unconstrained
for a in [1, 2, 3, 4, 5, 6]:
    g = nrng.standard_normal((1, 30)) + 1j * nrng.standard_normal((1, 30))
    wT = g @ Re(sig_triv, 1, a)
    t1 = np.zeros((1, len(cons) + 1, 1), dtype=complex); t1[0, -1, 0] = 1
    A1, b1 = rows_for(sig_triv, 1, wT, t1); blocks.append(A1); rhs.append(b1); labels.append(('triv', a, A1.shape[0]))
At = np.vstack(blocks[1:]); bt = np.concatenate(rhs[1:])
xt, *_ = np.linalg.lstsq(At, bt, rcond=None)
print(f'[kfree] trivial irrep alone, 6 characters: relative residual {np.linalg.norm(At @ xt - bt) / np.linalg.norm(bt):.2e}', flush=True)
Aall = np.vstack(blocks); ball = np.concatenate(rhs)
x, *_ = np.linalg.lstsq(Aall, ball, rcond=None)
res = Aall @ x - ball
print(f'[kfree] JOINT (V at a = {a_ch} + trivial at all 6): {Aall.shape[0]} equations, {Aall.shape[1]} unknowns; relative residual '
      f'{np.linalg.norm(res) / np.linalg.norm(ball):.2e}; V part {np.linalg.norm(res[:nV]) / np.linalg.norm(ball[:nV]):.2e}, '
      f'trivial part {np.linalg.norm(res[nV:]) / np.linalg.norm(ball[nV:]):.2e}  {time.time() - t0:.0f}s', flush=True)
s = np.linalg.svd(Aall, compute_uv=False)
print(f'[kfree] joint matrix rank {int((s > 1e-9 * s[0]).sum())} of {Aall.shape[0]} rows (singular values: top {s[0]:.2e}, '
      f'smallest {s[-1]:.2e})', flush=True)
np.save(f'kfree_y_{mode}_a{a_ch}.npy', x)
print('done %.0fs' % (time.time() - t0))
