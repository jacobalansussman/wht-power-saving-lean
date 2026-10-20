"""run1.py -- e3-bpexact: lift the survivor, re-design its group-algebra content Y on the lift, check the frame-sum
relaxation with w5-blueprint's relax.Eng.decide (independent machinery, top content free), then the EXACT distance-4
test of d4test.py in the irreps sigma = trivial (calibration), V (natural, 16-dim), conj(V).
argv: seed psi nBm"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
from relax import Eng
t0 = time.time()
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 21
psi = int(sys.argv[2]) if len(sys.argv) > 2 else 1
nBm = int(sys.argv[3]) if len(sys.argv) > 3 else 400
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 1); nrng = np.random.default_rng(seed + 2)
# more letters of B_m: random products of the sampled ones (B_m is a group)
base = [b for (b, s) in L.Bm]
Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < nBm:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
mi = minv(L.m)
ks = [mm(L.m, b, mi) for b in Bm]
assert all(is_letter(k) for k in ks)
D = D4(L)
# trivial-irrep permutation matrices of k_b and p_j
Rt = [D.R(k, sig_triv, 1).real for k in ks]
Rp = [D.R(p, sig_triv, 1).real for p in L.p]
span = np.linalg.matrix_rank(np.array([R.reshape(-1) for R in Rt]))
print(f'[run1] |B_m letters| {len(Bm)}, span of their 30x30 slot permutations {span} (A has dim 184)  {time.time() - t0:.0f}s', flush=True)
# free slots (relaxation's F): common neighbours of Lam_1 and Lam_F, plus slots fixed by every k_b
P = L.P
fixed = [d for d in range(30) if all(R[d, d] == 1 for R in Rt)]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(fixed))
print('[run1] fixed slots', fixed, ' free slots', Fr, ' d0', D.d0, ' S', D.S, flush=True)
# design y: for each nontrivial effective character a, W_a = gamma R(P_a) (gamma random), W_a Y = e_d0 + f_a (f_a on Fr)
own = 0
req = [a for a in range(7) if a != own]
gamma = nrng.standard_normal(30) + 1j * nrng.standard_normal(30)
RPa = {a: sum(L.chars[a][j] * Rp[j] for j in range(7)) for a in range(7)}
nb = len(ks); rowsM = []; rhs = []
for ia, a in enumerate(req):
    Wa = gamma @ RPa[a]
    blk = np.zeros((30, nb + len(Fr) * len(req)), dtype=complex)
    for i, R in enumerate(Rt): blk[:, i] = Wa @ R
    for j, d in enumerate(Fr): blk[d, nb + len(Fr) * ia + j] = -1
    e = np.zeros(30, dtype=complex); e[D.d0] = 1
    rowsM.append(blk); rhs.append(e)
Mbig = np.vstack(rowsM); bb = np.concatenate(rhs)
sol, *_ = np.linalg.lstsq(Mbig, bb, rcond=None)
ys = sol[:nb]
# make the lift generic: add a random null-space vector of the y-part (keeps the relaxation solution)
_, s, Vh = np.linalg.svd(Mbig); r = int((s > 1e-9 * s[0]).sum()); Nn = Vh[r:].conj().T
if Nn.shape[1]:
    z = Nn @ (nrng.standard_normal(Nn.shape[1]) + 1j * nrng.standard_normal(Nn.shape[1]))
    sol2 = sol + z * (np.linalg.norm(sol) / max(1e-12, np.linalg.norm(z)))
    ys = sol2[:nb]
print('[run1] reduced problem residual %.2e ; null dim %d' % (np.linalg.norm(Mbig @ sol - bb), Nn.shape[1]), flush=True)
# the lifted F1
F1 = L.F1(ys, Bm)
print('[run1] F1 terms', len(F1), ' absorption defect %.1e' % L.check_absorbed(F1), f' {time.time() - t0:.0f}s', flush=True)
# independent check: relax.Eng.decide on the lifted terms (cells from the actual matrices)
E = Eng('Z7', P=L.P, seed=seed)
def conv(u):
    out = {}
    for (l, kk), (c, w) in u.items():
        k2 = (l, inv(sigma(w))); out[k2] = out.get(k2, 0) + c
    return out
O1 = {(0, inv(sigma(L.m1))): 1.0 + 0j}; O2 = {(0, inv(sigma(L.g2))): 1.0 + 0j}
F1c = conv(F1)
ok, info = E.decide([E.inp(0), O1, O2], [F1c], E=[psi])
print('[run1] FRAME-SUM RELAXATION of the lifted skeleton (relax.Eng.decide, top free): passes =', ok,
      ' mu_ok', info.get('mu_ok'), ' sv_gap', info.get('sv_gap'), f' {time.time() - t0:.0f}s', flush=True)
okI, infoI = E.decide([E.inp(0), O1, O2], [conv(L.F1([1.0], [I16]))], E=[psi])
print('[run1] control Y = 1: passes =', okI, flush=True)
# exact distance-4 test
for name, sig, ds, full in (('trivial', sig_triv, 1, True), ('natural V', sig_nat, 16, True), ('conj V', sig_conj, 16, True)):
    res = D.test(sig, ds, ys, ks, req, None, nrng, full_rows=full)
    print(f'[d4] sigma = {name} (dim {ds}): null dim {res["dim_null"]}; rank of the kernel block per effective character a '
          f'(need {ds}): ' + ', '.join(f'a={a}:{res[a]}' for a in range(7)), f'  {time.time() - t0:.0f}s', flush=True)
# with no requirement at all (top serving one character at a time)
for name, sig, ds in (('natural V', sig_nat, 16),):
    for a in range(7):
        res = D.test(sig, ds, ys, ks, [a], None, nrng, full_rows=True)
        print(f'[d4] sigma = {name}, character a={a} ALONE: null dim {res["dim_null"]}, kernel-block rank {res[a]} (need {ds})', flush=True)
np.save(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'ys_seed{seed}_psi{psi}.npy'), ys)
print('done %.0fs' % (time.time() - t0))
