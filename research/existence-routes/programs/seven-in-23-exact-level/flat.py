"""flat.py -- e3-bpexact: sigma = V with the PAULI-CORRECTED ("flat") lift, and what F2 can add.

(1) Check R_V(k) = Delta (Pi_k (x) k) Delta^-1, Delta = blockdiag(g_d) (V extends to the whole group).
(2) Flat lift: replace each content letter k_b by sum_pi c_pi (pi k_b) with sum_pi c_pi pi = k_b^-1 (Paulis span M16);
    then R_V(kappa*) = Delta (Y (x) I) Delta^-1 with Y = R_triv(kappa) (the designed relaxation element), and the content is
    still letters of B_m (Paulis times the old letters) with complex coefficients.  Checked on a few terms.
(3) Top-alone kernel-block rank r_a per character a, for the generic lift and for the flat lift.
(4) F2 criterion [argued in NOTE section 4]: F2 absorbed with twist homomorphism pi_g (order 7) contributes in character a
    the block Y2 E_a(pi) phi (E_a = (1/7) sum_j a(j) V(pi)^j, rank m_a(pi)); then character a can be served iff
    r_a >= 16 - m_a(pi).  Printed for pi = p (F2 twisted like F1).
argv: seed psi nBm"""
import sys, os, time, random, itertools
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
t0 = time.time()
seed, psi, nBm = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
L = Lift(seed=seed, psi=psi)
rng = random.Random(seed + 1); nrng = np.random.default_rng(seed + 2)
ys = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), f'ys_seed{seed}_psi{psi}.npy'))
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < nBm:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
mi = minv(L.m); ks = [mm(L.m, b, mi) for b in Bm]
D = D4(L); P = L.P
Delta = np.zeros((480, 480), dtype=complex)
for d in range(30): Delta[16 * d:16 * d + 16, 16 * d:16 * d + 16] = D.g[d]
Di = np.linalg.inv(Delta)
def Pi(k):
    M = np.zeros((30, 30))
    for d, (dd, beta) in enumerate(D.act(k)): M[d, dd] = 1
    return M
for k in ks[:5] + L.p[1:3]:
    lhs = D.R(k, sig_nat, 16); rhs = Delta @ np.kron(Pi(k), k) @ Di
    assert np.allclose(lhs, rhs), 'R_V != Delta (Pi (x) V) Delta^-1'
print('[flat] (1) R_V(k) = Delta (Pi_k (x) k) Delta^-1 checked on 7 letters', flush=True)
# Paulis
def pauli(a, b):
    M = np.zeros((16, 16), dtype=complex)
    for y in range(16): M[y ^ a, y] = (-1) ** bin(b & y).count('1')
    return M
PA = [pauli(a, b) for a in range(16) for b in range(16)]
PAv = np.array([M.reshape(-1) for M in PA]).T          # 256 x 256, columns vec(pi)
PAinv = np.linalg.inv(PAv)
Y = sum(y * Pi(k) for y, k in zip(ys, ks))
Kflat = Delta @ np.kron(Y, np.eye(16)) @ Di
# check (2) on actual terms for two letters: sum_pi c_pi R_V(pi k) = Delta (Pi_k (x) I) Delta^-1
for k in ks[:2]:
    c = PAinv @ np.linalg.inv(k).reshape(-1)
    acc = sum(c[i] * D.R(mm(PA[i], k), sig_nat, 16) for i in range(256) if abs(c[i]) > 1e-12)
    assert np.allclose(acc, Delta @ np.kron(Pi(k), np.eye(16)) @ Di), 'flat lift identity fails'
    nz = int((np.abs(c) > 1e-12).sum())
print(f'[flat] (2) flat-lift identity checked on actual Pauli-multiplied terms (each letter -> {nz} Pauli terms)', flush=True)
Kgen = sum(y * D.R(k, sig_nat, 16) for y, k in zip(ys, ks))
RVp = [D.R(p, sig_nat, 16) for p in L.p]
Fr_fixed = [d for d in range(30) if all(Pi(k)[d, d] == 1 for k in ks[:200])]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(Fr_fixed))
modes = {'all30': [d for d in range(30) if d != D.d0 and d not in Fr], 'd4': D.Sx}
col0 = np.arange(16 * D.d0, 16 * D.d0 + 16)
def kernel_rank(M, cons):
    colx = np.concatenate([np.arange(16 * d, 16 * d + 16) for d in cons])
    C = M[:, colx]
    U, s, Vh = np.linalg.svd(C, full_matrices=True); r = int((s > 1e-9 * max(1, s[0])).sum())
    Nul = U[:, r:].conj().T
    B0 = Nul @ M[:, col0]
    if not B0.size: return 0, 0.
    sv = np.linalg.svd(B0, compute_uv=False)
    scale = np.abs(M).max()                       # absolute scale of the problem (Nul rows are orthonormal)
    return int((sv > 1e-8 * scale).sum()), float(sv[0] / scale)
# eigen-multiplicities of V(p): E_a(p) = (1/7) sum_j a(j) p^j
m_p = {a: int(round(np.trace(sum(L.chars[a][j] * L.p[j] for j in range(7)) / 7).real)) for a in range(7)}
print('[flat] multiplicities m_a(p) = rank E_a(V(p)):', m_p, ' sum', sum(m_p.values()), flush=True)
for mode, cons in modes.items():
    for name, K in (('generic lift', Kgen), ('flat lift', Kflat)):
        out = []
        for a in range(7):
            Ma = (sum(L.chars[a][j] * RVp[j] for j in range(7)) / 7) @ K
            r, top = kernel_rank(Ma, cons)
            out.append((a, r, 16 - m_p[a], r >= 16 - m_p[a], top))
        print(f'[flat] {mode:5s} {name:12s}: top-alone kernel rank r_a / needed with F2 twisted by p (16 - m_a): ' +
              ', '.join(f'a={a}: {r}/{need}{"" if ok else " FAIL"} (top sv {top:.0e})' for a, r, need, ok, top in out), flush=True)
print('done %.0fs' % (time.time() - t0))
