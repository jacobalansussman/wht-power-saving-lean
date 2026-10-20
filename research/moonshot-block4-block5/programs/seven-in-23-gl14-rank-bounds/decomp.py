"""decomp.py -- w8-bpgl14: the B'_m-decomposition of the slot module M_sigma (sigma = one irrep of the letters, or a
block-diagonal SUM of several, so that constituents shared between them are merged into one type: the coupling that
kappa is ONE element of C[B'_m] for every sigma).  Same objects as e3-bpexact/crit.py, faster commutant:
  X_(g.d, g.e) = sigma(beta(d,g))^-1 X_(d,e) sigma(beta(e,g))   (commutant of R_sigma(B'_m); orbitals of slot pairs)
solved per orbital AND per (sigma_i, sigma_j) sub-block by successive null-space restriction (exact up to round-off).
Then: random Hermitian commutant element -> irreducible copies; a second one -> types and consistent (intertwined) bases.
Saved (npz): for each type tau: d_tau, m_tau, the copies' orthonormal bases U_k (n x d_tau, consistent), and
n_(tau,a) = rank of rho_tau(e_a) for the twist p of F1 and the second twist p2 of O O F F O (ooffo.py's q2), a = 0..6.
argv: seed psi sigmas(comma: gl14,trivial,aff16,sq16,gl15,V) ngens out.npz"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
t0 = time.time()
seed, psi, names, ngens, outf = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3].split(','), int(sys.argv[4]), sys.argv[5]
L = Lift(seed=seed, psi=psi)
_v16 = np.ones(16) / 4; _U16, _, _ = np.linalg.svd(np.eye(16) - np.outer(_v16, _v16)); _Q15 = _U16[:, :15].T
def sig_15p(b): return _Q15 @ sig_aff16(b) @ _Q15.T          # aff16 = trivial + 15p (15p irreducible: irrchk.py)
def get_sig(name):
    if name == 'trivial': return sig_triv, 1
    if name == '15p': return sig_15p, 15
    if name == 'V': return sig_nat, 16
    return make_sig(name)
sigl = [get_sig(nm) for nm in names]
dss = [ds for (_, ds) in sigl]; offs = np.cumsum([0] + dss); DS = int(offs[-1])
def sigJ(b):
    M = np.zeros((DS, DS), dtype=complex)
    for (s, ds), o in zip(sigl, offs[:-1]): M[o:o + ds, o:o + ds] = s(b)
    return M
n = 30 * DS
mi = minv(L.m); D = D4(L); P = L.P
# ---- second twist q2 exactly as in e3-bpexact/ooffo.py (seed+1 rng, 400 letters of B_m, first suitable)
rng = random.Random(seed + 1)
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < 400:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
Pi = lambda k: D.R(k, sig_triv, 1).real
def exact7(b):
    s = sigma(b); o = L.sp_order(s)
    if o <= 0 or o % 7: return None
    qq = snap(np.linalg.matrix_power(b, o // 7)); M7 = np.linalg.matrix_power(qq, 7); a2 = 0
    while not np.allclose(M7, I16):
        M7 = snap(M7 @ M7); a2 += 1
        if a2 > 10: return None
    e = (2 ** a2) * pow(2 ** a2, -1, 7)
    qf = snap(np.linalg.matrix_power(qq, e))
    return qf if np.allclose(np.linalg.matrix_power(qf, 7), I16) and not np.allclose(qf, I16) else None
q1 = L.q; pow1 = set(key(x) for x in L.qp); q2 = None
for b in Bm:
    qf = exact7(b)
    if qf is not None and key(qf) not in pow1 and not any(np.allclose(sigma(qf), sigma(x)) for x in L.qp):
        P1 = Pi(mm(L.m, q1, mi)); P2 = Pi(mm(L.m, qf, mi))
        if not any(np.allclose(np.linalg.matrix_power(P1, j), P2) for j in range(7)):
            q2 = qf; break
assert q2 is not None
q2p = [I16]
for j in range(1, 7): q2p.append(mm(q2p[-1], q2))
p2 = [mm(L.m, x, mi) for x in q2p]
# ---- generators of B'_m = m B_m m^-1 and test letters
rr = random.Random(seed + 99)
gens = [mm(L.m, b, mi) for b in base[:ngens]] + [L.p[1], p2[1]]
tests = [mm(L.m, rr.choice(Bm), mi) for _ in range(12)]
gact = [D.act(k) for k in gens]
gs = [[sigJ(beta) for (dd, beta) in act] for act in gact]
gsi = [[x.conj().T for x in row] for row in gs]          # sigma unitary
print(f'[decomp] {names} DS {DS} n {n}; {len(gens)} generators; setup {time.time() - t0:.0f}s', flush=True)
# ---- propagation over slot pairs: rep_of, trans (A, C) with X_pair = A X_rep C; each revisit gives a constraint
# X_rep = A X_rep C, applied at once to the current solution basis of every (sigma_i, sigma_j) sub-block
def restrict(B, A, C, i, j):
    si, sj = slice(offs[i], offs[i + 1]), slice(offs[j], offs[j + 1]); di, dj = dss[i], dss[j]
    if B.shape[0] == 0: return B
    Xs = B.reshape(-1, di, dj)
    Rm = (np.einsum('ab,kbc,cd->kad', A[si, si], Xs, C[sj, sj]) - Xs).reshape(B.shape[0], -1)
    if np.abs(Rm).max() < 1e-11: return B
    U_, s_, Vh_ = np.linalg.svd(Rm.T, full_matrices=True)
    rk = int((s_ > 1e-9).sum())
    B = Vh_[rk:].conj() @ B
    if B.shape[0] == 0: return B
    q, _ = np.linalg.qr(B.T); return q.T
rep_of = {}; trans = {}; sol = {}; ncons = 0
nb = len(dss)
for d in range(30):
    for e in range(30):
        if (d, e) in rep_of: continue
        r = (d, e)
        rep_of[r] = r; trans[r] = (np.eye(DS, dtype=complex), np.eye(DS, dtype=complex))
        Bs = {(i, j): np.eye(dss[i] * dss[j], dtype=complex) for i in range(nb) for j in range(nb)}
        fr = [r]
        while fr:
            nf = []
            for (x, y) in fr:
                A, C = trans[(x, y)]
                for gi, act in enumerate(gact):
                    xx, yy = act[x][0], act[y][0]
                    A2 = gsi[gi][x] @ A; C2 = C @ gs[gi][y]
                    if (xx, yy) not in rep_of:
                        rep_of[(xx, yy)] = r; trans[(xx, yy)] = (A2, C2); nf.append((xx, yy))
                    else:
                        A1, C1 = trans[(xx, yy)]
                        Ac, Cc = A1.conj().T @ A2, C2 @ C1.conj().T      # X = Ac X Cc  (A1, C1 unitary)
                        ncons += 1
                        for ij in Bs: Bs[ij] = restrict(Bs[ij], Ac, Cc, *ij)
            fr = nf
        sol[r] = Bs
reps = list(sol)
print(f'[decomp] {len(reps)} orbitals of slot pairs; {ncons} constraints  {time.time() - t0:.0f}s', flush=True)
basis = []          # (rep, i, j, X sub-block)
for r in reps:
    for (i, j), Bb in sol[r].items():
        for v in Bb: basis.append((r, i, j, v.reshape(dss[i], dss[j])))
print(f'[decomp] commutant dim {len(basis)}  {time.time() - t0:.0f}s', flush=True)
pairs_of = {r: [pr for pr, rr_ in rep_of.items() if rr_ == r] for r in reps}
def build(coefs):
    X = np.zeros((n, n), dtype=complex)
    for c, (r, i, j, Xr) in zip(coefs, basis):
        if abs(c) < 1e-15: continue
        Xf = np.zeros((DS, DS), dtype=complex); Xf[offs[i]:offs[i + 1], offs[j]:offs[j + 1]] = Xr
        for pr in pairs_of[r]:
            A, C = trans[pr]
            X[pr[0] * DS:(pr[0] + 1) * DS, pr[1] * DS:(pr[1] + 1) * DS] += c * (A @ Xf @ C)
    return X
nrng = np.random.default_rng(seed + 14)
X1 = build(nrng.standard_normal(len(basis)) + 1j * nrng.standard_normal(len(basis)))
Rt = [D.R(k, sigJ, DS) for k in tests]
err = max(np.abs(X1 @ R - R @ X1).max() for R in Rt) / np.abs(X1).max()
print(f'[decomp] commutant check on {len(tests)} other letters of B\'_m: relative defect {err:.1e}', flush=True)
H = X1 + X1.conj().T
ev, U = np.linalg.eigh(H)
scale = max(1, np.abs(ev).max())
copies = []; i = 0; gaps = []
while i < n:
    j = i + 1
    while j < n and abs(ev[j] - ev[i]) < 1e-7 * scale: j += 1
    if j < n: gaps.append((ev[j] - ev[j - 1]) / scale)
    copies.append(U[:, i:j]); i = j
spread = max(abs(ev[i + c.shape[1] - 1] - ev[i]) for i, c in zip(np.cumsum([0] + [c.shape[1] for c in copies[:-1]]), copies)) / scale
print(f'[decomp] {len(copies)} copies; min gap between clusters {min(gaps):.1e}, max spread in a cluster {spread:.1e}', flush=True)
X2 = build(nrng.standard_normal(len(basis)) + 1j * nrng.standard_normal(len(basis)))
types = []
for ci, Uc in enumerate(copies):
    for t in types:
        U0 = copies[t[0]]
        if U0.shape[1] == Uc.shape[1] and np.abs(Uc.conj().T @ X2 @ U0).max() > 1e-8 * np.abs(X2).max():
            t.append(ci); break
    else:
        types.append([ci])
Rp = [D.R(p, sigJ, DS) for p in L.p]; Rp2 = [D.R(p, sigJ, DS) for p in p2]
out = {'names': np.array(names), 'offs': offs, 'DS': DS, 'd0': D.d0, 'S': np.array(D.S), 'Sx': np.array(D.Sx),
       'free_d4': np.array([d for d in range(30) if d not in D.S])}
Fr_fixed = [d for d in range(30) if all(act[d][0] == d for act in gact)]
out['Fr_all30'] = np.array(sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(Fr_fixed)))
summ = []; chk = 0.0
for ti, t in enumerate(types):
    U0 = copies[t[0]]; Us = [U0]
    for ck in t[1:]:
        Uk = copies[ck]; Mk = Uk @ (Uk.conj().T @ X2 @ U0)
        Q_, s_, Vh_ = np.linalg.svd(Mk, full_matrices=False); Us.append(Q_ @ Vh_)
    # consistency check: rho_k(g) = Uk^* R(g) Uk must equal rho_0(g) for a test letter
    for R in Rt[:3]:
        r0 = U0.conj().T @ R @ U0
        for Uk in Us[1:]: chk = max(chk, np.abs(Uk.conj().T @ R @ Uk - r0).max())
    na = []; na2 = []
    for a in range(7):
        Ea = sum(L.chars[a][j] * Rp[j] for j in range(7)) / 7
        Ea2 = sum(L.chars[a][j] * Rp2[j] for j in range(7)) / 7
        na.append(int(np.linalg.matrix_rank(U0.conj().T @ Ea @ U0, tol=1e-8)))
        na2.append(int(np.linalg.matrix_rank(U0.conj().T @ Ea2 @ U0, tol=1e-8)))
    # which sigma blocks does the type live in (weight of the copies on each block)
    wts = []
    for i_ in range(len(dss)):
        idx = np.concatenate([np.arange(d * DS + offs[i_], d * DS + offs[i_ + 1]) for d in range(30)])
        wts.append(int(round(sum(np.linalg.norm(Uk[idx]) ** 2 for Uk in Us) / U0.shape[1])))
    out[f'U{ti}'] = np.stack(Us)
    out[f'n{ti}'] = np.array(na); out[f'n2_{ti}'] = np.array(na2)
    summ.append((U0.shape[1], len(t), tuple(na), tuple(na2), tuple(wts)))
out['ntypes'] = len(types)
print(f'[decomp] {len(types)} types; sum d m = {sum(d * m for d, m, *_ in summ)} (n = {n}); sum m^2 = {sum(m * m for _, m, *_ in summ)} '
      f'(commutant dim {len(basis)}); intertwined-basis defect {chk:.1e}', flush=True)
for ti, (d, m, na, na2, wts) in enumerate(summ):
    print(f'  type {ti}: d {d} m {m} n_a(p) {na} n_a(p2) {na2} copies per sigma block {wts}')
np.savez_compressed(outf, **out)
print('done %.0fs' % (time.time() - t0))
