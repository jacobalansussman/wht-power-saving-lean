"""crit.py -- e3-bpexact: the EXACT kappa-free criterion of NOTE s.4.4 in one representation sigma, one character a:
    T = w' R(kappa) solvable (rows free, kappa in C[B'_m])  iff  rank(T_tau) <= n_(tau,a) for every B'_m-constituent tau.
Steps: (1) commutant C of R_sigma(B'_m) on C^(30 d) from the block relations X_(g.d, g.d') = s(b(d,g))^-1 X_(d,d') s(b(d',g))
(propagation over slot pairs, consistency constraints solved per orbital); (2) a random Hermitian element of C splits
C^(30 d) into irreducible copies; copies are grouped by isomorphism (a random element of C links them); (3) n_(tau,a) =
rank of R(e_a) on one copy; (4) with F = 0 (no help from the free slots) and Z = I: e_tau = rank(T_tau).  e_tau <= n_(tau,a)
for all tau is SUFFICIENT for the character to be served (kappa free); a violation is decisive only if the free slots
cannot lower the rank (reported separately: the rank modulo everything the free slots can add, a lower bound).
argv: seed psi sigma a mode"""
import sys, os, time, random
import numpy as np
import scipy.linalg as _sl
_np_svd = np.linalg.svd
def _svd(a, full_matrices=True, compute_uv=True):      # skeptic patch: LAPACK gesdd sometimes fails; fall back to gesvd
    try: return _np_svd(a, full_matrices=full_matrices, compute_uv=compute_uv)
    except np.linalg.LinAlgError:
        return _sl.svd(a, full_matrices=full_matrices, compute_uv=compute_uv, lapack_driver='gesvd')
np.linalg.svd = _svd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4gen import *
t0 = time.time()
seed, psi, sname, a, mode = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5]
L = Lift(seed=seed, psi=psi)
sig, ds = (sig_nat, 16) if sname == 'V' else make_sig(sname)
rng = random.Random(seed + 13); nrng = np.random.default_rng(seed + 14)
base = [b for (b, s) in L.Bm]
mi = minv(L.m); D = D4(L); P = L.P
gens = [mm(L.m, b, mi) for b in base[:12]] + [L.p[1]]
gact = [D.act(k) for k in gens]
gs = [[sig(beta) for (dd, beta) in act] for act in gact]
gsi = [[np.linalg.inv(x) for x in row] for row in gs]
n = 30 * ds
# (1) commutant by propagation over slot pairs
rep_of = {}; trans = {}          # pair -> (orbital rep, A, C) with X_pair = A X_rep C
cons = {}                        # rep -> list of (A, C) with A X C = X  (consistency)
for d in range(30):
    for e in range(30):
        if (d, e) in rep_of: continue
        rep_of[(d, e)] = (d, e); trans[(d, e)] = (np.eye(ds), np.eye(ds)); cons[(d, e)] = None
        fr = [(d, e)]
        while fr:
            nf = []
            for (x, y) in fr:
                r, A, C = rep_of[(x, y)], *trans[(x, y)]
                for gi, act in enumerate(gact):
                    xx, yy = act[x][0], act[y][0]
                    A2 = gsi[gi][x] @ A; C2 = C @ gs[gi][y]          # X_(g.x, g.y) = s(b(x,g))^-1 X_(x,y) s(b(y,g))
                    if (xx, yy) not in rep_of:
                        rep_of[(xx, yy)] = r; trans[(xx, yy)] = (A2, C2); nf.append((xx, yy))
                    else:
                        A1, C1 = trans[(xx, yy)]
                        rowc = np.kron(A1, C1.T) - np.kron(A2, C2.T)  # A1 X C1 = A2 X C2  (row-major vec)
                        cur = cons[r]
                        cur = rowc if cur is None else np.vstack([cur, rowc])
                        if cur.shape[0] > 2 * ds * ds:
                            _, s_, Vh_ = np.linalg.svd(cur, full_matrices=False)
                            cur = Vh_[:int((s_ > 1e-10 * max(1, s_[0])).sum())] * s_[:int((s_ > 1e-10 * max(1, s_[0])).sum())][:, None]
                        cons[r] = cur
            fr = nf
reps = list(cons)
basis = []                       # list of (rep, X_rep) solutions
for r in reps:
    if cons[r] is not None:
        M = cons[r]; _, s, Vh = np.linalg.svd(M); k = int((s > 1e-9 * max(1, s[0])).sum())
        sols = Vh[k:].conj()
    else:
        sols = np.eye(ds * ds)
    for v in sols: basis.append((r, v.reshape(ds, ds)))
print(f'[crit] {sname} (dim {ds}) a={a}: {len(reps)} orbitals of slot pairs, commutant dim {len(basis)}  {time.time() - t0:.0f}s', flush=True)
def build(coefs):
    X = np.zeros((n, n), dtype=complex)
    for c, (r, Xr) in zip(coefs, basis):
        if abs(c) < 1e-15: continue
        for pr, rr in rep_of.items():
            if rr != r: continue
            A, C = trans[pr]
            X[pr[0] * ds:(pr[0] + 1) * ds, pr[1] * ds:(pr[1] + 1) * ds] += c * (A @ Xr @ C)
    return X
X1 = build(nrng.standard_normal(len(basis)) + 1j * nrng.standard_normal(len(basis)))
# check commutation with a few other letters
Rchk = [D.R(mm(L.m, b, mi), sig, ds) for b in base[20:24]]
err = max(np.abs(X1 @ R - R @ X1).max() for R in Rchk) / np.abs(X1).max()
print(f'[crit] commutant check on 4 other letters: relative defect {err:.1e}', flush=True)
H = X1 + X1.conj().T
ev, U = np.linalg.eigh(H)
# group eigenvalues into copies (clusters)
copies = []; i = 0
while i < n:
    j = i + 1
    while j < n and abs(ev[j] - ev[i]) < 1e-7 * max(1, abs(ev).max()): j += 1
    copies.append(U[:, i:j]); i = j
X2 = build(nrng.standard_normal(len(basis)) + 1j * nrng.standard_normal(len(basis)))
types = []
for ci, Uc in enumerate(copies):
    for t in types:
        U0 = copies[t[0]]
        if U0.shape[1] == Uc.shape[1] and np.abs(Uc.conj().T @ X2 @ U0).max() > 1e-8 * np.abs(X2).max():
            t.append(ci); break
    else:
        types.append([ci])
Rp = [D.R(p, sig, ds) for p in L.p]
Ea = sum(L.chars[a][j] * Rp[j] for j in range(7)) / 7
Fr_fixed = [d for d in range(30) if all(act[d][0] == d for act in gact)]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(Fr_fixed))
free = Fr if mode == 'all30' else [d for d in range(30) if d not in D.S]
ok_all = True; lines = []
for t in types:
    dt = copies[t[0]].shape[1]; m = len(t)
    # consistent bases: U_k = normalise( P_k X2 U_0 )  (an intertwiner from copy 0 to copy k)
    U0 = copies[t[0]]; Us = [U0]
    for ck in t[1:]:
        Uk = copies[ck]; Mk = Uk @ (Uk.conj().T @ X2 @ U0)
        # make it an isometry onto copy k: polar factor
        Q_, s_, Vh_ = np.linalg.svd(Mk, full_matrices=False); Us.append(Q_ @ Vh_)
    nta = int(np.linalg.matrix_rank(U0.conj().T @ Ea @ U0, tol=1e-8))
    # T_tau with F = 0, Z = I: columns (r, k) = (E_d0 row r) U_k  in C^dt
    cols = [U[D.d0 * ds + r, :] for U in Us for r in range(ds)]
    Tt = np.array(cols).T                       # dt x (m ds)
    et = int(np.linalg.matrix_rank(Tt, tol=1e-8))
    # what the free slots can add (union over all free coordinates): lower bound on the rank after any F
    Fsp = np.array([U[s * ds + r, :] for U in Us for s in free for r in range(ds)]).T if free else np.zeros((dt, 0))
    if Fsp.shape[1]:
        Qf, sf, _ = np.linalg.svd(Fsp, full_matrices=False); Qf = Qf[:, :int((sf > 1e-8 * max(1, sf[0])).sum())]
        Tq = Tt - Qf @ (Qf.conj().T @ Tt); lb = int(np.linalg.matrix_rank(Tq, tol=1e-8))
    else:
        lb = et
    lines.append((dt, m, nta, et, lb))
    if et > nta: ok_all = False
viol = [x for x in lines if x[3] > x[2]]
dec = [x for x in lines if x[4] > x[2]]
print(f'[crit] {sname} a={a} {mode}: {len(copies)} irreducible copies in {len(types)} types; (d_tau, m_tau, n_tau_a, e_tau(F=0), lower bound '
      f'with F) per type: ' + ' '.join(f'({dt},{m},{nt},{et},{lb})' for dt, m, nt, et, lb in lines), flush=True)
print(f'[crit] {sname} a={a} {mode}: F = 0 suffices: {ok_all}; types violating with F = 0: {len(viol)}; violating even after any F '
      f'(lower bound > n): {len(dec)}  {time.time() - t0:.0f}s', flush=True)
print('done %.0fs' % (time.time() - t0))
