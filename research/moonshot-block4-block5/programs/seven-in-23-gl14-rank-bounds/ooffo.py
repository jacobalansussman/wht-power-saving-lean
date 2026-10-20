"""ooffo.py -- e3-bpexact: the group-algebra construction in the order O O F F O (7 in 23: O1 O2 F1 F2 O3).
F1 and F2 are both group-algebra shared moves on the SAME O2 term (same move m, cell F), with different order-7 twists
q1, q2 (letters of B_m) and different characters psi1 != psi2; the top O3 reads both.  In array character chi the top's
F_i-part has effective character a_i = chi - psi_i, so F_i's own-character failure (a_i = 0) can be covered by F_j.
(1) Frame-sum relaxation of the actual lifted skeleton with relax.Eng.decide, ALL 7 characters required (E = []):
    no terminal shared move is needed if it passes.
(2) Exact distance-4 test in sigma = V (natural, 16-dim): per chi, the top-alone kernel-block rank with BOTH moves,
    for generic lifts and for flat (Pauli-corrected) lifts.  Need 16 in every chi (no F2 here).
argv: seed psi1 psi2 nBm"""
import sys, os, time, random
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from d4test import *
from relax import Eng
t0 = time.time()
seed, psi1, psi2, nBm = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
L = Lift(seed=seed, psi=psi1)
rng = random.Random(seed + 1); nrng = np.random.default_rng(seed + 3)
base = [b for (b, s) in L.Bm]; Bm = list(base); seen = set(key(b) for b in Bm)
while len(Bm) < nBm:
    w = mm(rng.choice(base), rng.choice(Bm))
    if key(w) not in seen: seen.add(key(w)); Bm.append(w)
mi = minv(L.m); ks = [mm(L.m, b, mi) for b in Bm]
D = D4(L); P = L.P
Pi = lambda k: D.R(k, sig_triv, 1).real
Rt = [Pi(k) for k in ks]
# a second order-7 twist q2 in B_m, not in <q1>
q1 = L.q; pow1 = set(key(x) for x in L.qp)
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
q2 = None
for b in Bm:
    qf = exact7(b)
    if qf is not None and key(qf) not in pow1 and not any(np.allclose(sigma(qf), sigma(x)) for x in L.qp):
        # different cyclic group at the slot level
        P1 = Pi(mm(L.m, q1, mi)); P2 = Pi(mm(L.m, qf, mi))
        if not any(np.allclose(np.linalg.matrix_power(P1, j), P2) for j in range(7)):
            q2 = qf; break
assert q2 is not None
q2p = [I16]
for j in range(1, 7): q2p.append(mm(q2p[-1], q2))
p2 = [mm(L.m, x, mi) for x in q2p]
assert all(is_letter(x) for x in p2)
print(f'[ooffo] second twist found; {time.time() - t0:.0f}s', flush=True)
fixed = [d for d in range(30) if all(R[d, d] == 1 for R in Rt)]
Fr = sorted(set([d for d in range(30) if P.D[D.slot[d], L.lam1] == 1]) | set(fixed))
def design(pp, own):
    """relaxation design of y for twist pp (list of p_j), all characters except the move's own effective one"""
    Rp = [Pi(x) for x in pp]
    gamma = nrng.standard_normal(30) + 1j * nrng.standard_normal(30)
    req = [a for a in range(7) if a != own]
    nb = len(ks); rowsM = []; rhs = []
    for ia, a in enumerate(req):
        Wa = gamma @ sum(L.chars[a][j] * Rp[j] for j in range(7))
        blk = np.zeros((30, nb + len(Fr) * len(req)), dtype=complex)
        for i, R in enumerate(Rt): blk[:, i] = Wa @ R
        for j, d in enumerate(Fr): blk[d, nb + len(Fr) * ia + j] = -1
        e = np.zeros(30, dtype=complex); e[D.d0] = 1
        rowsM.append(blk); rhs.append(e)
    M = np.vstack(rowsM); bb = np.concatenate(rhs)
    sol, *_ = np.linalg.lstsq(M, bb, rcond=None)
    _, s, Vh = np.linalg.svd(M); r = int((s > 1e-9 * s[0]).sum()); Nn = Vh[r:].conj().T
    z = Nn @ (nrng.standard_normal(Nn.shape[1]) + 1j * nrng.standard_normal(Nn.shape[1]))
    sol = sol + z * (np.linalg.norm(sol) / np.linalg.norm(z))
    return sol[:nb], float(np.linalg.norm(M @ sol - bb))
y1, r1 = design(L.p, 0); y2, r2 = design(p2, 0)
print(f'[ooffo] designs: residuals {r1:.1e}, {r2:.1e}', flush=True)
# actual shared moves (formal sums)
def F_of(ys, qp, psi):
    u = {}; ps = np.conj(L.chars[psi])
    for h in range(7):
        for y, b in zip(ys, Bm):
            if abs(y) < 1e-14: continue
            L.fs_add(u, h, mm(L.m, qp[h], b, L.g2, L.cp[h]), ps[h] * y)
    return L.clean(u)
F1 = F_of(y1, L.qp, psi1); F2 = F_of(y2, q2p, psi2)
# absorption check for F2: tau_g(F2) = psi2(g) c_g p2_g^-1 F2
worst = 0
for g in range(7):
    lhs = L.tau(g, F2); rhs = L.lmul(mm(L.c[g], minv(p2[g])), F2)
    rhs = {k: [L.chars[psi2][g] * v[0], v[1]] for k, v in rhs.items()}
    worst = max(worst, L.fs_diff(lhs, rhs))
print(f'[ooffo] F1 {len(F1)} terms (absorption {L.check_absorbed(F1):.0e}), F2 {len(F2)} terms (absorption {worst:.0e})  {time.time() - t0:.0f}s', flush=True)
E = Eng('Z7', P=L.P, seed=seed)
def conv(u):
    out = {}
    for (l, kk), (c, w) in u.items():
        k2 = (l, inv(sigma(w))); out[k2] = out.get(k2, 0) + c
    return out
O1 = {(0, inv(sigma(L.m1))): 1.0 + 0j}; O2 = {(0, inv(sigma(L.g2))): 1.0 + 0j}
F1c, F2c = conv(F1), conv(F2)
for label, shared in (('O O F F O (both shared moves read by the top)', [F1c, F2c]), ('F1 only (= O O F O)', [F1c]), ('F2 only', [F2c])):
    ok, info = E.decide([E.inp(0), O1, O2], shared, E=[])
    print(f'[ooffo] FRAME-SUM, ALL 7 characters, {label}: passes = {ok}; mu_ok {info.get("mu_ok")}', flush=True)
# exact sigma = V, top alone, per array character chi: rows w1 M1 + w2 M2
Delta = np.zeros((480, 480), dtype=complex)
for d in range(30): Delta[16 * d:16 * d + 16, 16 * d:16 * d + 16] = D.g[d]
Di = np.linalg.inv(Delta)
def Kgen(ys): return sum(y * D.R(k, sig_nat, 16) for y, k in zip(ys, ks))
def Kflat(ys): return Delta @ np.kron(sum(y * Pi(k) for y, k in zip(ys, ks)), np.eye(16)) @ Di
RV1 = [D.R(x, sig_nat, 16) for x in L.p]; RV2 = [D.R(x, sig_nat, 16) for x in p2]
cons = [d for d in range(30) if d != D.d0 and d not in Fr]
colx = np.concatenate([np.arange(16 * d, 16 * d + 16) for d in cons]); col0 = np.arange(16 * D.d0, 16 * D.d0 + 16)
def krank(M):
    C = M[:, colx]; U, s, Vh = np.linalg.svd(C, full_matrices=True); r = int((s > 1e-9 * max(1, s[0])).sum())
    Nul = U[:, r:].conj().T; B0 = Nul @ M[:, col0]
    sv = np.linalg.svd(B0, compute_uv=False) if B0.size else np.zeros(1)
    return int((sv > 1e-8 * np.abs(M).max()).sum())
for lift, K1, K2 in (('generic', Kgen(y1), Kgen(y2)), ('flat', Kflat(y1), Kflat(y2))):
    out = []
    for chi in range(7):
        a1, a2 = (chi - psi1) % 7, (chi - psi2) % 7
        M1 = (sum(L.chars[a1][j] * RV1[j] for j in range(7)) / 7) @ K1
        M2 = (sum(L.chars[a2][j] * RV2[j] for j in range(7)) / 7) @ K2
        out.append((chi, krank(M1), krank(M2), krank(np.vstack([M1, M2]))))
    print(f'[ooffo] sigma = V, {lift} lifts, top alone: per chi (rank via F1, via F2, via both; need 16): ' +
          ', '.join(f'chi={c}: {a}/{b}/{ab}' for c, a, b, ab in out), flush=True)
print('done %.0fs' % (time.time() - t0))
