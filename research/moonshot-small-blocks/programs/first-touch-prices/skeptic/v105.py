# v105.py -- v-w7-mixed3: independent rebuild of w7-mixed3's exhaustive single-term defect table in Theorem M's 105.
# Model: W(E7)^+ = Sp(6,2) from e7model.py (written by v-meta2-finite-group-compute; nothing from c2-finite-group or
# w7-mixed3 is used).  X = the 105a inside W1 = 7 (x) Lambda^2 7 (eigenvalue -21 of the transvection class sum).
# N = sum of A-weight spaces over the 28-orbit (Theorem M's N).  E_T = N cap T^-1 N.  defect_T(g) = rank of
# (X -> X/N) o rho(h g) on E_T, for ALL g in P0 = A . GL(3,2) (10752), EXACT mod two primes; stable ranks in floats in an
# orthonormalised model.  Also: number of distinct defect images among defects <= 4.
import itertools, time, sys, collections, numpy as np
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from e7model import *
t0 = time.time()
A_ws = [[z_(1)], [z_(2)], [z_(3)], [z_(1) | z_(2)], [z_(1) | z_(3)], [z_(2) | z_(3)]]
L_ws = [[1 << i, 1 << (3 + j), (1 << i) | (1 << (3 + j))] for i, j in itertools.permutations(range(3), 2)]
h_w, h2_w, k_w = [ee(1)], [ee(1), ee(2)], [ee(1), ee(2), ee(3)]
L0 = frozenset(v for v in range(64) if v & 7 == 0)
for w in A_ws + L_ws: assert frozenset(word_fn(w)(v) for v in L0) == L0      # all in P0 = Stab(L0)
pairs = list(itertools.combinations(range(7), 2))
def wedge2(R):
    out = np.zeros((21, 21), dtype=R.dtype)
    for a, (i, j) in enumerate(pairs):
        for b, (k, l) in enumerate(pairs):
            out[a, b] = R[i, k] * R[j, l] - R[i, l] * R[j, k]
    return out
def W1(R): return np.kron(R, wedge2(R))
# ---- P0 as 7x7 integer matrices: A (64) and Levi (BFS, expect 168) ----
A7 = [word_R7(w) for w in A_ws]; L7 = [word_R7(w) for w in L_ws]
Aall7 = []
for bits in range(64):
    a = np.eye(7, dtype=np.int64)
    for j in range(6):
        if bits >> j & 1: a = a @ A7[j]
    Aall7.append(a)
assert len({a.tobytes() for a in Aall7}) == 64
Lev7 = [np.eye(7, dtype=np.int64)]; keys = {Lev7[0].tobytes()}; fr = list(Lev7)
while fr:
    nx = []
    for g in fr:
        for l in L7:
            gg = l @ g
            if gg.tobytes() not in keys: keys.add(gg.tobytes()); Lev7.append(gg); nx.append(gg)
    fr = nx
P0keys = {(a @ l).tobytes() for a in Aall7 for l in Lev7}
print('|A| 64, |Levi|', len(Lev7), '|A.Levi| distinct', len(P0keys), flush=True)
# closure check: P0 = A.Levi is a group (products of generators stay inside)
for g in [a @ l for a in Aall7[:8] for l in Lev7[:20]]:
    for s in A7 + L7: assert (s @ g).tobytes() in P0keys and (g @ s).tobytes() in P0keys
h7, h27, k7 = word_R7(h_w), word_R7(h2_w), word_R7(k_w)
D = 147
Z = sum(W1(R7_t(w)) for w in range(1, 64))
# ---------------- exact part, mod p ----------------
def rref(M, p):
    M = np.array(M, dtype=np.int64) % p; rows, cols = M.shape; r = 0; piv = []
    for c in range(cols):
        if r == rows: break
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]
        if pr != r: M[[r, pr]] = M[[pr, r]]
        M[r] = (M[r] * pow(int(M[r, c]), p - 2, p)) % p
        col = M[:, c].copy(); col[r] = 0; nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        piv.append(c); r += 1
    return M[:r], piv
def rank(M, p): return len(rref(M, p)[1])
def nullspace(M, p):        # columns spanning {x : M x = 0}
    R, piv = rref(M, p); n = M.shape[1]; free = [c for c in range(n) if c not in piv]
    B = np.zeros((n, len(free)), dtype=np.int64)
    for j, f in enumerate(free):
        B[f, j] = 1
        for i, c in enumerate(piv): B[c, j] = (-R[i, f]) % p
    return B
def mm(A, B, p): return ((np.asarray(A) % p) @ (np.asarray(B) % p)) % p
def isprime(n): return n > 1 and all(n % d for d in range(2, int(n ** 0.5) + 1))
PR = [q for q in range(2**25 - 1, 2**25 - 400, -1) if isprime(q)][:2]
Arho = [W1(a) for a in A7]
for p in PR:
    I = np.eye(D, dtype=np.int64)
    Pp = I.copy()
    for c in (-45, -9): Pp = mm(Pp, (Z - c * I) % p, p)
    _, pv = rref(Pp.T.copy(), p)       # independent columns of Pp
    rX = rank(Pp, p); assert rX == 105, rX
    Rr, pvc = rref(Pp, p); BX = Pp[:, [c for c in pvc]] if False else None
    # columns of Pp spanning X: pivot columns of Pp
    _, pcol = rref(Pp, p); BX = Pp[:, pcol]; assert BX.shape[1] == 105
    def Pi(s, Y):
        for i in range(6): Y = (Y + (-1 if s >> i & 1 else 1) * mm(Arho[i], Y, p)) % p
        return Y
    # characters of A and the Levi action: psi_s(a_j) = (-1)^{s_j}; l.psi(a) = psi(l^-1 a l); find via weight spaces
    wsp = {s: Pi(s, BX) for s in range(64)}
    wd = {s: rank(wsp[s], p) for s in range(64)}; assert sum(wd.values()) == 105
    # orbits of the Levi on characters, computed from the 7-dim matrices: l^-1 a_j l is an element of A; read its bits
    Aidx = {a.tobytes(): b for b, a in enumerate(Aall7)}
    def act(l, s):
        li = np.round(np.linalg.inv(l)).astype(np.int64); out = 0
        for j in range(6):
            c = Aidx[(li @ A7[j] @ l).tobytes()]
            out |= (bin(s & c).count('1') & 1) << j
        return out
    orbs, seen = [], set()
    for s in range(64):
        if s in seen: continue
        orb = {s}; frn = [s]
        while frn:
            nx = []
            for x in frn:
                for l in L7:
                    y = act(l, x)
                    if y not in orb: orb.add(y); nx.append(y)
            frn = nx
        seen |= orb; orbs.append(sorted(orb))
    O28 = [o for o in orbs if len(o) == 28]; assert len(O28) == 1; O28 = O28[0]
    Ncols = np.concatenate([wsp[s] for s in O28], axis=1); _, pn = rref(Ncols, p); NB = Ncols[:, pn]
    dN = NB.shape[1]
    # annihilator of N in the 147-dim ambient: rows y with y^T NB = 0
    Q = nullspace(NB.T.copy(), p).T          # (147 - dN) x 147
    assert rank(Q, p) == D - dN
    Lrho = [W1(l) for l in Lev7]
    for M in Arho + [W1(l) for l in L7]: assert rank(mm(Q, mm(M, NB, p), p), p) == 0     # N is P0-stable
    hW, h2W, kW = W1(h7), W1(h27), W1(k7)
    prof = [rank(mm(Q, mm(M, NB, p), p), p) for M in (hW, h2W, kW)]
    print(f'p={p}: dim X 105, dim N {dN}, profile {tuple(prof)}', flush=True)
    E = {}
    for nm, T in (('kernel', kW), ('two-row', h2W)):
        cf = nullspace(mm(Q, mm(T, NB, p), p), p); E[nm] = mm(NB, cf, p)
        print(f'  dim E_{nm} = {E[nm].shape[1]}', flush=True)
    QhA = [mm(Q, mm(hW, W1(a), p), p) for a in Aall7]
    for nm in ('kernel', 'two-row'):
        LE = [mm(W1(l), E[nm], p) for l in Lev7]
        cnt = collections.Counter(); imgs = set()
        for ai in range(64):
            for li in range(len(Lev7)):
                M = mm(QhA[ai], LE[li], p); R, piv = rref(M.T.copy(), p)   # row space of M^T = image of M (in Q coords)
                r = len(piv); cnt[r] += 1
                if r <= 4: imgs.add(R.tobytes())
        print(f'  {nm}: defects over all {64*len(Lev7)} elements of P0: {dict(sorted(cnt.items()))}; '
              f'distinct images (defect <= 4): {len(imgs)}  ({time.time()-t0:.0f}s)', flush=True)
# ---------------- float part: stable ranks, orthonormalised model ----------------
Lc = np.linalg.cholesky(C.astype(float))      # C = Lc Lc^T
def orth(R): return Lc.T @ R @ np.linalg.inv(Lc.T)
def W1f(R): Ro = orth(R.astype(float)); return np.kron(Ro, wedge2(Ro))
Zf = sum(W1f(R7_t(w)) for w in range(1, 64)); assert np.abs(Zf - Zf.T).max() < 1e-8
ev, U = np.linalg.eigh(Zf); BXf = U[:, np.abs(ev + 21) < 1e-6]; assert BXf.shape[1] == 105
Af = [W1f(a) for a in A7]
for M in Af: assert np.abs(M @ M.T - np.eye(D)).max() < 1e-9
def Pif(s):
    P = np.eye(D)
    for i in range(6): P = P @ (np.eye(D) + (-1 if s >> i & 1 else 1) * Af[i]) / 2
    return P
PNf = sum(Pif(s) for s in O28) @ (BXf @ BXf.T)
w_, V_ = np.linalg.eigh((PNf + PNf.T) / 2); Nf = V_[:, w_ > 0.5]; assert Nf.shape[1] == 56
PN = Nf @ Nf.T
def Ef(T):
    G = T @ Nf; R = G - PN @ G; u, sv, vh = np.linalg.svd(R); r = int((sv > 1e-8).sum()); return Nf @ vh[r:].T
hf = W1f(h7); Ef_ = {'kernel': Ef(W1f(k7)), 'two-row': Ef(W1f(h27))}
Aall_f = [W1f(a) for a in Aall7]; Lev_f = [W1f(l) for l in Lev7]
for nm, E_ in Ef_.items():
    by = collections.defaultdict(lambda: [99.0, 0.0, 0]); gap = [1e9, 0.0]
    for a in Aall_f:
        ha = hf @ a
        for l in Lev_f:
            W = ha @ (l @ E_); W = W - PN @ W
            s = np.linalg.svd(W, compute_uv=False); r = int((s > 1e-7).sum())
            if r: gap[0] = min(gap[0], s[r - 1])
            if r < len(s): gap[1] = max(gap[1], s[r])
            if r == 0: continue
            sr = float((s ** 2).sum() / s[0] ** 2); st = by[r]; st[0] = min(st[0], sr); st[1] = max(st[1], sr); st[2] += 1
    print(f'float {nm}: dim E {E_.shape[1]}; (rank: min srank, max srank, count)',
          {r: (round(v[0], 4), round(v[1], 4), v[2]) for r, v in sorted(by.items())},
          f'; smallest kept sv {gap[0]:.3e}, largest dropped sv {gap[1]:.3e}  ({time.time()-t0:.0f}s)', flush=True)
