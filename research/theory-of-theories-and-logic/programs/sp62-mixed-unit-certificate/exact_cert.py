# exact_cert.py -- exact (integer) certification of a cut certificate N = X^{O} (the sum of the A-weight spaces of X for the
# characters psi of A = Sym(3,2) in one GL(3,2)-orbit O) in an irreducible X of Sp(6,2) sitting with multiplicity 1 in the
# monomial module M = Ind_B^G(lambda).  Everything is an integer vector of Z^2835 acted on by signed permutation matrices.
#  1. u = sum_g chi(g) rho(g) e_F is read off the float projection (times |G|/chi(1)) and ROUNDED; then u in X is checked
#     EXACTLY: (Z2a - z2a) u = 0 and (Z2b - z2b) u = 0 (class sums of 2a, 2b; integer central characters), and the pair
#     (z2a, z2b) is checked to occur for no other constituent of M.
#  2. N is spanned by 64 Pi_psi rho(g) u (psi in O, g in G): integer vectors in X^{O}; X is spanned by rho(g) u.
#  3. 105 rows R on which a basis of X is invertible (mod p => over Q): restriction to R is injective on X, so
#     rank_Q [N, gN] = rank_Q [N, gN]_R, computed EXACTLY by fraction-free (Bareiss) elimination in Python integers.
# Usage: exact_cert.py lam J orbit_size [J2 orbit2 ...]     e.g. exact_cert.py 3 8 28 11 28
import sys, json, time, random
import numpy as np
args = sys.argv[1:]; sys.argv = sys.argv[:1]
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fgcore import *
lam = int(args[0]); jobs = [(int(args[i]), int(args[i + 1])) for i in range(1, len(args), 2)]
P31 = 2147483647
t0 = time.time()
def apply(Ps, y):
    P, s = Ps; out = np.zeros_like(y); out[P] = s.astype(np.int64)[:, None] * y if y.ndim == 2 else s.astype(np.int64) * y; return out
def rank_mod(M, p=P31, want_pivots=False):
    M = np.array(M, dtype=np.int64) % p; r = 0; piv = []; rows, cols = M.shape
    for c in range(cols):
        nz = np.nonzero(M[r:, c])[0]
        if len(nz) == 0: continue
        pr = r + nz[0]; M[[r, pr]] = M[[pr, r]]
        iv = pow(int(M[r, c]), p - 2, p); M[r] = (M[r] * iv) % p
        col = M[:, c].copy(); col[r] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): M[nzr] = (M[nzr] - (col[nzr, None] * M[r][None, :]) % p) % p
        piv.append(c); r += 1
        if r == rows: break
    return (r, piv) if want_pivots else r
def bareiss_rank(rows):
    A = [list(map(int, r)) for r in rows]; n = len(A); m = len(A[0]); r = 0; prev = 1
    for c in range(m):
        pr = next((i for i in range(r, n) if A[i][c] != 0), None)
        if pr is None: continue
        A[r], A[pr] = A[pr], A[r]
        for i in range(r + 1, n):
            if A[i][c] == 0 and True:
                row = []
                for k in range(m):
                    if k < c: row.append(0); continue
                    num = A[r][c] * A[i][k]
                    assert num % prev == 0, "Bareiss division not exact"
                    row.append(num // prev)
                A[i] = row; continue
            row = []
            for k in range(m):
                if k < c: row.append(0); continue
                num = A[r][c] * A[i][k] - A[i][c] * A[r][k]
                assert num % prev == 0, "Bareiss division not exact"
                row.append(num // prev)
            A[i] = row
        prev = A[r][c]; r += 1
        if r == n: break
    return r
MA = [mono(g, lam) for g in TR]; MB = [mono(g, lam) for g in C2B]
Za = np.zeros((NF, NF), np.int64); Zb = np.zeros((NF, NF), np.int64)
for P, s in MA: Za[P, np.arange(NF)] += s.astype(np.int64)
for P, s in MB: Zb[P, np.arange(NF)] += s.astype(np.int64)
alpha = 0.6180339887
w, V = np.linalg.eigh(Za + alpha * Zb)
cl = []; st = 0
for i in range(1, NF + 1):
    if i == NF or w[i] - w[i - 1] > 1e-6: cl.append((st, i)); st = i
consts = []
for (a, b) in cl:
    Q = V[:, a:b]; ra = np.trace(Q.T @ Za @ Q) / (b - a); rb = (w[a:b].mean() - ra) / alpha
    consts.append((round(ra, 6), round(rb, 6), b - a, a, b))
print("lambda", lam, "constituents (z2a, z2b, dim):", [(c[0], c[1], c[2]) for c in consts], f"({time.time()-t0:.0f}s)", flush=True)
# A-characters and GL(3,2)-orbits on them
Ael = {}
for bits in range(64):
    g = np.arange(NV)
    for i in range(6):
        if bits >> i & 1: g = comp(g, A_GENS[i])
    Ael[key(g)] = bits
def abits(g): return Ael[key(g)]
def act(l, s):   # (l.psi)(a) = psi(l^-1 a l): new bit j = s . bits(l^-1 a_j l)
    li = inv(l); return sum(((bin(s & abits(comp(li, comp(A_GENS[j], l)))).count('1') & 1) << j) for j in range(6))
orbits = []; seen = set()
for s in range(64):
    if s in seen: continue
    orb = {s}; fr = [s]
    while fr:
        nx = []
        for x in fr:
            for l in L_GENS:
                y = act(l, x)
                if y not in orb: orb.add(y); nx.append(y)
        fr = nx
    seen |= orb; orbits.append(sorted(orb))
print("GL(3,2)-orbits on characters of A, sizes:", sorted(len(o) for o in orbits), flush=True)
MAg = [mono(g, lam) for g in A_GENS]
def Pi(s, y):     # 64 * projector onto the psi_s weight space, psi_s(a_i) = (-1)^{s_i}
    for i in range(6):
        y = y + (-1 if s >> i & 1 else 1) * apply(MAg[i], y)
    return y
rng = random.Random(1)
def rand_elt(n=12):
    g = np.arange(NV)
    for _ in range(n): g = comp(GENS[rng.randrange(len(GENS))], g)
    return g
results = []; SAVE = {}
for J, osz in jobs:
    d, va, vb = TABLE[J]; z2a, z2b = 63 * va // d, 315 * vb // d
    assert 63 * va % d == 0 and 315 * vb % d == 0
    hit = [c for c in consts if abs(c[0] - z2a) < 1e-6 and abs(c[1] - z2b) < 1e-6]
    assert len(hit) == 1 and hit[0][2] == d, ("not multiplicity one / not separated", hit)
    a, b = hit[0][3], hit[0][4]; Q = V[:, a:b]
    F = int(np.argmax((Q * Q).sum(axis=1))); v0 = Q @ Q[F]; c = 1451520 // d
    u = np.rint(c * v0).astype(np.int64); err = np.abs(c * v0 - u).max()
    assert err < 1e-6, err
    ok_member = bool((Za @ u == z2a * u).all() and (Zb @ u == z2b * u).all()) and u.any()
    print(f"\nirreducible {d}({va:+d}): u = sum_g chi(g) rho(g) e_F rounded (max rounding error {err:.1e}); exact membership u in X: {ok_member}; "
          f"(z2a, z2b) = ({z2a}, {z2b}) unique among constituents: True", flush=True)
    assert ok_member
    # basis of X
    Bcols = [u]; rB = 1
    while rB < d:
        g = rand_elt(); y = apply(mono(g, lam), u)
        r2 = rank_mod(np.array(Bcols + [y]))
        if r2 > rB: Bcols.append(y); rB = r2
    Bm = np.array(Bcols).T
    rB, piv = rank_mod(Bm.T, want_pivots=True); R = piv; assert rB == d and len(R) == d
    # N = X^O = sum over psi in O of Pi_psi X, spanned by the columns of Pi_psi B (B a basis of X)
    orb = next(o for o in orbits if len(o) == osz)
    wdim = {}
    for s_ in range(64):
        wdim[s_] = rank_mod(Pi(s_, Bm).T)
    assert sum(wdim.values()) == d, ("weight dims do not add up", sum(wdim.values()))
    Ncols = []; rN = 0
    for s_ in orb:
        Y = Pi(s_, Bm)
        for c in range(Y.shape[1]):
            y = Y[:, c]
            if not y.any(): continue
            r2 = rank_mod(np.array(Ncols + [y]))
            if r2 > rN: Ncols.append(y); rN = r2
    target = sum(wdim[s_] for s_ in orb)
    assert rN == target, (rN, target)
    print(f"  A-weight space dims (mod p, summing to {d}): orbit-wise", {len(o): sum(wdim[x] for x in o) for o in orbits}, flush=True)
    Nm = np.array(Ncols).T; dimN = rN
    NR = Nm[R]
    dN_exact = bareiss_rank(NR.T.tolist())
    prof = []
    for gi, g in enumerate(HS):
        gN = apply(mono(g, lam), Nm)
        assert (gN[np.setdiff1d(np.arange(NF), R)] is not None)
        M = np.concatenate([Nm, gN], axis=1)[R]
        SAVE.setdefault(f'{d}({va:+d})', {})[['h', 'h2', 'k'][gi]] = M.tolist(); SAVE[f'{d}({va:+d})']['N'] = NR.tolist(); SAVE[f'{d}({va:+d})']['B_R'] = Bm[R].tolist()
        rex = bareiss_rank(M.T.tolist()); rp = rank_mod(np.concatenate([Nm, gN], axis=1).T)
        prof.append((rex - dN_exact, rp - dimN))
    psi = [p[0] for p in prof]
    print(f"  dim N = {dN_exact} (exact; mod p {dimN}); orbit size {osz}; profile exact (psi1, psi2, psi3) = {tuple(psi)}; "
          f"mod-p profile {tuple(p[1] for p in prof)}; psi3/psi1 = {psi[2]}/{psi[0]}; psi3 + 4 psi2 - 10 psi1 = {psi[2] + 4*psi[1] - 10*psi[0]}; "
          f"7 psi3 + 38 psi2 - 88 psi1 = {7*psi[2] + 38*psi[1] - 88*psi[0]}  ({time.time()-t0:.0f}s)", flush=True)
    results.append(dict(J=J, deg=d, chi2a=va, lam=lam, orbit=osz, dimN=dN_exact, psi=psi, max_abs_entry=int(np.abs(Nm).max())))
json.dump(results, open(f"out_exact_cert_lam{lam}.json", "w"))
json.dump(SAVE, open(f"exact_mats_lam{lam}.json", "w"))
