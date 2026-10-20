# b56-ceilings (wave 16).  Ceilings of TWO-FACTOR Venn certificates in W = Lambda^a V (x) Lambda^b V, V = C^n, n = 2^m.
#   N_J = span{ e_S (x) e_T : |S cap T| in J }  (kept by every monomial map).
# dir(e1) is monomial in the pair basis D+-, O+- (b4-cert/venn2cut.py), so cap(J) = sum_{c in J, c' not in J} M[c][c'],
# M[c][c'] = #{D/O basis vectors b : c(b) = c, c(Rb) = c'}.  Generating function (this file): per pair
#   empty (S has 0):  (1+t)^2;   S-full (S has 2): (1+tw)^2;   S-half (S has 1):
#   phi = 2 + 2 t^2 w + t w (1 + z^-1) + t (1 + z)     [D-: w, D+: w z^-1, O-: 1, O+: z ;  t = T-count, w = c, z = c'-c]
#   M[c][c+d] = [s^a t^b w^c z^d] ( (1+t)^2 + s phi + s^2 (1+tw)^2 )^(n/2).
# Upper bounds on dem used here (all [proved]):
#   (screen)  dem <= min(dim N, codim N)
#   (refined) Lambda^a (x) Lambda^b = (+)_j S_j, S_j = S_(2^j 1^(a+b-2j)), multiplicity free; N_J is a sum of whole weight
#             spaces, so N_J = (+)_j Y_j with Y_j = S_j cap N_J, and K keeps every S_j:
#             dem(N_J) = sum_j dem(Y_j) <= sum_j min(dim Y_j, dim S_j - dim Y_j),
#             dim Y_j = sum_{c in J} C(n,c) C(n-c,a+b-2c) [C(a+b-2c, j-c) - C(a+b-2c, j-c-1)]   (Kostka numbers)
#   (trivial) dem <= m cap  (K = product of m directions, each conjugate to dir(e1) by a free map)
# Usage: python3 -I venn2ceil.py M [exact]     (exact = Python integers; default float64)
import sys, itertools, time
# (needs numpy; the original loaded a private copy of numpy and scipy here)
import numpy as np, warnings; warnings.filterwarnings("ignore")
from math import comb
from fractions import Fraction

m = int(sys.argv[1]); exact = len(sys.argv) > 2 and sys.argv[2] == 'exact'  # argv[3] = largest a (default n/2)
n = 2 ** m; P = n // 2
dt = object if exact else np.float64

def C(x, y):
    return comb(x, y) if 0 <= y <= x else 0

def chain_for_a(a):
    """A[s=a][t][w][z] for all t: returns array shape (n+1, a+1, 2a+1) (t, c, d+a)"""
    Z = a  # z offset
    A = np.zeros((a + 1, n + 1, a + 1, 2 * a + 1), dtype=dt)
    A[0, 0, 0, Z] = 1
    terms = [  # (coef, ds, dt, dw, dz)
        (1, 0, 0, 0, 0), (2, 0, 1, 0, 0), (1, 0, 2, 0, 0),
        (2, 1, 0, 0, 0), (2, 1, 2, 1, 0), (1, 1, 1, 1, 0), (1, 1, 1, 1, -1), (1, 1, 1, 0, 0), (1, 1, 1, 0, 1),
        (1, 2, 0, 0, 0), (2, 2, 1, 1, 0), (1, 2, 2, 2, 0)]
    for _ in range(P):
        B = np.zeros_like(A)
        for (cf, ds, dtt, dw, dz) in terms:
            src = A[:a + 1 - ds, :n + 1 - dtt, :a + 1 - dw, max(0, -dz):2 * a + 1 - max(0, dz)]
            B[ds:, dtt:, dw:, max(0, dz):2 * a + 1 - max(0, -dz)] += cf * src
        A = B
    return A[a]

def module_dims(a, b):
    N = a + b
    c0 = max(0, N - n); jmax = min(a, b)
    js = list(range(c0, jmax + 1)); cs = list(range(c0, min(a, b) + 1))
    dimS = {j: C(n, N - j) * C(n, j) - C(n, N - j + 1) * C(n, j - 1) for j in js}
    D = {}  # D[j][c] = dim of weight-c part of S_j
    for j in js:
        D[j] = {c: (C(n, c) * C(n - c, N - 2 * c) * (C(N - 2 * c, j - c) - C(N - 2 * c, j - c - 1)) if c <= j else 0) for c in cs}
        assert sum(D[j].values()) == dimS[j], (a, b, j)
    return js, cs, dimS, D

def analyse(a, b, Mab, maxvals_exhaustive=18):
    js, cs, dimS, D = module_dims(a, b)
    cidx = {c: i for i, c in enumerate(cs)}
    k = len(cs)
    M = np.zeros((k, k), dtype=object)
    for c in cs:
        for d in range(-a, a + 1):
            c2 = c + d
            if c2 in cidx:
                M[cidx[c], cidx[c2]] = Mab[c, d + a] if not exact else int(Mab[c, d + a])
    dimc = [sum(D[j][c] for j in js) for c in cs]
    tot = C(n, a) * C(n, b)
    assert sum(dimc) == tot
    if not exact:
        Mf = M.astype(float)
        rows = Mf.sum(axis=1)
        assert np.allclose(rows, np.array(dimc, float), rtol=1e-9), (a, b)
        # chain eigenvalues
        dv = np.array(dimc, float)
        S = Mf / np.sqrt(np.outer(dv, dv))
        ev = np.sort(np.linalg.eigvalsh((S + S.T) / 2))[::-1]
    else:
        for i, c in enumerate(cs):
            assert sum(M[i]) == dimc[i]
        ev = None
    best = None
    if not exact:
        Mf = M.astype(float); dv = np.array(dimc, float)
        Dm = np.array([[float(D[j][c]) for j in js] for c in cs])      # k x nj
        dS = np.array([float(dimS[j]) for j in js])
        if k <= maxvals_exhaustive:
            bits = np.arange(1, (1 << k) - 1, dtype=np.int64)
            X = ((bits[:, None] >> np.arange(k)[None, :]) & 1).astype(float)
        else:
            X = np.array([[1.0 if lo <= i < hi else 0.0 for i in range(k)] for lo in range(k) for hi in range(lo + 1, k + 1)
                          if not (lo == 0 and hi == k)])
        cap = ((X @ Mf) * (1 - X)).sum(1)
        dimN = X @ dv
        Y = X @ Dm
        refined = np.minimum(Y, (1 - X) @ Dm).sum(1)   # complement summed directly: no cancellation
        scr = np.minimum(dimN, (1 - X) @ dv)
        ok = cap > 0.5
        val = np.where(ok, np.minimum(refined, m * cap) / np.where(ok, cap, 1), -1)
        sv = np.where(ok, np.minimum(scr, m * cap) / np.where(ok, cap, 1), -1)
        i = int(np.argmax(val))
        J = tuple(cs[t] for t in range(k) if X[i, t] > 0.5)
        best = (val[i], sv.max(), J, dimN[i], cap[i], refined[i])
        return best, ev
    # J families: exhaustive if few values, else thresholds and intervals
    return best, ev

t0 = time.time()
print(f"m={m} n={n} {'EXACT' if exact else 'float64'}", flush=True)
overall = None
rows = []
AMAX = int(sys.argv[3]) if len(sys.argv) > 3 else P
for a in range(1, AMAX + 1):
    T = chain_for_a(a)
    for b in range(a, n - a + 1):
        Mab = T[b]
        best, ev = analyse(a, b, Mab)
        if best is None: continue
        lam2 = None if ev is None else ev[1]
        rows.append((a, b, best, lam2))
        if overall is None or best[0] > overall[0][0]:
            overall = (best, a, b)
    print(f"  a={a} done ({time.time()-t0:.1f}s)", flush=True)
rows.sort(key=lambda r: -float(r[2][0]))
mu1 = (n - 2) / (2 * (n - 1))
print(f"mu_1 = (n-2)/(2(n-1)) = {mu1:.6f};  spectral ceiling 2/(1-mu_1) = 4(n-1)/n = {4*(n-1)/n:.6f};  hook ceiling 4(n-1)/(n+2) = {4*(n-1)/(n+2):.6f}")
if rows and rows[0][3] is not None:
    lmax = max(r[3] for r in rows)
    print(f"max over (a,b) of the chain's second eigenvalue: {lmax:.9f}  (<= mu_1: {lmax <= mu1 + 1e-12})")
print("top (a,b) by refined ceiling = min(sum_j min(dimY_j, dimS_j-dimY_j), m cap)/cap   [screen = min(dimN,codimN,m cap)/cap]")
for (a, b, best, lam2) in rows[:12]:
    val, sv, J, dimN, cap, refined = best
    print(f"  L^{a} (x) L^{b}: refined {float(val):.4f}  screen {float(sv):.4f}  J={J} dimN {dimN} cap {cap} refined-dem-bound {refined}"
          + ("" if lam2 is None else f"  lambda2 {lam2:.6f}"))
nonhook = [r for r in rows if r[0] >= 2]
if nonhook:
    a, b, best, lam2 = nonhook[0]
    print(f"best with a >= 2 (not a hook): L^{a} (x) L^{b}: refined {float(best[0]):.4f}  screen {float(best[1]):.4f}  J={best[2]}")
    scr = max(nonhook, key=lambda r: float(r[2][1]))
    print(f"largest screen with a >= 2: L^{scr[0]} (x) L^{scr[1]}: screen {float(scr[2][1]):.4f}")
print(f"done {time.time()-t0:.1f}s")
