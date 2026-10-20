"""v-rm-theory: own recomputation of rm-theory 5c (block-4 transfer, row 3 of wedge4's cert_4col.txt).
N = span of b(P) = sum over the 4 cosets t+P of e_{t+p1} ^ e_{t+p2} ^ e_{t+p3} ^ e_{t+p4} (P's listed order;
well defined because translating by p in P is an even permutation of P). Planes are read from wedge4's file
directly (own parser). For every unordered basis b of F_2^4 (840): dem(b) = rank_p[N | Lambda^4(K_b) N] - dim N,
K_b = prod 2dir(b_j), entries mod p with i -> sqrt(-1). Rank mod p <= true rank (lower bound on demand: safe).
Also caps per direction mod p (calibration against the file's exact caps) and that shifts keep N.
usage: vk_cert4.py ROW(1-4) p
"""
import itertools, sys, re, time
import numpy as np

row = int(sys.argv[1]); p = int(sys.argv[2])
txt = open(sys.argv[3] if len(sys.argv) > 3 else 'cert_4col.txt').read()
blocks = re.findall(r"planes: (\[\[.*?\]\])\s*\n\s*cap for directions 1..15: (\[[^\]]*\])", txt)
assert len(blocks) == 4, len(blocks)
planes = eval(blocks[row - 1][0]); filecaps = eval(blocks[row - 1][1])
n = 16
quads = list(itertools.combinations(range(n), 4)); qidx = {s: k for k, s in enumerate(quads)}
def perm_sign(seq):
    s = 1; a = list(seq)
    for i in range(len(a)):
        for j in range(i + 1, len(a)):
            if a[i] > a[j]: s = -s
    return s
def bvec(P):
    v = {}
    seen = set()
    for t in range(n):
        coset = frozenset(t ^ x for x in P)
        if coset in seen: continue
        seen.add(coset)
        seq = [t ^ x for x in P]
        k = qidx[tuple(sorted(seq))]
        v[k] = v.get(k, 0) + perm_sign(seq)
    return v
Nv = [bvec(P) for P in planes]
supp = sorted({k for v in Nv for k in v})
r = next(a for a in range(2, p) if pow(a, (p - 1) // 4, p) ** 2 % p == p - 1)
r = pow(r, (p - 1) // 4, p); assert r * r % p == p - 1
def two_dir(w):
    M = np.zeros((n, n), dtype=np.int64)
    for x in range(n):
        M[x, x] = (M[x, x] + 1 + r) % p
        M[x ^ w, x] = (M[x ^ w, x] + 1 - r) % p
    return M % p
def shiftm(a):
    M = np.zeros((n, n), dtype=np.int64)
    for x in range(n): M[x ^ a, x] = 1
    return M
Q = np.array(quads)  # 1820 x 4
PERMS = list(itertools.permutations(range(4))); SG = [perm_sign(s) for s in PERMS]
def wedge_cols(g):
    """Lambda^4(g) columns for e_S, S in supp: matrix 1820 x len(supp), entries det g[T,S] mod p."""
    S = np.array([quads[k] for k in supp])            # s x 4
    out = np.zeros((len(quads), len(supp)), dtype=np.int64)
    for sg, pm in zip(SG, PERMS):
        term = np.ones((len(quads), len(supp)), dtype=np.int64)
        for i in range(4):
            term = (term * g[Q[:, i][:, None], S[:, pm[i]][None, :]]) % p
        out = (out + sg * term) % p
    return out
def Nmat():
    A = np.zeros((len(quads), len(Nv)), dtype=np.int64)
    for j, v in enumerate(Nv):
        for k, c in v.items(): A[k, j] = c % p
    return A
N0 = Nmat()
sidx = {k: i for i, k in enumerate(supp)}
def apply(g):
    WC = wedge_cols(g)
    B = np.zeros((len(quads), len(Nv)), dtype=np.int64)
    for j, v in enumerate(Nv):
        for k, c in v.items(): B[:, j] = (B[:, j] + c * WC[:, sidx[k]]) % p
    return B
def rank_mod(A):
    A = (A.T.copy()) % p; rk = 0
    rows, cols = A.shape
    for c in range(cols):
        nz = np.nonzero(A[rk:, c])[0]
        if len(nz) == 0: continue
        pv = rk + nz[0]
        if pv != rk: A[[rk, pv]] = A[[pv, rk]]
        A[rk] = (A[rk] * pow(int(A[rk, c]), p - 2, p)) % p
        col = A[:, c].copy(); col[rk] = 0
        nzr = np.nonzero(col)[0]
        if len(nzr): A[nzr] = (A[nzr] - np.outer(col[nzr], A[rk])) % p
        rk += 1
        if rk == rows: break
    return rk
t0 = time.time()
dimN = rank_mod(N0)
leak = lambda g: rank_mod(np.hstack([N0, apply(g)])) - dimN
print(f"row {row} p={p}: {len(planes)} planes, dim N = {dimN}; shifts leak: {[leak(shiftm(a)) for a in range(1, n)]}", flush=True)
caps = [leak(two_dir(w)) for w in range(1, n)]
print(f"caps mod p: {caps}; file caps {filecaps}; equal {caps == filecaps}; sum {sum(caps)}", flush=True)
def indep(b):
    span = {0}
    for v in b:
        if v in span: return False
        span = span | {s ^ v for s in span}
    return True
bases = [b for b in itertools.combinations(range(1, n), 4) if indep(b)]
dems = []
for i, b in enumerate(bases):
    g = np.eye(n, dtype=np.int64)
    for w in b: g = (g @ two_dir(w)) % p
    dems.append(leak(g))
    if i % 120 == 0: print(f"  basis {i} {b}: dem {dems[-1]} ({time.time()-t0:.0f}s)", flush=True)
from fractions import Fraction
avg_dem = Fraction(sum(dems), len(dems)); avg_cap = Fraction(sum(filecaps), 15)
hist = {d: dems.count(d) for d in sorted(set(dems))}
print(f"bases {len(bases)}; dem(standard (1,2,4,8)) = {dems[bases.index((1,2,4,8))]}; histogram {hist}")
print(f"avg dem {avg_dem} = {float(avg_dem):.4f}; avg cap {avg_cap}; ratio {avg_dem/avg_cap} = {float(avg_dem/avg_cap):.4f}; {time.time()-t0:.0f}s")
