# Exact rational f(H_0) for Lambda^a (x) Lambda^b, J (skeptic's own code). H_0 = unnormalised Hadamard on address bit 0;
# dir(e1) = I - iX_e1 = D H_0 D with D = diag(1,-i) on bit 0 (a free phase), so f(dir e1) = f(H_0) exactly.
# The off-block (rows outside N, cols in N) is block diagonal after permutation (H_0 keeps the pair pattern), so its exact
# rank = sum of exact ranks (Fractions) over connected components of its nonzero pattern.
import sys, itertools
from fractions import Fraction
a, b = int(sys.argv[1]), int(sys.argv[2]); J = {int(x) for x in sys.argv[3].split(',')}
n = 8
H0 = [[0]*n for _ in range(n)]
for x in range(n):
    y = x ^ 1
    H0[x][x] = -1 if (x & 1) else 1
    H0[y][x] = 1
def sgn(pm):
    s = 1
    for i in range(len(pm)):
        for j in range(i+1, len(pm)):
            if pm[i] > pm[j]: s = -s
    return s
PER = {k: [(pm, sgn(pm)) for pm in itertools.permutations(range(k))] for k in range(1, 8)}
def det(R, C):
    t = 0
    for pm, s in PER[len(R)]:
        pr = s
        for i in range(len(R)):
            pr *= H0[R[i]][C[pm[i]]]
            if pr == 0: break
        t += pr
    return t
def compound(k):
    S = list(itertools.combinations(range(n), k))
    return S, {(i, j): v for i, R in enumerate(S) for j, Cc in enumerate(S) if (v := det(R, Cc)) != 0}
Sa, Ca = compound(a); Sb, Cb = compound(b)
nb = len(Sb)
inN = lambda i, j: len(set(Sa[i]) & set(Sb[j])) in J
# nonzeros of Lambda^a(H0) (x) Lambda^b(H0) in (out-of-N row, in-N col)
rowsA = {}; 
for (i, j), v in Ca.items(): rowsA.setdefault(j, []).append((i, v))
rowsB = {}
for (i, j), v in Cb.items(): rowsB.setdefault(j, []).append((i, v))
entries = {}
for ca in range(len(Sa)):
    for cb in range(nb):
        if not inN(ca, cb): continue
        for (ra, va) in rowsA.get(ca, []):
            for (rb, vb) in rowsB.get(cb, []):
                if not inN(ra, rb):
                    entries[((ra, rb), (ca, cb))] = va * vb
# union-find components
par = {}
def find(x):
    while par.setdefault(x, x) != x:
        par[x] = par[par[x]]; x = par[x]
    return x
for (r, c) in entries: 
    rr, cc = find(('r', r)), find(('c', c))
    if rr != cc: par[rr] = cc
comps = {}
for (r, c), v in entries.items(): comps.setdefault(find(('c', c)), []).append((r, c, v))
def exact_rank(ent):
    rows = sorted({r for r, c, v in ent}); cols = sorted({c for r, c, v in ent})
    ri = {r: i for i, r in enumerate(rows)}; ci = {c: i for i, c in enumerate(cols)}
    M = [[Fraction(0)] * len(cols) for _ in rows]
    for r, c, v in ent: M[ri[r]][ci[c]] = Fraction(v)
    rk = 0
    for c in range(len(cols)):
        piv = next((i for i in range(rk, len(rows)) if M[i][c] != 0), None)
        if piv is None: continue
        M[rk], M[piv] = M[piv], M[rk]
        for i in range(len(rows)):
            if i != rk and M[i][c] != 0:
                f = M[i][c] / M[rk][c]
                M[i] = [x - f * y for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk, len(rows), len(cols)
tot = 0; mx = 0
for k, ent in comps.items():
    rk, nr, nc = exact_rank(ent); tot += rk; mx = max(mx, nr, nc)
dimN = sum(1 for i in range(len(Sa)) for j in range(nb) if inN(i, j))
print(f"Lambda^{a} (x) Lambda^{b} J={sorted(J)}: dim N {dimN}; components {len(comps)} (largest side {mx}); EXACT f(H_0) = cap = {tot}")
