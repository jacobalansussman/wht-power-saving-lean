#!/usr/bin/env python3
"""certify.py (w8-rmexact): exact, independent check of a cell-partition multi-target LP bound (mtlp.py --dump).
Rebuilds every row from its label with its own form/target code (no import of forms.py); takes only the float duals
from the dump.  Rows: ('cover',F): sum_{x in F} c_x - sum_{x in F} p_xx - 2 sum_{x<y in F} p_xy <= a_F;
('del',G): -sum_{x in G} c_x <= -a_G; ('piece<=slack',i,j,G): p_ij - sum_{x in G} c_x <= -a_G (needs i,j in G);
('piece<=c',i,j): p_ij - c_i <= 0.  Variables: c_x (x = 1..2^m-1, index x-1), then p_ij (i<=j) in lexicographic order.
lam = -dual rationalised; bound = -lam.b + U * (sum of negative reduced costs), U = m (valid for every x <= U).
H: 'trivial', 'singer:k' (multiplication by alpha^k in F_2[a]/(poly)), 'orth', 'GL'."""
import sys, json, itertools
from fractions import Fraction as Fr
dump, Hspec, m = sys.argv[1], sys.argv[2], int(sys.argv[3])
poly = int(sys.argv[4], 0) if len(sys.argv) > 4 else (0b10011 if m == 4 else 0b1011)
d = json.load(open(dump)); y = d['dual']; labs = [r[2] for r in d['rows']]
U = Fr(m)
idx = [(i, j) for i in range(m) for j in range(i, m)]
def Qv(Q, x):
    v = 0
    for b, (i, j) in enumerate(idx):
        if (Q >> b) & 1 and (x >> i) & 1 and (x >> j) & 1: v ^= 1
    return v
def gfmul(a, b):
    r = 0
    while b:
        if b & 1: r ^= a
        b >>= 1; a <<= 1
        if a >> m: a ^= poly
    return r
def cols_to_gram(cols):     # h given by images of e_j (columns); return h h^T as form-pairing bits
    S = 0
    for b, (i, j) in enumerate(idx):
        v = 0
        for k in range(m): v ^= ((cols[k] >> i) & 1) & ((cols[k] >> j) & 1)
        if v: S |= 1 << b
    return S
if Hspec == 'trivial': Hcols = [[1 << j for j in range(m)]]
elif Hspec.startswith('singer'):
    k = int(Hspec.split(':')[1]); g = 1
    for _ in range(k): g = gfmul(g, 2)
    Hcols = []; p = 1
    while True:
        Hcols.append([gfmul(p, 1 << j) for j in range(m)]); p = gfmul(p, g)
        if p == 1: break
elif Hspec in ('GL', 'orth'):
    Hcols = []
    for cols in itertools.product(range(1, 1 << m), repeat=m):
        span = {0}
        for c in cols: span |= {s ^ c for s in span}
        if len(span) == 1 << m:
            if Hspec == 'orth' and cols_to_gram(list(cols)) != cols_to_gram([1 << j for j in range(m)]): continue
            Hcols.append(list(cols))
elif Hspec.startswith('rows:'):
    rows = [int(v) for v in Hspec.split(':')[1].split(',')]
    cols0 = [sum(((rows[i] >> j) & 1) << i for i in range(m)) for j in range(m)]   # images of e_j
    def apply_cols(cols, x):
        r = 0
        for j in range(m):
            if (x >> j) & 1: r ^= cols[j]
        return r
    Hcols = []; cur = [1 << j for j in range(m)]
    while True:
        Hcols.append(cur)
        cur = [apply_cols(cols0, c) for c in cur]          # g * cur
        if cur == [1 << j for j in range(m)]: break
n = len(Hcols); T = [cols_to_gram(c) for c in Hcols]
def a(F): return Fr(sum(bin(F & S).count('1') & 1 for S in T), n)
nd = (1 << m) - 1
pairs = [(i, j) for i in range(nd) for j in range(i, nd)]; pk = {p: nd + k for k, p in enumerate(pairs)}
nv = nd + len(pairs)
def supp(F): return [x - 1 for x in range(1, 1 << m) if Qv(F, x)]
lam = [Fr(-v).limit_denominator(10 ** 7) if -v > 1e-11 else Fr(0) for v in y]
rc = [Fr(1)] * nd + [Fr(0)] * len(pairs); bound = Fr(0); used = 0
for lab, l in zip(labs, lam):
    if l == 0: continue
    used += 1; row = {}
    kind = lab[0]
    if kind == 'cover':
        F = int(lab[1]); P = supp(F); rhs = a(F)
        for i in P: row[i] = row.get(i, 0) + 1; row[pk[(i, i)]] = row.get(pk[(i, i)], 0) - 1
        for i, j in itertools.combinations(P, 2): row[pk[(i, j)]] = row.get(pk[(i, j)], 0) - 2
    elif kind == 'del':
        G = int(lab[1]); P = supp(G); rhs = -a(G); row = {i: -1 for i in P}
    elif kind == 'piece<=slack':
        i, j, G = int(lab[1]), int(lab[2]), int(lab[3]); P = supp(G); assert i in P and j in P
        rhs = -a(G); row = {x: -1 for x in P}; row[pk[(min(i, j), max(i, j))]] = 1
    elif kind == 'piece<=c':
        i, j = int(lab[1]), int(lab[2]); rhs = Fr(0); row = {pk[(min(i, j), max(i, j))]: 1, i: -1}
    for k, v in row.items(): rc[k] += l * v
    bound -= l * rhs
neg = sum(min(Fr(0), r) for r in rc)
cert = bound + U * neg
print("%s m=%d |H|=%d targets=%d: dual rows used %d; -lam.b = %s (%.9f); negative reduced costs %s; CERTIFIED M/W >= %s = %.9f" % (
    Hspec, m, n, len(set(T)), used, bound, float(bound), neg, cert, float(cert)))
