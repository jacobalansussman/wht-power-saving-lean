# vshannon2.py -- v-nc-gain: my own Shannon LP on nc-gain's 3-cell contraction of Gamma_2 (a line through frame 0 /
# frames opposite frame 0 / the rest), built from MY Gamma_2. Grouped elements (valid: a group's h is at most the sum
# of its members'): one source element per ordered pair of cells (rate = number of opposite pairs from one to the
# other), one arc element per ordered pair of adjacent cells. Rows: elemental Shannon, independent sources with their
# rates, closure for EVERY union of cells. Objective min sum of arc-group entropies; claim (nc-gain job_cells): 80.
import sys, itertools, collections
import numpy as np, scipy.sparse as sp
from scipy.optimize import linprog
class G: pass
g = G(); g.D = np.load('gamma2_D_dpgidx.npy').astype(np.int16); n = g.n = g.D.shape[0]; m = 2
v = 0; nb = [int(x) for x in np.nonzero(g.D[v] == 1)[0]]
x = next(a for a in nb if any(g.D[a, b] == 1 for b in nb if b != a)); y = next(b for b in nb if b != x and g.D[x, b] == 1)
line = {v, x, y}
part = [0 if u in line else (1 if g.D[v, u] == 2 else 2) for u in range(n)]
k = 3
arcsG = collections.Counter((part[a], part[b]) for a in range(n) for b in range(n) if g.D[a, b] == 1 and part[a] != part[b])
dem = collections.Counter((part[s], part[t]) for s in range(n) for t in range(n) if g.D[s, t] == m and part[s] != part[t])
sess = sorted(dem); arcs = sorted(arcsG); ns = len(sess); na = len(arcs); N = 1 << (ns + na)
print(flush=True); print("cells sizes", collections.Counter(part), "arc groups", dict(arcsG), "session groups", dict(dem))
Y = lambda i: 1 << i; X = lambda i: 1 << (ns + i); ne = ns + na; full = N - 1
I, J, V = [], [], []; r = 0
def row(d):
    global r
    for msk, cf in d.items():
        if msk and cf: I.append(r); J.append(msk - 1); V.append(cf)
    r += 1
for i in range(ne): row({full: 1, full ^ (1 << i): -1})
for i, j in itertools.combinations(range(ne), 2):
    rest = [t for t in range(ne) if t not in (i, j)]
    for K in range(1 << len(rest)):
        kk = 0
        for b in range(len(rest)):
            if K >> b & 1: kk |= 1 << rest[b]
        d = collections.Counter({kk | 1 << i: 1, kk | 1 << j: 1}); d[kk | 1 << i | 1 << j] -= 1; d[kk] -= 1
        row(d)
Aub = -sp.csr_matrix((V, (I, J)), shape=(r, N - 1)); bub = np.zeros(r)
eqs = []
allY = sum(Y(i) for i in range(ns))
for i, s in enumerate(sess): eqs.append(({Y(i): 1}, float(dem[s])))
eqs.append(({allY: 1}, float(sum(dem.values()))))
for Qm in range(1, (1 << k) - 1):
    Q = {c for c in range(k) if Qm >> c & 1}
    In = sum(Y(i) for i, (p, q) in enumerate(sess) if p in Q) | sum(X(i) for i, (p, q) in enumerate(arcs) if q in Q and p not in Q)
    Out = sum(X(i) for i, (p, q) in enumerate(arcs) if p in Q and q not in Q)
    Dec = sum(Y(i) for i, (p, q) in enumerate(sess) if q in Q)
    Bg = In | Out | Dec
    if Bg != In: eqs.append(({Bg: 1, In: -1} if In else {Bg: 1}, 0.0))
I, J, V = [], [], []
for rr, (d, _) in enumerate(eqs):
    for msk, cf in d.items(): I.append(rr); J.append(msk - 1); V.append(cf)
Aeq = sp.csr_matrix((V, (I, J)), shape=(len(eqs), N - 1)); beq = np.array([b for _, b in eqs])
c = np.zeros(N - 1)
for i in range(na): c[X(i) - 1] = 1.0
res = linprog(c, A_ub=Aub, b_ub=bub, A_eq=Aeq, b_eq=beq, bounds=(0, None), method='highs')
nA = int((g.D == 1).sum()); nU = int((g.D == m).sum()); cross = sum(arcsG.values())
print(f"my Shannon LP on the 3-cell contraction: {res.fun:.6f} (nc-gain: 80); per session after averaging {res.fun * nA / cross / nU:.5f} (nc-gain 1.07143)")
