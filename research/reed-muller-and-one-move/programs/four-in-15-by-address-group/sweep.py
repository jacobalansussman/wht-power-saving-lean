#!/usr/bin/env python3
"""sweep.py (w8-rmexact): the cell multi-target LP (mtlp.py) for EVERY cyclic subgroup of GL(4,2) of the given
orders (deduplicated; the LP depends only on the multiset of virtual targets h h^T).  The constraint matrix is built
once; only the right-hand side (the odd fractions a_F) changes.  Output: per subgroup the LP value (a lower bound on
M/W for every unit of the no-phase game with that address group), then a histogram, and how many are < 3.75."""
import sys, itertools, numpy as np
sys.path.insert(0, '.')
from forms import setup, pair
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
orders = [int(v) for v in sys.argv[1].split(',')]
part = sys.argv[2] if len(sys.argv) > 2 else '0/1'
pi_, pn_ = map(int, part.split('/'))
m = 4; D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']
I = [1 << i for i in range(m)]
def order_of(g):
    h = g; k = 1
    while h != I: h = mm(g, h); k += 1
    return k
subs = {}
for g in GL:
    k = order_of(g)
    if k not in orders: continue
    grp = [I]; h = g
    while h != I: grp.append(h); h = mm(g, h)
    key = frozenset(tuple(x) for x in grp)
    subs[key] = grp
subs = sorted(subs.values(), key=lambda gr: sorted(tuple(x) for x in gr))[pi_::pn_]
nd = len(dirs); di = {x: i for i, x in enumerate(dirs)}
pairs = [(i, j) for i in range(nd) for j in range(i, nd)]; pk = {p: nd + k for k, p in enumerate(pairs)}
nv = nd + len(pairs)
ri, ci, vv, rhs_of = [], [], [], []
r = 0
def addrow(row, which):
    global r
    for k, v in row.items(): ri.append(r); ci.append(k); vv.append(v)
    rhs_of.append(which); r += 1
for F in range(1, nf):
    P = sorted(di[x] for x in supp[F]); row = {}
    for i in P: row[i] = 1; row[pk[(i, i)]] = -1
    for i, j in itertools.combinations(P, 2): row[pk[(i, j)]] = -2
    addrow(row, ('+', F))
for G in range(1, nf):
    P = sorted(di[x] for x in supp[G])
    addrow({i: -1 for i in P}, ('-', G))
    for i, j in itertools.combinations_with_replacement(P, 2):
        row = {x: -1 for x in P}; row[pk[(i, j)]] = 1; addrow(row, ('-', G))
for (i, j), k in pk.items():
    addrow({k: 1, i: -1}, ('0', 0))
    if i != j: addrow({k: 1, j: -1}, ('0', 0))
A = coo_matrix((vv, (ri, ci)), shape=(r, nv)).tocsr()
cobj = np.zeros(nv); cobj[:nd] = 1
from collections import Counter
hist = Counter(); below = 0
for gr in subs:
    T = [D['enc_sym'](mm(h, tr(h))) for h in gr]; n = len(T)
    a = np.array([sum(pair(F, S) for S in T) / n for F in range(nf)])
    b = np.array([a[F] if s == '+' else (-a[F] if s == '-' else 0.0) for (s, F) in rhs_of])
    res = linprog(cobj, A_ub=A, b_ub=b, bounds=[(0, None)] * nv, method='highs')
    v = res.fun
    gen = gr[1]
    orth = all(D['enc_sym'](mm(h, tr(h))) == D['enc_sym'](I) for h in gr)
    print("order %2d gen %s targets %2d orthogonal %s: LP M/W >= %.6f" % (n, gen, len(set(T)), orth, v), flush=True)
    hist[round(v, 4)] += 1; below += v < 3.75 - 1e-9
print("histogram of LP values:", sorted(hist.items()))
print("subgroups with LP < 3.75 (4 in 15 not excluded by this LP): %d of %d" % (below, len(subs)))
