#!/usr/bin/env python3
"""lp_triples.py (t2-block4): the multi-target rank-cover LP (mtlp_t.py's rows, built once) for random target sets
{I, T2, ..., TW} at block 4; only the right-hand side a_F changes.  Usage: lp_triples.py W count seed"""
import sys, itertools, numpy as np
sys.path.insert(0, '.')
from forms import setup, pair
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
W, cnt, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
m = 4; D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']
ALLT = sorted({D['enc_sym'](mm(h, tr(h))) for h in GL}); Ienc = D['enc_sym']([1 << i for i in range(m)])
nd = len(dirs); di = {x: i for i, x in enumerate(dirs)}
pairs = [(i, j) for i in range(nd) for j in range(i, nd)]; pk = {p: nd + k for k, p in enumerate(pairs)}
nv = nd + len(pairs)
ri, ci, vv, kind = [], [], [], []
r = 0
def add(row, which):
    global r
    for k, v in row.items(): ri.append(r); ci.append(k); vv.append(v)
    kind.append(which); r += 1
for F in range(1, nf):
    P = sorted(di[x] for x in supp[F]); row = {}
    for i in P: row[i] = 1; row[pk[(i, i)]] = -1
    for i, j in itertools.combinations(P, 2): row[pk[(i, j)]] = -2
    add(row, ('+', F))
for G in range(1, nf):
    P = sorted(di[x] for x in supp[G])
    add({i: -1 for i in P}, ('-', G))
    for i, j in itertools.combinations_with_replacement(P, 2):
        row = {i_: -1 for i_ in P}; row[pk[(i, j)]] = row.get(pk[(i, j)], 0) + 1
        add(row, ('-', G))
for (i, j), k in pk.items():
    add({k: 1, i: -1}, ('0',))
    if j != i: add({k: 1, j: -1}, ('0',))
A = coo_matrix((vv, (ri, ci)), shape=(r, nv)).tocsr()
cobj = np.zeros(nv); cobj[:nd] = 1
rng = np.random.default_rng(seed)
res_all = []
for c in range(cnt):
    T = [Ienc] + [ALLT[rng.integers(len(ALLT))] for _ in range(W - 1)]
    a = [sum(pair(F, S) for S in T) / len(T) for F in range(nf)]
    b = np.array([a[k[1]] if k[0] == '+' else (-a[k[1]] if k[0] == '-' else 0.0) for k in kind])
    res = linprog(cobj, A_ub=A, b_ub=b, bounds=[(0, None)] * nv, method='highs')
    res_all.append(res.fun)
    print('W=%d T=%s LP %.6f' % (W, T, res.fun), flush=True)
v = np.array(res_all); print('W=%d: %d sets, LP min %.4f median %.4f max %.4f; below %.4f (one move saved): %d' % (W, cnt, v.min(), np.median(v), v.max(), 4 - 1 / W, (v <= 4 - 1 / W + 1e-9).sum()))
