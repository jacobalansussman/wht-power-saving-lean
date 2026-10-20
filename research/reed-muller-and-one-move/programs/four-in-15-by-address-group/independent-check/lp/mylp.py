# mylp.py (v-w8-rmexact): own build of the multi-target cell LP (Proposition L) with scipy HiGHS; float values.
import sys, itertools, numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
sys_argv = sys.argv[:]; sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
def lp(H):
    n = len(H); T = [mul(h, tp(h)) for h in H]
    a = [sum(pairing(Q, S) for S in T) / n for Q in range(NF)]
    X = DIRS; xi = {x: i for i, x in enumerate(X)}
    cells = [(x, y) for i, x in enumerate(X) for y in X[i:]]; pi = {c: 15 + k for k, c in enumerate(cells)}
    nv = 15 + len(cells)
    rows = []; rhs = []
    def add(coef, b): rows.append(coef); rhs.append(b)   # sum coef . v <= b
    for G in range(1, NF):
        add({xi[x]: -1.0 for x in SUPP[G]}, -a[G])                       # deletion
        co = {}
        for x in SUPP[G]:
            for y in SUPP[G]:
                k = pi[(x, y) if xi[x] <= xi[y] else (y, x)]; co[k] = co.get(k, 0) - 1.0
        for x in SUPP[G]: co[xi[x]] = co.get(xi[x], 0) + 1.0
        add(co, a[G])                                                     # cover: slack(G) <= sum p
    for (x, y) in cells:
        k = pi[(x, y)]
        add({k: 1.0, xi[x]: -1.0}, 0.0)
        if y != x: add({k: 1.0, xi[y]: -1.0}, 0.0)
        for G in range(1, NF):
            if x in SUPP[G] and y in SUPP[G]:
                co = {k: 1.0}
                for z in SUPP[G]: co[xi[z]] = co.get(xi[z], 0) - 1.0
                add(co, -a[G])                                            # p <= slack(G)
    A = lil_matrix((len(rows), nv))
    for r, co in enumerate(rows):
        for k, v in co.items(): A[r, k] = v
    c = np.zeros(nv); c[:15] = 1.0
    res = linprog(c, A_ub=A.tocsr(), b_ub=np.array(rhs), bounds=[(0, None)] * nv, method='highs')
    return res.fun, len(rows)
mode = sys_argv[1]
if mode == 'gens':
    for s in sys_argv[2:]:
        H = cyclic(from_rows([int(v) for v in s.split(',')])); v, nr = lp(H)
        print('gen', s, '|H|', len(H), 'rows', nr, 'LP', '%.6f' % v, flush=True)
elif mode == 'singer5':
    def gmul(a, b):
        r = 0
        for i in range(4):
            if (b >> i) & 1: r ^= a << i
        for d in (6, 5, 4):
            if (r >> d) & 1: r ^= 0b10011 << (d - 4)
        return r
    for k in (5, 3, 1):
        ak = 1
        for _ in range(k): ak = gmul(ak, 2)
        cols = [gmul(ak, 1 << j) for j in range(4)]
        A_ = tuple(tuple((cols[j] >> i) & 1 for j in range(4)) for i in range(4))
        v, nr = lp(cyclic(A_)); print('Singer power', k, 'rows', nr, 'LP %.6f' % v, flush=True)
elif mode == 'classes3':
    GL = all_GL(); O4 = [g for g in GL if mul(g, tp(g)) == ID]
    inv = {g: next(h for h in GL if mul(g, h) == ID) for g in O4}
    subs = {}
    for g in GL:
        cg = cyclic(g)
        if len(cg) == 3: subs.setdefault(frozenset(cg), g)
    seen = set(); tot = 0; big = 0
    for s in sorted(subs, key=lambda s: sorted(s)):
        if s in seen: continue
        orbit = {frozenset(mul(mul(o, h), inv[o]) for h in s) for o in O4}; seen |= orbit
        v, nr = lp(sorted(s)); tot += len(orbit)
        if v >= 3.75 - 1e-9: big += len(orbit)
        print('class size %d gen %s LP %.6f' % (len(orbit), to_rows(subs[s]), v), flush=True)
    print('order-3 groups with LP >= 3.75: %d of %d' % (big, tot))
