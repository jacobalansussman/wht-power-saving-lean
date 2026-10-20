# mycert.py (v-w8-rmexact): exact LP lower bound from my own dual. Build the cell LP (as mylp.py), solve with HiGHS,
# rationalise the duals, and check weak duality EXACTLY with Fractions: y <= 0, A^T y <= c componentwise, bound = b.y.
import sys, numpy as np
from fractions import Fraction as Fr
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
sys_argv = sys.argv[:]; sys.argv = ['x']
exec(open('vchain.py').read().split("if __name__ == '__main__':")[0])
def build(H):
    n = len(H); T = [mul(h, tp(h)) for h in H]
    a = [Fr(sum(pairing(Q, S) for S in T), n) for Q in range(NF)]
    X = DIRS; xi = {x: i for i, x in enumerate(X)}
    cells = [(x, y) for i, x in enumerate(X) for y in X[i:]]; pi = {c: 15 + k for k, c in enumerate(cells)}
    rows = []; rhs = []
    def add(coef, b): rows.append(coef); rhs.append(b)
    for G in range(1, NF):
        add({xi[x]: -1 for x in SUPP[G]}, -a[G])
        co = {}
        for x in SUPP[G]:
            for y in SUPP[G]:
                k = pi[(x, y) if xi[x] <= xi[y] else (y, x)]; co[k] = co.get(k, 0) - 1
        for x in SUPP[G]: co[xi[x]] = co.get(xi[x], 0) + 1
        add(co, a[G])
    for (x, y) in cells:
        k = pi[(x, y)]
        add({k: 1, xi[x]: -1}, Fr(0))
        if y != x: add({k: 1, xi[y]: -1}, Fr(0))
        for G in range(1, NF):
            if x in SUPP[G] and y in SUPP[G]:
                co = {k: 1}
                for z in SUPP[G]: co[xi[z]] = co.get(xi[z], 0) - 1
                add(co, -a[G])
    return rows, rhs, 15 + len(cells)
for s in sys_argv[1:]:
    H = cyclic(from_rows([int(v) for v in s.split(',')]))
    rows, rhs, nv = build(H)
    A = lil_matrix((len(rows), nv))
    for r, co in enumerate(rows):
        for k, v in co.items(): A[r, k] = v
    c = np.zeros(nv); c[:15] = 1
    res = linprog(c, A_ub=A.tocsr(), b_ub=np.array([float(v) for v in rhs]), bounds=[(0, None)] * nv, method='highs')
    yf = res.ineqlin.marginals
    for den in (10 ** 4, 10 ** 6, 10 ** 8):
        y = [min(Fr(float(v)).limit_denominator(den), Fr(0)) if abs(v) > 1e-11 else Fr(0) for v in yf]
        col = [Fr(0)] * nv
        for r, co in enumerate(rows):
            if y[r]:
                for k, v in co.items(): col[k] += v * y[r]
        ok = all(col[k] <= (1 if k < 15 else 0) for k in range(nv))
        bound = sum(y[r] * rhs[r] for r in range(len(rows)) if y[r])
        if ok: break
    print('gen', s, '|H|', len(H), 'float LP %.6f' % res.fun, 'nonzero duals', sum(1 for v in y if v), 'dual feasible exactly:', ok,
          'exact bound', bound, '= %.6f' % float(bound), '>= 15/4:', bound >= Fr(15, 4), flush=True)
