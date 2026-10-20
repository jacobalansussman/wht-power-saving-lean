# Independent rebuild of the rank-cover LP (ex2-families Theorems 17, 18) from the statement only.
# usage: my_rclp.py m schemes [deletion] [allodd]
#   schemes: comma list of 'cells' (direction-by-direction cells), 'ten' (halves + quarter rectangles)
# LP:  min sum_z c_z  s.t.
#   for every even (m-2)-flat F and every scheme partition of F x F into pieces:
#       sum_{x in F} c_x <= sum_pieces p_piece          (rank of nonsingular E_F <= sum of piece ranks)
#   p_(R,C) <= sum_R c, p_(R,C) <= sum_C c, p_(R,C) <= sum_G c - 1 for each odd G with R u C inside G
#   (one variable per unordered {R,C}: rank E[R,C] and rank E[C,R] obey the same bounds; p = max of both)
#   c, p >= 0.   [deletion: sum_G c >= 1 for each odd G]
# Then an exact rational lower bound from the HiGHS dual (repair + scaling, Fractions).
import sys, json, time
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
import geo

m = int(sys.argv[1])
schemes = sys.argv[2].split(",")
use_del = "deletion" in sys.argv[3:]
allodd = "allodd" in sys.argv[3:]
t0 = time.time()
ev, od = geo.even_odd_flats(m, m - 2)
print("block", m, "even (m-2)-flats", len(ev), "odd", len(od))
if allodd:
    od_supports = [S for (S, tr) in geo.all_form_supports(m) if tr == 1]
else:
    od_supports = od
pts = list(range(1, 1 << m))
nc = len(pts)

pieces = {}   # key -> index (key = frozenset({R, C}) with R, C frozensets)
def piece(R, C):
    key = frozenset([frozenset(R), frozenset(C)])
    if key not in pieces:
        pieces[key] = len(pieces)
    return pieces[key]

def hyperplanes(F):
    F = sorted(F)
    d = len(F).bit_length() - 1
    a = F[0]
    V = [x ^ a for x in F]
    out = set()
    # hyperplanes of the flat F: (d-1)-subflats
    from itertools import combinations
    for basis in combinations([v for v in V if v], d - 1):
        S = geo.span(basis)
        if len(S) != 1 << (d - 1):
            continue
        for b in F:
            out.add(frozenset(b ^ s for s in S))
    return [H for H in out if H <= set(F)]

rows_F = []   # (F, list of piece indices with multiplicity)
for F in ev:
    if "cells" in schemes:
        lst = []
        for x in F:
            for y in F:
                lst.append(piece({x}, {y}))
        rows_F.append((F, lst, "cells"))
    if "ten" in schemes:
        Hs = hyperplanes(F)
        seen = set()
        for H in Hs:
            Hp = frozenset(F) - H
            for J in Hs:
                if J == H or J == Hp:
                    continue
                K, L = H & J, H - J
                Kp, Lp = Hp & J, Hp - J
                part = frozenset([(H, H), (Hp, Hp), (K, Kp), (K, Lp), (L, Kp), (L, Lp),
                                  (Kp, K), (Lp, K), (Kp, L), (Lp, L)])
                if part in seen:
                    continue
                seen.add(part)
                # check exact tiling of F x F
                cells = {}
                for (R, C) in part:
                    for x in R:
                        for y in C:
                            cells[(x, y)] = cells.get((x, y), 0) + 1
                assert len(cells) == len(F) ** 2 and all(v == 1 for v in cells.values())
                lst = [piece(R, C) for (R, C) in part]
                rows_F.append((F, lst, "ten"))
print("F rows", len(rows_F), "pieces", len(pieces))

nv = nc + len(pieces)
A = []  # sparse rows: dict var->coef, rhs ; meaning sum coef*x <= rhs
labels = []
def cidx(z):
    return z - 1
for (F, lst, sch) in rows_F:
    r = {}
    for x in F:
        r[cidx(x)] = r.get(cidx(x), 0) + 1
    for p in lst:
        r[nc + p] = r.get(nc + p, 0) - 1
    A.append((r, 0)); labels.append(("F", sorted(F), sch))
od_sets = [frozenset(G) for G in od_supports]
for key, p in pieces.items():
    RC = list(key)
    R = RC[0]; C = RC[1] if len(RC) > 1 else RC[0]
    for side in (R, C):
        r = {nc + p: 1}
        for x in side:
            r[cidx(x)] = r.get(cidx(x), 0) - 1
        A.append((r, 0)); labels.append(("size", sorted(R), sorted(C)))
    U = R | C
    nG = 0
    for G in od_sets:
        if U <= G:
            nG += 1
            r = {nc + p: 1}
            for x in G:
                r[cidx(x)] = r.get(cidx(x), 0) - 1
            A.append((r, -1)); labels.append(("slack", sorted(R), sorted(C), sorted(G)))
    assert nG > 0, "piece with no odd form around it"
if use_del:
    for G in od_sets:
        r = {cidx(x): -1 for x in G}
        A.append((r, -1)); labels.append(("del", sorted(G)))
print("rows", len(A), "vars", nv, "build %.1fs" % (time.time() - t0))

Am = lil_matrix((len(A), nv))
b = np.zeros(len(A))
for i, (r, rhs) in enumerate(A):
    for k, v in r.items():
        Am[i, k] = v
    b[i] = rhs
cost = np.zeros(nv); cost[:nc] = 1
res = linprog(cost, A_ub=Am.tocsr(), b_ub=b, bounds=[(0, None)] * nv, method="highs")
print("status", res.status, res.message)
print("LP optimum (float) %.12f" % res.fun)
cz = res.x[:nc]
print("c at optimum (rounded):", [round(v, 6) for v in cz])
y = -res.ineqlin.marginals   # y >= 0 for <= rows in min problem
print("dual objective (float) %.12f" % (-(b @ y)))

# exact rational bound: for x >= 0 feasible, (A^T y)^T x <= b^T y.
def exact_bound(den):
    ys = {}
    for i, v in enumerate(y):
        if v > 1e-12:
            f = Fr(v).limit_denominator(den)
            if f > 0:
                ys[i] = f
    col = {}
    for i, f in ys.items():
        for k, v in A[i][0].items():
            col[k] = col.get(k, Fr(0)) + f * v
    # repair: piece columns must be >= 0; add size-rows (p - sum_R c <= 0, rhs 0) as needed
    size_row = {}
    for i, lab in enumerate(labels):
        if lab[0] == "size":
            p = [k for k, v in A[i][0].items() if k >= nc][0]
            size_row.setdefault(p, i)
    nrep = 0
    for k in range(nc, nv):
        v = col.get(k, Fr(0))
        if v < 0:
            i = size_row[k]
            ys[i] = ys.get(i, Fr(0)) - v
            for kk, vv in A[i][0].items():
                col[kk] = col.get(kk, Fr(0)) + (-v) * vv
            nrep += 1
    assert all(col.get(k, Fr(0)) >= 0 for k in range(nc, nv))
    bty = sum(f * A[i][1] for i, f in ys.items())
    alpha = max(-col.get(k, Fr(0)) for k in range(nc))
    assert alpha > 0
    return (-bty) / alpha, len(ys), nrep, ys

best = None
for den in (10**2, 10**3, 10**4, 10**6, 10**9, 10**12):
    bd, nrows, nrep, ys = exact_bound(den)
    print("limit_denominator %g: exact bound %s = %.12f  (rows used %d, repairs %d)" % (den, bd, float(bd), nrows, nrep))
    if best is None or bd > best[0]:
        best = (bd, den, ys)
print("BEST exact lower bound on M/W:", best[0], "= %.12f" % float(best[0]))
tag = "m%d_%s%s%s" % (m, "-".join(schemes), "_del" if use_del else "", "_allodd" if allodd else "")
with open("my_cert_%s.json" % tag, "w") as fh:
    json.dump({"bound": str(best[0]), "rows": [[labels[i], str(f)] for i, f in best[2].items()]}, fh)
print("time %.1fs" % (time.time() - t0))
