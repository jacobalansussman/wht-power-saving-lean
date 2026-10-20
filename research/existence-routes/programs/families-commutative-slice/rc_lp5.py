"""rc_lp5.py: block-5 rank-cover LP with plane pieces (the ten-piece scheme of Corollary 15 for every even 3-flat F and
every choice of halves and splitting direction), plus the singleton/cross-cell scheme.  W = 1.
Piece variables (bounds independent of F):
  y_x   : rank of E[x, x]                 <= c_x,             <= c(G) - 1  (odd G containing x)
  z_xy  : rank of E[x, y] and E[y, x]     <= c_x, <= c_y,     <= c(G) - 1  (odd G containing x, y)
  w_U   : rank of E[U, U], U a 2-flat     <= c(U),            <= c(G) - 1  (odd G containing U)
  v_KK' : rank of E[K, K'] and E[K', K], K, K' disjoint parallel lines with K u K' a 2-flat U
                                          <= c(K), <= c(K'),  <= c(G) - 1  (odd G containing U)
Rows: deletion; for every even 3-flat F: (36-piece) c(F) <= sum y + 2 sum z; (10-piece) for every pair of complementary
planes H, H' of F and every line direction d inside them: c(F) <= w_H + w_H' + 2 (v_{K K'} + v_{K L'} + v_{L K'} + v_{L L'}).
Odd forms used for the bounds: weight <= maxG.   usage: rc_lp5.py maxG
"""
import sys, itertools
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import csr_matrix
from refl import quad_forms

maxG = int(sys.argv[1]) if len(sys.argv) > 1 else 8
m = 5; n = 31
vals, tr = quad_forms(m)
V = vals[:, 1:].astype(np.int8); wt = V.sum(1)
oddS = [frozenset(np.nonzero(V[q])[0].tolist()) for q in np.nonzero(tr == 1)[0]]
boundS = [S for S in oddS if len(S) <= maxG]
def isflat(L):                      # directions are z-1 indices; flat test on the actual vectors
    pts = [z + 1 for z in L]
    s = set(pts); a = pts[0]
    D = {p ^ a for p in pts}
    return all((u ^ v) in D for u in D for v in D)
evenF = [frozenset(np.nonzero(V[q])[0].tolist()) for q in np.nonzero((tr == 0) & (wt == 8))[0]]
evenF = [F for F in evenF if isflat(sorted(F))]
var = {}
def vid(k):
    if k not in var: var[k] = len(var)
    return var[k]
for z in range(n): vid(("c", z))
R = []; LAB = []; CUR = ["?"]
def row(d, rhs): R.append((d, rhs)); LAB.append(CUR[0])
CUR[0] = "deletion"
for S in oddS: row({vid(("c", z)): -1.0 for z in S}, -1.0)
def bounded(key, U, sizes):
    """piece variable with bounds: <= c(T) for T in sizes, <= c(G) - 1 for odd G containing U (from boundS)"""
    if key in var: return var[key]
    j = vid(key); old = CUR[0]; CUR[0] = "bound:" + key[0]
    for T in sizes:
        d = {j: 1.0}
        for z in T: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) - 1.0
        row(d, 0.0)
    for S in boundS:
        if U <= S:
            d = {j: 1.0}
            for z in S: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) - 1.0
            row(d, -1.0)
    CUR[0] = old
    return j
nF = 0
for F in evenF:
    L = sorted(F); nF += 1
    # 36-piece
    CUR[0] = "R36"
    d = {}
    for z in L: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) + 1.0
    for x in L: j = bounded(("y", x), frozenset([x]), [[x]]); d[j] = d.get(j, 0.0) - 1.0
    for a, b in itertools.combinations(L, 2):
        j = bounded(("z", a, b), frozenset([a, b]), [[a], [b]]); d[j] = d.get(j, 0.0) - 2.0
    row(d, 0.0)
    # 10-piece: complementary planes H, H' of F (affine hyperplanes of the 3-flat) and a line direction
    pts = [z + 1 for z in L]
    a0 = pts[0]; Dir = sorted({p ^ a0 for p in pts} - {0})          # the 7 nonzero directions of F
    # planes of F: choose a 2-dim direction subspace (7 of them) -> its two cosets in F
    subs2 = set()
    for u, w in itertools.combinations(Dir, 2):
        subs2.add(frozenset([0, u, w, u ^ w]))
    for S2 in subs2:
        H = frozenset(p for p in pts if (p ^ a0) in S2); Hp = frozenset(pts) - H
        for t in sorted(S2 - {0}):                                   # line direction t inside S2
            def lines(Pl):
                out = []; seen = set()
                for p in sorted(Pl):
                    if p in seen: continue
                    ln = frozenset([p, p ^ t]); seen |= ln; out.append(ln)
                return out
            K, Lh = lines(H); Kp, Lp = lines(Hp)
            idx = lambda X: frozenset(z - 1 for z in X)
            CUR[0] = "R10"
            d = {}
            for z in L: d[vid(("c", z))] = d.get(vid(("c", z)), 0.0) + 1.0
            for P in (H, Hp):
                j = bounded(("w", idx(P)), idx(P), [sorted(idx(P))]); d[j] = d.get(j, 0.0) - 1.0
            for X in (K, Lh):
                for Y in (Kp, Lp):
                    key = ("v",) + tuple(sorted([tuple(sorted(idx(X))), tuple(sorted(idx(Y)))]))
                    j = bounded(key, idx(X) | idx(Y), [sorted(idx(X)), sorted(idx(Y))]); d[j] = d.get(j, 0.0) - 2.0
            row(d, 0.0)
nv = len(var)
ri, ci, vv, b = [], [], [], []
for i, (d, r) in enumerate(R):
    for j, cft in d.items(): ri.append(i); ci.append(j); vv.append(cft)
    b.append(r)
A = csr_matrix((vv, (ri, ci)), shape=(len(R), nv))
cost = np.zeros(nv); cost[:n] = 1
res = linprog(cost, A_ub=A, b_ub=np.array(b), bounds=[(0, None)] * n + [(None, None)] * (nv - n), method="highs")
print(f"block 5: even 3-flats {nF}; bound forms weight <= {maxG}: {len(boundS)}; variables {nv}; rows {len(R)}; "
      f"LP status {res.status}; value {res.fun if res.status == 0 else None} (plain 5, deletion LP 3.875)", flush=True)
if res.status == 0:
    c = res.x[:n]
    print("optimal c:", {z + 1: round(float(c[z]), 4) for z in range(n) if c[z] > 1e-9})
    from fractions import Fraction as Fr
    from exact_cert import certify
    tgt = Fr(res.fun).limit_denominator(1000)
    cr = certify(A, np.array(b), cost, n, res, tgt)
    print({k: (str(v) if k == "value" else v) for k, v in cr.items() if k not in ("sol", "sup")}, flush=True)
    import collections, json
    kinds = collections.Counter(LAB[i] for k, i in enumerate(cr["sup"]) if cr["sol"][k] != 0)
    print("certificate rows by kind:", dict(kinds))
    json.dump([{"row": i, "kind": LAB[i], "lambda": str(cr["sol"][k]), "coeffs": {str(j): v for j, v in R[i][0].items()}, "rhs": R[i][1]}
               for k, i in enumerate(cr["sup"]) if cr["sol"][k] != 0], open("rc_cert_b5.json", "w"))
    def jk(k):
        return [sorted(x) if isinstance(x, (set, frozenset)) else (list(x) if isinstance(x, tuple) else x) for x in k]
    json.dump({str(v): jk(k) for k, v in var.items()}, open("rc_vars_b5.json", "w"))
    if cr["certificate"]:
        print(f"EXACT: every commutative-slice unit at block 5 has M >= {cr['value']} W = {float(cr['value']):.6f} W")
