"""defilp.py (r3-classes): deficiency ILP, a strengthening of the clique test (CLQ).
Fact: L = E + K is unit lower triangular in time order and rank K <= W.  So for EVERY set X of moves
    |X| = rank L[X,X] <= rank E[X,X] + rank K[X,X] <= rank E[X,X] + W.
Upper bound on rank E[X,X] (order-free, coefficient-free): choose disjoint 'parts' P = D_G n X for odd forms G
(rank E[P,P] <= s_G, since E[P,P] is a principal block of E[D_G,D_G], of rank s_G by T10); every other cell of
X x X that is not a derived zero must be covered by a row or column of a vertex cover (term rank >= rank).  So
    def(X) := |X| - sum_parts s_G - |cover| <= rank K[X,X] <= W,
and def(X) >= W + 1 is a contradiction.  Singletons-only gives back CLQ.  The ILP maximises def(X).
usage: defilp.py tables.bin k1,..,kW "hex Z masks" [maxs]   (maxs: largest slack of a candidate part, default 8)"""
import struct, sys, numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
raw = open(sys.argv[1], 'rb').read(); nf = struct.unpack_from('<i', raw, 0)[0]
mask = struct.unpack_from('<%dI' % nf, raw, 4); wt = struct.unpack_from('<%di' % nf, raw, 4 + 4 * nf)
vid = struct.unpack_from('<%di' % nf, raw, 4 + 8 * nf); o = 4 + 12 * nf
nv = struct.unpack_from('<i', raw, o)[0]; pk = struct.unpack_from('<%dI' % nv, raw, o + 4)
ks = [int(x) for x in sys.argv[2].split(',')]; W = len(ks)
Z = [0] + [int(h, 16) for h in sys.argv[3].split()]
maxs = int(sys.argv[4]) if len(sys.argv) > 4 else 8
w = [sum(pk[v] >> k & 1 for k in ks) for v in range(nv)]
n = 31
def isz(x, y): return (Z[x] >> (y - 1)) & 1
cands = []
for q in range(1, nf):
    r = w[vid[q]]; s = wt[q] - r
    if r >= 2 and s <= maxs: cands.append((q, s, [x for x in range(1, 32) if mask[q] >> (x - 1) & 1]))
nG = len(cands)
# variable layout: x_v (0..30), rc_v (31..61), cc_v (62..92), p_G (93..)
NX, NR, NC, NPp = 0, 31, 62, 93
nvar = NPp + nG
c = np.zeros(nvar); c[NX:NX + 31] = -1; c[NR:NR + 31] = 1; c[NC:NC + 31] = 1
for g, (q, s, D) in enumerate(cands): c[NPp + g] = s
cover = {}
for g, (q, s, D) in enumerate(cands):
    for x in D:
        for y in D: cover.setdefault((x, y), []).append(g)
rows = []
# cell constraints: x_x + x_y - rc_x - cc_y - sum_{G: x,y in D_G} p_G <= 1   (non-zero cells, x != y)
# diagonal non-zero cells: x_x - rc_x - cc_x - sum_{G ni x} p_G <= 0
cells = [(x, y) for x in range(1, 32) for y in range(1, 32) if not isz(x, y)]
A = lil_matrix((len(cells) + 31, nvar)); lb = []; ub = []
for i, (x, y) in enumerate(cells):
    A[i, NX + x - 1] += 1
    if y != x: A[i, NX + y - 1] += 1
    A[i, NR + x - 1] -= 1; A[i, NC + y - 1] -= 1
    for g in cover.get((x, y), []): A[i, NPp + g] -= 1
    lb.append(-np.inf); ub.append(1 if y != x else 0)
# disjoint parts (as supports)
for v in range(1, 32):
    i = len(cells) + v - 1
    for g, (q, s, D) in enumerate(cands):
        if v in D: A[i, NPp + g] = 1
    lb.append(-np.inf); ub.append(1)
res = milp(c, constraints=LinearConstraint(A.tocsr(), lb, ub), integrality=np.ones(nvar), bounds=Bounds(0, 1),
           options={"time_limit": float(sys.argv[5]) if len(sys.argv) > 5 else 120})
val = -res.fun if res.x is not None else None
xs = res.x
X = [v for v in range(1, 32) if xs[NX + v - 1] > .5]
parts = [(cands[g][0], cands[g][1], [v for v in cands[g][2] if v in X]) for g in range(nG) if xs[NPp + g] > .5]
cov = ([v for v in range(1, 32) if xs[NR + v - 1] > .5], [v for v in range(1, 32) if xs[NC + v - 1] > .5])
print("targets %s candidates %d status %s deficiency %s (W = %d) |X|=%d parts=%s cover rows %s cols %s bound %s" % (
    sys.argv[2], nG, res.status, round(val) if val is not None else None, W, len(X), [(p[0], p[1], len(p[2])) for p in parts],
    cov[0], cov[1], getattr(res, 'mip_dual_bound', None)))
if val is not None and round(val) > W: print("EXCLUDED (deficiency %d > W = %d)" % (round(val), W))
