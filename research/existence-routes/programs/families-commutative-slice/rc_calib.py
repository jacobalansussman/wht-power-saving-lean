# calibration of the rank-cover inequalities (rc_lp2.py) on EXISTING units at block 4: for each unit (random disjoint
# unions of W representations of I_4, optionally with idle moves; plus the calibration list), fix c = counts / W and check
# that every rank-cover row can be satisfied, i.e. for every even 2-flat F:
#   c(F) <= sum_x min(c_x, min_{odd G ∋ x} (c(G) - 1)) + 2 sum_{x<y} min(c_x, c_y, min_{odd G ⊇ {x,y}} (c(G) - 1))
# (that is the R row with every y, z at its largest allowed value). A violation on an existing unit would refute the
# inequality (Theorem 17).  Also reports the minimum margin.
import sys, itertools
from fractions import Fraction as Fr
import numpy as np
from refl import quad_forms
m = 4
vals, tr = quad_forms(m)
V = vals[:, 1:]; n = 15
wt = V.sum(1)
odd = [q for q in range(len(V)) if tr[q] == 1]
oddsupp = sorted(set(frozenset(np.nonzero(V[q])[0].tolist()) for q in odd), key=len)
evenF = [frozenset(np.nonzero(V[f])[0].tolist()) for f in range(len(V))
         if tr[f] == 0 and wt[f] == 4 and np.bitwise_xor.reduce(np.nonzero(V[f])[0] + 1) == 0]
def margin(cnt, W):
    c = [Fr(cnt.get(z + 1, 0), W) for z in range(n)]
    worst = None
    for F in evenF:
        P = sorted(F)
        rhs = Fr(0)
        for x in P:
            sl = min(sum(c[z] for z in S) - 1 for S in oddsupp if x in S)
            rhs += min(c[x], sl)
        for a, b in itertools.combinations(P, 2):
            sl = min(sum(c[z] for z in S) - 1 for S in oddsupp if a in S and b in S)
            rhs += 2 * min(c[a], c[b], sl)
        mg = rhs - sum(c[x] for x in P)
        worst = mg if worst is None or mg < worst else worst
    return worst
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 1)
def cls(z): return tuple(((z >> j) & 1) * ((z >> k) & 1) for j in range(m) for k in range(j, m))
I = tuple(1 if j == k else 0 for j in range(m) for k in range(j, m))
pool = []
while len(pool) < 300:
    k = int(rng.integers(m, m + 5)); zs = [int(z) for z in rng.integers(1, 16, size=k)]
    s = np.zeros(len(I), dtype=int)
    for z in zs: s = (s + np.array(cls(z))) % 2
    if tuple(s.tolist()) == I: pool.append(zs)
bad = 0; mins = []
cases = [([1, 2, 4, 8], 1), ([7, 11, 13, 14], 1), ([1, 7, 2, 11, 4, 13, 8, 14], 2), ([1, 2, 4, 8] * 3, 3)]
for t in range(400):
    W = int(rng.integers(1, 6))
    S = sum((pool[i] for i in rng.integers(0, len(pool), size=W)), [])
    S += [int(z) for z in rng.integers(1, 16, size=int(rng.integers(0, 4)))]
    cases.append((S, W))
for S, W in cases:
    cnt = {}
    for z in S: cnt[z] = cnt.get(z, 0) + 1
    mg = margin(cnt, W); mins.append(mg)
    if mg < 0: bad += 1; print("VIOLATED on an existing unit:", W, sorted(S), mg)
print(f"{len(cases)} existing block-4 units: rank-cover rows violated {bad}; smallest margin {min(mins)}")
# and the excluded U_4 shapes must violate it
for r, W in [(1, 4), (3, 12), (5, 19), (3, 11)]:
    cnt = {z: r for z in range(1, 16)}
    print(f"U_4^({r}) with W={W}: margin {margin(cnt, W)}")
