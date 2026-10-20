"""box_enum.py: ALL move multisets for a commutative-slice target (m, W, M) that pass the deletion test, then the
tight-pair clique test (Theorem 1) on each.  Exhaustive: the LP {p >= 0, sum p <= M/W, sum_z p_z Q(z) >= 1 for odd Q}
bounds every count c_z = W p_z to [W lo_z, W hi_z]; every integer vector in that box with sum M is tested exactly.
usage: box_enum.py m W M [maxbox]"""
import sys, itertools, math, time
import numpy as np
from scipy.optimize import linprog
from refl import quad_forms
from clique_test import clique_test

m, W, M = map(int, sys.argv[1:4]); maxbox = float(sys.argv[4]) if len(sys.argv) > 4 else 5e7
vals, tr = quad_forms(m)
A = vals[tr == 1][:, 1:].astype(float)
n = A.shape[1]
Aub = np.vstack([-A, np.ones((1, n))]); bub = np.concatenate([-np.ones(len(A)), [M / W + 1e-12]])
lo, hi = [], []
for z in range(n):
    c = np.zeros(n); c[z] = 1
    r1 = linprog(c, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
    r2 = linprog(-c, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
    if r1.status != 0:
        print(f"m={m} W={W} M={M}: LP infeasible -> no multiset passes the deletion test"); sys.exit()
    lo.append(math.ceil(W * r1.fun - 1e-9)); hi.append(math.floor(-W * r2.fun + 1e-9))
box = 1
for a, b in zip(lo, hi): box *= (b - a + 1)
print(f"m={m} W={W} M={M}: count ranges {list(zip(lo, hi))}; box size {box:.3g}", flush=True)
if box > maxbox:
    print("box too large; stop"); sys.exit()
Ai = vals[tr == 1][:, 1:].astype(np.int32)
t0 = time.time(); npass = 0; nexcl = 0; survivors = []
# enumerate vectors with sum M inside the box (recursive, prune on remaining sum)
suf_lo = np.cumsum(lo[::-1])[::-1].tolist() + [0]; suf_hi = np.cumsum(hi[::-1])[::-1].tolist() + [0]
c = [0] * n
def rec(i, rem):
    global npass, nexcl
    if i == n:
        if rem != 0: return
        w = Ai @ np.array(c, dtype=np.int32)
        if w.min() < W: return
        npass += 1
        S = [z + 1 for z in range(n) for _ in range(c[z])]
        v = clique_test(m, S, W, verbose=False)
        if v[0] == "EXCLUDED": nexcl += 1
        else: survivors.append((S, v))
        return
    for x in range(lo[i], hi[i] + 1):
        r = rem - x
        if r < suf_lo[i + 1] or r > suf_hi[i + 1]: continue
        c[i] = x; rec(i + 1, r)
    c[i] = 0
rec(0, M)
print(f"  deletion-passing multisets {npass}; EXCLUDED by the clique test {nexcl}; survivors {len(survivors)}; "
      f"t={time.time()-t0:.0f}s", flush=True)
for S, v in survivors[:20]:
    print("  SURVIVOR", S, v, flush=True)
