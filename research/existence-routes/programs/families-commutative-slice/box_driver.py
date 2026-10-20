"""box_driver.py: LP box (as box_enum.py) -> pruned exhaustive C search box_enum2 -> clique test on every
deletion-passing multiset.  usage: box_driver.py m W M"""
import sys, math, subprocess, time
import numpy as np
from scipy.optimize import linprog
from refl import quad_forms
from clique_test import clique_test
m, W, M = map(int, sys.argv[1:4])
vals, tr = quad_forms(m)
A = vals[tr == 1][:, 1:].astype(float); n = A.shape[1]
Aub = np.vstack([-A, np.ones((1, n))]); bub = np.concatenate([-np.ones(len(A)), [M / W + 1e-12]])
lo, hi = [], []
for z in range(n):
    cz = np.zeros(n); cz[z] = 1
    r1 = linprog(cz, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
    if r1.status != 0:
        print(f"m={m} W={W} M={M}: LP infeasible -> nothing passes the deletion test", flush=True); sys.exit()
    r2 = linprog(-cz, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
    lo.append(max(0, math.ceil(W * r1.fun - 1e-7))); hi.append(math.floor(-W * r2.fun + 1e-7))
if any(a > b for a, b in zip(lo, hi)):
    print(f"m={m} W={W} M={M}: empty LP box -> nothing passes the deletion test", flush=True); sys.exit()
t0 = time.time()
inp = f"{m} {W} {M}\n" + "\n".join(f"{a} {b}" for a, b in zip(lo, hi)) + "\n"
p = subprocess.run(["./box_enum2"], input=inp, capture_output=True, text=True)
vecs = [list(map(int, l.split())) for l in p.stdout.strip().splitlines() if l.strip()]
nex = 0; surv = []
for cvec in vecs:
    S = [z + 1 for z in range(n) for _ in range(cvec[z])]
    v = clique_test(m, S, W, verbose=False)
    if v[0] == "EXCLUDED": nex += 1
    else: surv.append((S, v))
print(f"m={m} W={W} M={M}: box {list(zip(lo, hi))}; {p.stderr.strip()}; EXCLUDED by clique {nex}; survivors {len(surv)}; "
      f"t={time.time()-t0:.0f}s", flush=True)
for S, v in surv[:10]: print("  SURVIVOR", S, v, flush=True)
