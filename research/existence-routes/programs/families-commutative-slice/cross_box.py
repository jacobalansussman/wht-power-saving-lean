# cross-check of the session-2 tests against the coefficient solver at block 4, small W: list every deletion-passing
# multiset of (4, W, M) (LP box + box_enum2, as box_driver.py), classify by clique / even+odd pattern / cover(s=1,2),
# then run the LM solver (full order-2 system, 3 move orders) on (i) up to N multisets excluded ONLY by the new tests
# (a convergence would refute Theorem 10 or 12) and (ii) up to N multisets passing everything (control).
# usage: cross_box.py W M N restarts tmax
import sys, math, subprocess, time
import numpy as np, warnings
warnings.filterwarnings("ignore")
from scipy.optimize import linprog
from refl import quad_forms
from clique_test import clique_test
from even_test import epattern_test
from cover_test import cover_test
from refl2 import System2, random_x
m = 4
W, M, N, R, tmax = map(int, sys.argv[1:6])
vals, tr = quad_forms(m)
A = vals[tr == 1][:, 1:].astype(float); n = A.shape[1]
Aub = np.vstack([-A, np.ones((1, n))]); bub = np.concatenate([-np.ones(len(A)), [M / W + 1e-12]])
lo, hi = [], []
for z in range(n):
    cz = np.zeros(n); cz[z] = 1
    r1 = linprog(cz, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
    r2 = linprog(-cz, A_ub=Aub, b_ub=bub, bounds=[(0, None)] * n, method="highs")
    lo.append(max(0, math.ceil(W * r1.fun - 1e-7))); hi.append(math.floor(-W * r2.fun + 1e-7))
inp = f"{m} {W} {M}\n" + "\n".join(f"{a} {b}" for a, b in zip(lo, hi)) + "\n"
p = subprocess.run(["./box_enum2"], input=inp, capture_output=True, text=True)
vecs = [list(map(int, l.split())) for l in p.stdout.strip().splitlines() if l.strip()]
newly, allpass = [], []
ncl = 0
for cvec in vecs:
    S = [z + 1 for z in range(n) for _ in range(cvec[z])]
    if clique_test(m, S, W, verbose=False)[0] == "EXCLUDED": ncl += 1; continue
    e = epattern_test(m, S, W)[0]; c1 = cover_test(m, S, W, 1)[0]; c2 = cover_test(m, S, W, 2)[0]
    if "EXCLUDED" in (e, c1, c2): newly.append((S, e, c1, c2))
    else: allpass.append(S)
print(f"(4,{W},{M}): deletion-passing {len(vecs)}; clique-excluded {ncl}; excluded only by new tests {len(newly)}; pass all {len(allpass)}", flush=True)

def best_res(S):
    sy = System2(m, S, W); rng = np.random.default_rng(3); best = np.inf; t0 = time.time()
    for k in range(R):
        x, f, it = sy.solve(random_x(sy, rng), iters=500, tmax=tmax)
        best = min(best, f)
        if best < 1e-18 or time.time() - t0 > tmax: break
    return best

for tag, L in (("NEW-EXCLUDED (must stall)", newly[:N]), ("PASS-ALL (control)", [(S,) for S in allpass[:N]])):
    for item in L:
        S = item[0]; rng = np.random.default_rng(11)
        orders = [S, S[::-1], [int(v) for v in rng.permutation(S)]]
        bs = [best_res(o) for o in orders]
        print(f"  {tag}: S={S} {item[1:] if len(item) > 1 else ''} -> best residual {min(bs):.3e} {['%.1e' % b for b in bs]}", flush=True)
