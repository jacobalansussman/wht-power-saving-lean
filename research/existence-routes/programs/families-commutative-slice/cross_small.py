# cross-check of the new exact tests against the coefficient solver on SMALL cases (block 3):
# every multiset of M moves (M = 3W .. 3W + 3) passing the deletion test is classified by clique / even / odd / triple
# tests; the LM solver (refl2.System2, full order-2 system, several restarts) is run on (i) every multiset the NEW tests
# exclude but the clique test passes (it must stall: a convergence would refute Theorem 10 or 12) and (ii) up to 6
# multisets every test passes (control: the solver should find solutions for real units).
# usage: cross_small.py m W Mmax restarts
import sys, itertools, numpy as np, warnings, time
warnings.filterwarnings("ignore")
from refl import quad_forms
from clique_test import clique_test
from even_test import epattern_test
from triple_test import triple_test
from refl2 import System2, random_x
m, W, Mmax, R = map(int, sys.argv[1:5])
vals, tr = quad_forms(m); odd = tr == 1
n = (1 << m) - 1
newly = []; allpass = []
for M in range(m * W, Mmax + 1):
    cnt_new = cnt_pass = cnt_del = cnt_cl = 0
    for combo in itertools.combinations_with_replacement(range(1, n + 1), M):
        S = list(combo)
        cnt = np.bincount(S, minlength=1 << m)
        wt = (vals * cnt[None, :]).sum(1)
        if wt[odd].min() < W: continue
        cnt_del += 1
        c = clique_test(m, S, W, verbose=False)
        if c[0] == "EXCLUDED": cnt_cl += 1; continue
        e = epattern_test(m, S, W); t = triple_test(m, S, W)
        if e[0] == "EXCLUDED" or t[0] == "EXCLUDED":
            cnt_new += 1; newly.append((M, S, e[0], t[0]))
        else:
            cnt_pass += 1; allpass.append((M, S))
    print(f"m={m} W={W} M={M}: deletion-passing {cnt_del}, clique-excluded {cnt_cl}, newly excluded {cnt_new}, all pass {cnt_pass}", flush=True)

def solve_best(S, R, seed=1, iters=400, tmax=60):
    sy = System2(m, S, W); rng = np.random.default_rng(seed); best = np.inf; t0 = time.time()
    for k in range(R):
        x, f, it = sy.solve(random_x(sy, rng), iters=iters, tmax=tmax)
        best = min(best, f)
        if best < 1e-18 or time.time() - t0 > tmax: break
    return best

print("--- solver on NEWLY EXCLUDED multisets (must stall) ---", flush=True)
for (M, S, ev, tv) in newly[:12]:
    # try a few move orders: sorted, reversed, random
    rng = np.random.default_rng(7)
    orders = [S, S[::-1], list(rng.permutation(S))]
    bests = [solve_best(list(map(int, o)), R) for o in orders]
    print(f"  M={M} S={S} even:{ev} triple:{tv} -> solver best residual over 3 orders {min(bests):.3e} {['%.1e' % b for b in bests]}", flush=True)
print("--- solver on multisets that pass every test (control) ---", flush=True)
for (M, S) in allpass[:6]:
    rng = np.random.default_rng(7)
    orders = [S, S[::-1], list(rng.permutation(S))]
    bests = [solve_best(list(map(int, o)), R) for o in orders]
    print(f"  M={M} S={S} -> solver best residual over 3 orders {min(bests):.3e}", flush=True)
