# Falsification run for Theorem 12 / Corollary 13 on U_5 with W = 7 (7 in 31, block 5).
# The proof of the exclusion uses ONLY these characters: the even form Q(z) = z_4 z_5 (support F = dirs 24..31) and the
# odd forms of slack 1 (the 240 odd 3-flats).  So the order-2 system restricted to these 241 forms must be infeasible.
# The LM solver (refl2 machinery) is run on that restricted system for several move orders and restarts; a residual
# reaching ~0 would refute the theorem.  Control: the same restricted system WITHOUT the even form F (no claim either way),
# and the solver on known units (calib.py, earlier) shows it does converge when solutions exist.
# usage: t_thm12_solver.py orders restarts iters tmax_per_run
import sys, time, json
import numpy as np, warnings
warnings.filterwarnings("ignore")
from refl import quad_forms
from refl2 import System2, random_x

nord, R, iters, tmax = map(int, sys.argv[1:5])
m, W = 5, 7
vals, tr = quad_forms(m)


def restricted(S, with_even=True):
    sy = System2(m, S, W, subset=None)
    sub = vals[:, S]; wt = sub.sum(1)
    keep_D, keep_phi = [], []
    for f in range(vals.shape[0]):
        D = np.nonzero(sub[f])[0]
        if len(D) == 0: continue
        isF = (tr[f] == 0 and set(S[i] for i in D) == set(range(24, 32)))
        if (tr[f] == 1 and wt[f] == 8) or (with_even and isF):
            keep_D.append(D); keep_phi.append(int(tr[f]))
    # rebuild the batched groups for the kept forms
    sy.Ds, sy.phis = keep_D, np.array(keep_phi)
    sy.nres = len(keep_D) * W * W
    iu = sy.iu; M = sy.M
    pidx = -np.ones((M, M), dtype=int); pidx[iu] = np.arange(len(iu[0]))
    sy.groups = []
    sizes = np.array([len(D) for D in keep_D])
    for d in sorted(set(sizes.tolist())):
        q = np.nonzero(sizes == d)[0]
        Dg = np.array([keep_D[i] for i in q]); base = q * W * W
        aa, bb = np.triu_indices(d, 1)
        hcols = 2 * sy.nA + pidx[Dg[:, aa], Dg[:, bb]]
        sy.groups.append((d, q, Dg, base, aa, bb, hcols))
    return sy


rng0 = np.random.default_rng(2026)
base = list(range(1, 32))
orders = [base, base[::-1]] + [[int(v) for v in rng0.permutation(base)] for _ in range(max(0, nord - 2))]
for with_even in (True, False):
    for oi, S in enumerate(orders[:nord]):
        sy = restricted(S, with_even)
        rng = np.random.default_rng(oi + 10)
        best = np.inf; t0 = time.time()
        for k in range(R):
            x, f, it = sy.solve(random_x(sy, rng), iters=iters, tmax=tmax)
            best = min(best, f)
            print(json.dumps({"with_even_F": with_even, "order": oi, "restart": k, "forms": len(sy.Ds),
                              "res2": float(f), "its": int(it), "t": round(time.time() - t0, 1)}), flush=True)
            if best < 1e-16: break
