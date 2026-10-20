"""lm_batch.py (rm-theory): Levenberg-Marquardt (mt_lm.MTSystem) on a file of small-lift survivors 'dirs|targets'
(t_small.py DUMP format; targets are 'gl:' indices). Each line: orders lex and rand<k>, `restarts` restarts each.
Evidence only: a residual below 1e-20 must go to exact checks before it counts.
Usage: python3 lm_batch.py m file restarts iters seed [stride offset]"""
import sys, time
import numpy as np
from mt_lm import Setup, MTSystem

m = int(sys.argv[1]); fn = sys.argv[2]; R = int(sys.argv[3]); iters = int(sys.argv[4]); seed = int(sys.argv[5])
stride = int(sys.argv[6]) if len(sys.argv) > 6 else 1; off = int(sys.argv[7]) if len(sys.argv) > 7 else 0
lines = [l.strip() for l in open(fn) if l.strip()]
best_all = []
for li, line in enumerate(lines):
    if li % stride != off: continue
    dirs, tg = line.split("|")
    for order in ([f"rand{seed + li}"] if __import__("os").environ.get("ORDERS") == "rand" else ["lex", f"rand{seed + li}"]):
        st = Setup(m, "gl:" + tg, "list:" + dirs, order, seed + li)
        et = st.exact_tests()
        if et["excluded"]:
            print(f"{line} {order}: EXCLUDED by mt_lm exact tests {et}", flush=True); continue
        sy = MTSystem(st); rng = np.random.default_rng(seed + 1000 * li + (0 if order == 'lex' else 1))
        best = None; t0 = time.time()
        for rs in range(R):
            x0 = (rng.standard_normal(sy.nvar) + 1j * rng.standard_normal(sy.nvar)) / np.sqrt(2)
            x0[2 * sy.nA:] *= 0.3
            x, f, it = sy.solve(x0, iters=iters)
            if best is None or f < best[1]: best = (x, f)
            if f < 1e-20:
                np.save(f"hit_m{m}_{li}_{order}_{rs}.npy", x); break
        print(f"{line} {order}: W={st.W} M={st.M} tight={et['tight']} best {best[1]:.3e} ({time.time()-t0:.0f} s)",
              flush=True)
        best_all.append(best[1])
print("SUMMARY: runs", len(best_all), "min", min(best_all) if best_all else None,
      "below 1e-20:", sum(1 for b in best_all if b < 1e-20))
