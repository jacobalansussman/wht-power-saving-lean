"""full_test.py -- mw2-absorb4: full exact-structure test of a recipe (pickled by limit_skel.py / rs2.py):
closure over frames for every array (bm.Game.closure), then variable-projection LM over the free content coefficients
(dense; only for closures up to a few thousand columns), many starts; reports min residual^2 per start, and the relaxation
residual at the same coefficients.  Usage: full_test.py recipe.pkl starts [cap] [fix a=v ...]"""
import sys, pickle, time
from relax_search import *

d = pickle.load(open(sys.argv[1], 'rb'))
m, W, H = d.get('m', 4), d.get('W', 2), d.get('H', 'Z2swap')
starts = int(sys.argv[2]); cap = int(sys.argv[3]) if len(sys.argv) > 3 else 20000
sk = build_from_recipe(m, W, H, d['recipe'], d['skseed'])
t0 = time.time()
ok = sk.setup_final(cap)
print("closure dims", [(nr, nc) for (nr, nc, E) in sk.final] if ok else "cap", f"({time.time()-t0:.0f}s)", flush=True)
if not ok:
    sys.exit()
prep_fast(sk)
N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
print("free alphas", free_idx)
rng = random.Random(7)
res = []
for st in range(starts):
    x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
    x, c = vp_lm(sk, N, free_idx, x0, iters=150)
    al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]
    res.append((c, al))
    print(f"start {st}: full res^2 {c:.4e} alpha {np.round(al[free_idx], 4)}", flush=True)
res.sort(key=lambda t: t[0])
pickle.dump({"best": res[:5]}, open(sys.argv[1].replace('.pkl', '_full.pkl'), 'wb'))
