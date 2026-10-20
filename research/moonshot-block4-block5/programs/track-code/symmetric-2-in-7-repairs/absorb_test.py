"""absorb_test.py -- mw2-absorb4: replace the shared (orbit-sum) move of a recipe by TWISTED absorbed versions
(twist.py) for every admissible twist beta, and run the full exact-structure test on each.
Usage: absorb_test.py recipe.pkl move_index maxbeta starts [cap]"""
import sys, pickle, time
from twist import *

d = pickle.load(open(sys.argv[1], 'rb'))
t_sh = int(sys.argv[2]); maxb = int(sys.argv[3]); starts = int(sys.argv[4]); cap = int(sys.argv[5]) if len(sys.argv) > 5 else 4000
m, W, H = d.get('m', 4), d.get('W', 2), d.get('H', 'Z2swap')
rec = d['recipe']
base = build_T(m, W, H, rec, d['skseed'], {})
G, R = base.G, base.R
kind, summ, Dm = rec[t_sh]
print("shared move", t_sh, "kind", kind, "summands", len(summ))
rng = random.Random(5)
cands = beta_candidates(G, R, Dm, rng)
print("admissible twists:", len(cands), "orders of beta_tot:", sorted(set(o for (_, o, _) in cands)), flush=True)
rng.shuffle(cands)
cands.sort(key=lambda t: t[1] % 2 == 0)       # odd-order total factors first
results = []
t0 = time.time()
for k, (bt, o, btot) in enumerate(cands[:maxb]):
    rec2 = list(rec); rec2[t_sh] = ('T', summ, Dm)
    sk = build_T(m, W, H, rec2, d['skseed'], {t_sh: bt}, R)
    ab = sk.moves[t_sh]['absorbed']
    if not ab:
        print(f"twist {k}: not absorbed (orbit check failed)"); continue
    ok = sk.setup_final(cap)
    if not ok:
        print(f"twist {k}: order {o}: closure cap"); continue
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    best = np.inf
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        if nf == 0:
            r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False); best = float(np.sum(np.abs(r) ** 2)); break
        x, c = vp_lm(sk, N, free_idx, x0, iters=100)
        best = min(best, c)
    dims = [(nr, nc) for (nr, nc, E) in sk.final]
    results.append((best, o, k))
    print(f"twist {k}: order(beta_tot) {o}: closure {dims}: full res^2 {best:.4f} ({time.time()-t0:.0f}s)", flush=True)
    if best < 1e-12:
        pickle.dump({"recipe": rec2, "twist": {t_sh: bt}, "skseed": d['skseed'], "m": m, "W": W, "H": H}, open(f"HIT_absorb_{k}.pkl", "wb"))
        print("*** FULL HIT ***")
results.sort()
print("best:", results[:10])
