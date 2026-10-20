"""rs2.py -- mw2-absorb4: replay mw-sym4's relax_search for ONE (seed, start) exactly (same code path, same rng use),
but checkpoint the current recipe (matrices included) and the best relaxation coefficients after every improvement.
Usage: rs2.py m W H shape menu seed start steps outprefix"""
import sys, time, pickle
from relax_search import *

def relax_res_x(sk, iters=25, starts=2):
    """identical to relax_search.relax_res (same rng draws) but also returns the best x"""
    G = sk.G
    fidx = {}
    for i, T in enumerate(sk.gens):
        for (l, g) in T:
            key = (G.frm(g), l)
            if key not in fidx:
                fidx[key] = len(fidx)
    Kkey = (G.frm(G.KID), 0)
    if Kkey not in fidx:
        return 1.0, None
    nrow = len(fidx)
    rowmap = [np.array([fidx[(G.frm(g), l)] for (l, g) in T]) for T in sk.gens]
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    def rr(al):
        cv = sk.coefs(al)
        A = np.zeros((nrow, len(sk.gens)), dtype=complex)
        for i in range(len(sk.gens)):
            np.add.at(A[:, i], rowmap[i], cv[i])
        t = np.zeros(nrow, dtype=complex); t[fidx[Kkey]] = 1
        beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
        return A @ beta - t
    def full(x):
        al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]; return al
    best = np.inf; bx = None
    rng = random.Random(len(fidx))
    sk._last_amax = 0.0
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        if nf == 0:
            return float(np.sum(np.abs(rr(np.ones(N, dtype=complex))) ** 2)), None
        x, c = lm(lambda x: rr(full(x)), x0, iters=iters)
        if c < best:
            best = c; bx = full(x); sk._last_amax = float(np.abs(full(x)).max())
    return best, bx

def main():
    m, W, H = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    shape = sys.argv[4].split(','); menu = int(sys.argv[5]); seed = int(sys.argv[6]); s = int(sys.argv[7])
    steps = int(sys.argv[8]); outp = sys.argv[9]
    t0 = time.time()
    sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
    if sk is None:
        print("build failed"); return
    orbits = [1 + len(mv["images"]) for mv in sk.moves]
    cur, curx = relax_res_x(sk)
    hist = [float("%.3g" % cur)]
    rng = random.Random(seed * 7 + s)
    def save(tag):
        pickle.dump({"m": m, "W": W, "H": H, "seed": seed, "start": s, "skseed": seed * 100 + s,
                     "recipe": rec, "alpha": curx, "res": cur, "hist": hist, "orbits": orbits},
                    open(f"{outp}_{tag}.pkl", "wb"))
    save("init")
    for it in range(steps):
        rec2 = mutate(sk, rec, rng)
        if rec2 is None:
            break
        try:
            sk2 = build_from_recipe(m, W, H, rec2, seed * 100 + s, sk.R)
        except Exception:
            continue
        if [1 + len(mv["images"]) for mv in sk2.moves] != orbits:
            continue
        if sum(len(T) for T in sk2.gens) > 400:
            continue
        r2, x2 = relax_res_x(sk2)
        if r2 < cur - 1e-6 or (abs(r2 - cur) < 1e-6 and rng.random() < 0.2):
            if r2 < cur - 1e-6:
                hist.append(float("%.3g" % r2))
            sk, rec, cur, curx = sk2, rec2, r2, x2
            if r2 < cur + 1e-12:
                save("cur")
                print(f"it {it}: res {cur:.3e} amax {np.abs(curx).max() if curx is not None else 0:.3g} ({time.time()-t0:.0f}s)", flush=True)
        if cur < 1e-10:
            break
    save("final")
    print(f"final: cost {cost} orbits {orbits} history {hist} ({time.time()-t0:.0f}s)", flush=True)

if __name__ == "__main__":
    main()
