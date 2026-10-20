"""relax_search.py -- mw-sym4: local search in the FRAME-SUM RELAXATION (augmentation C[B] -> C), then the full test.
Start from a recorded symmetric skeleton (orbit structure fixed); moves: replace one summand's free letter by a random
free letter (sampled neighbour of the same source), add a summand (source = any earlier generator, letter sampled), or
drop a summand.  Objective: min over content coefficients (LM, few iterations) of the relaxation residual for every
array (by symmetry array 0 suffices).  Accept improvements (and sideways moves with prob 0.2).  When the relaxation
reaches 0, run the full exact-structure solver (symsolve variable projection) on that skeleton.
Usage: relax_search.py m W H shape menu seed nstarts steps"""
import sys, time
from repair import *

def relax_res(sk, iters=25, starts=2):
    G = sk.G
    fidx = {}
    for i, T in enumerate(sk.gens):
        for (l, g) in T:
            key = (G.frm(g), l)
            if key not in fidx:
                fidx[key] = len(fidx)
    Kkey = (G.frm(G.KID), 0)
    if Kkey not in fidx:
        return 1.0
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
    best = np.inf
    rng = random.Random(len(fidx))
    sk._last_amax = 0.0
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        if nf == 0:
            return float(np.sum(np.abs(rr(np.ones(N, dtype=complex))) ** 2))
        x, c = lm(lambda x: rr(full(x)), x0, iters=iters)
        if c < best:
            best = c; sk._last_amax = float(np.abs(full(x)).max())
    return best

def mutate(sk, rec, rng):
    rec2 = [(k, list(sm), d) for (k, sm, d) in rec]
    free_moves = [t for t, (k, sm, d) in enumerate(rec2) if k in ('F', 'P')]
    if not free_moves:
        return None
    t = rng.choice(free_moves)
    kind, summ, Dm = rec2[t]
    rep_t = sk.moves[t]["rep"]
    op = rng.random()
    G = sk.G
    if op < 0.5 and summ:                       # replace a letter
        k = rng.randrange(len(summ))
        src, b = summ[k]
        summ[k] = (src, rand_free_b(G, rng, sk.m))
    elif op < 0.85 or not summ:                 # add a summand
        src = rng.randrange(rep_t)
        summ.append((src, rand_free_b(G, rng, sk.m)))
    else:
        if len(summ) > 1:
            summ.pop(rng.randrange(len(summ)))
    rec2[t] = (kind, summ, Dm)
    return rec2

def main():
    m, W, H, shape, menu, seed, nstarts, steps = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
    t0 = time.time()
    for s in range(nstarts):
        sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
        if sk is None:
            print(f"start {s}: build failed", flush=True); continue
        orbits = [1 + len(mv["images"]) for mv in sk.moves]
        cur = relax_res(sk)
        hist = [float("%.3g" % cur)]
        rng = random.Random(seed * 7 + s)
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
            r2 = relax_res(sk2)
            if r2 < cur - 1e-6 or (abs(r2 - cur) < 1e-6 and rng.random() < 0.2):
                if r2 < cur - 1e-6:
                    hist.append(float("%.3g" % r2))
                sk, rec, cur = sk2, rec2, r2
            if cur < 1e-10:
                break
        line = f"start {s}: cost {cost} orbits {orbits} relaxation history {hist} max|alpha| {getattr(sk, '_last_amax', 0):.3g}"
        if cur < 1e-6:
            ok = sk.setup_final(6000)
            if ok:
                prep_fast(sk)
                N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
                bestf = np.inf
                for st in range(3):
                    x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
                    x, c = vp_lm(sk, N, free_idx, x0)
                    bestf = min(bestf, c)
                line += f" -> RELAXATION PASSED; full problem res^2 = {bestf:.3e}"
                if bestf < 1e-12:
                    line += "  *** FULL HIT ***"
                    import pickle
                    pickle.dump({"recipe": [(k, [(src, sk.G.key(b)) for (src, b) in sm], None if d is None else sk.G.key(d)) for (k, sm, d) in rec]}, open(f"RELAX_HIT_{m}_{W}_{H}_{seed}_{s}.pkl", "wb"))
            else:
                line += " -> relaxation passed, closure cap"
        print(line + f" ({time.time()-t0:.0f}s)", flush=True)

if __name__ == "__main__":
    G = None
    main()
