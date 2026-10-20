"""repair.py -- mw-sym4: closure-driven repair of symmetric skeletons.
Observed failure of every symmetric skeleton (STATUS D): the best combination for [K]x_0 leaves junk rows that NO other
generator can cancel (the closure has only W columns).  Repair loop: solve (variable projection LM), take the junk rows
with the largest residual, and for each try to create a PARTNER: a summand b u_s added to the content of an earlier
free orbit-move t (s built before t) such that D_t b g_s lands on the junk term's frame (for some tau-image of the junk
term, so that some orbit image of move t carries it).  Neighbour frames of g_s are enumerated by sampling free letters.
Rebuild, re-solve; keep the change if res^2 drops.  The orbit structure (and so the cost) never changes.
Usage: repair.py m W H shape menu seed nskel rounds"""
import sys, time
from symsolve import *

def build_from_recipe(m, W, H, recipe, seed, R=None):
    sk = Skel(m, W, H, seed, R)
    for (kind, summ, Dm) in recipe:
        sk.add_move(kind, summ, Dm)
    return sk

def record_recipe(m, W, H, shape, menu, seed):
    """run build_skeleton while recording each add_move call"""
    rec = []
    orig = Skel.add_move
    def am(self, kind, summ, Dm_override=None):
        rec.append((kind, list(summ), Dm_override))
        return orig(self, kind, summ, Dm_override)
    Skel.add_move = am
    try:
        sk, cost = build_skeleton(m, W, H, shape, menu, seed)
    finally:
        Skel.add_move = orig
    return sk, cost, rec

def solve(sk, starts=2):
    if not sk.setup_final(6000):
        return None, None
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    best = (np.inf, None)
    for st in range(starts):
        x0 = np.array([sk.rng.gauss(0, 1) for _ in range(2 * nf)])
        if nf == 0:
            r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False)
            best = (float(np.sum(np.abs(r) ** 2)), x0); break
        x, c = vp_lm(sk, N, free_idx, x0)
        if c < best[0]:
            best = (c, x)
    al = np.ones(N, dtype=complex)
    if nf and best[1] is not None:
        al[free_idx] = best[1][:nf] + 1j * best[1][nf:]
    return best[0], al

def junk_rows(sk, al, top=6):
    """largest-residual non-target rows of array 0's system: list of (l, gid, weight)"""
    G = sk.G
    tok = []; gd = []
    for i, T in enumerate(sk.gens):
        d = {}
        for j, t in enumerate(T):
            d[t] = len(tok) + 1; tok.append((i, j))
        gd.append(d)
    rows, cols, ent = G.closure(gd, (0, G.KID), 6000)
    inv = {r: t for t, r in rows.items()}
    fl = flat(sk, al)
    nr, nc, R_, C_, F_ = sk.fin[0]
    A = np.zeros((nr, nc), dtype=complex); np.add.at(A, (R_, C_), fl[F_])
    t = np.zeros(nr, dtype=complex); t[0] = 1
    beta = np.linalg.lstsq(A, t, rcond=1e-10)[0]
    res = A @ beta - t
    order = np.argsort(-np.abs(res))
    out = []
    for ri in order[:top + 1]:
        if ri == 0 or abs(res[ri]) < 1e-9:
            continue
        l, gid = inv[ri]
        out.append((l, gid, abs(res[ri]) ** 2))
    return out[:top]

def neighbour_letters(sk, Dm, gid, rng, cache, samples=240):
    """frame_key(D b g) -> b for sampled free b"""
    key = (id(Dm), gid)
    if key in cache:
        return cache[key]
    G = sk.G
    out = {}
    g = G.MATS[gid]
    for _ in range(samples):
        b = rand_free_b(G, rng, sk.m)
        h = G.mul(Dm, b, g)
        fk = G.frame_key(h)
        if fk not in out:
            out[fk] = b
    cache[key] = out
    return out

def main():
    m, W, H, shape, menu, seed, nskel, rounds = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8])
    t0 = time.time()
    for s in range(nskel):
        sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
        if sk is None:
            print(f"skel {s}: build failed", flush=True); continue
        res, al = solve(sk)
        if res is None:
            print(f"skel {s}: closure cap", flush=True); continue
        hist = [round(res, 4)]
        rng = random.Random(seed + s)
        cache = {}
        for rd in range(rounds):
            G = sk.G; R = sk.R
            junk = junk_rows(sk, al)
            if not junk:
                break
            # rep generator index of each move in the recipe order
            rep_idx = [mv["rep"] for mv in sk.moves]
            improved = False
            tried = 0
            ncand_total = [0]
            for (l, gid, wgt) in junk:
                # tau-images of the junk term (any orbit image of a move may carry it)
                targets = {}
                for h in R.Hc:
                    l2, g2 = R.tau_term_cached(h, l, gid) if h is not R.Hc[0] else (l, gid)
                    targets[(l2, G.frm(g2))] = 1
                cands = []
                for t, (kind, summ, Dm) in enumerate(rec):
                    if kind not in ('F', 'P'):
                        continue
                    Dmat = G.D1 if Dm is None else Dm
                    for src in range(rep_idx[t]):
                        for (ls, gs) in sk.gens[src]:
                            nb = neighbour_letters(sk, Dmat, gs, rng, cache)
                            for fk, b in nb.items():
                                if (ls, fk) in targets:
                                    cands.append((t, src, b))
                rng.shuffle(cands)
                ncand_total[0] += len(cands)
                for (t, src, b) in cands[:4]:
                    tried += 1
                    rec2 = [(k, list(sm), d) for (k, sm, d) in rec]
                    rec2[t][1].append((src, b))
                    try:
                        sk2 = build_from_recipe(m, W, H, rec2, seed * 100 + s, sk.R)
                    except Exception as e:
                        continue
                    if [1 + len(mv["images"]) for mv in sk2.moves] != [1 + len(mv["images"]) for mv in sk.moves]:
                        continue
                    res2, al2 = solve(sk2)
                    if res2 is not None and res2 < res - 1e-6:
                        sk, rec, res, al = sk2, rec2, res2, al2
                        improved = True
                        break
                if improved:
                    break
            hist.append((round(float(res), 4), f"junk{len(junk)} cands{ncand_total[0]} tried{tried}"))
            if not improved or res < 1e-12:
                break
        print(f"skel {s}: cost {cost} orbits {[1 + len(mv['images']) for mv in sk.moves]} res^2 history {hist} ({time.time()-t0:.0f}s)", flush=True)
        if res < 1e-12:
            np.save(f"REPAIR_HIT_{m}_{W}_{H}_{seed}_{s}.npy", al)
            import pickle
            pickle.dump({"recipe_keys": [(k, [(src, G.key(b)) for (src, b) in sm], (G.key(d) if d is not None else None)) for (k, sm, d) in rec], "shape": shape}, open(f"REPAIR_HIT_{m}_{W}_{H}_{seed}_{s}.pkl", "wb"))
            print("HIT", flush=True)

if __name__ == "__main__":
    main()
