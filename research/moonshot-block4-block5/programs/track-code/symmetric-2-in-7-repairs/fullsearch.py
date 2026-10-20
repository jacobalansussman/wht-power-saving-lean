"""fullsearch.py -- mw2-absorb4: local search driven by the FULL exact-structure residual (not the frame-sum
relaxation, which mw2-absorb4 step 1 showed is passed at points the full problem rejects), on symmetric skeletons with
ONE ABSORBED shared move (twisted orbit sum, twist.py) and free orbits elsewhere.
Objective: closure per array (cap), variable-projection LM over the content coefficients (2 starts).  Mutations: replace
/ add / drop a summand of a free move, replace / add a summand of the shared move, change the twist.
Usage: fullsearch.py m W H shape menu seed nstarts steps cap [oddfirst]
shape uses mw-sym4 tokens; the 'L' move is converted to a twisted absorbed move 'T'."""
import sys, time, pickle, os
from twist import *

def full_res(sk, cap, rng, starts=2, iters=40):
    ok = sk.setup_final(cap)
    if not ok:
        return None, None
    prep_fast(sk)
    N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
    best = (np.inf, None)
    for st in range(starts):
        x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
        if nf == 0:
            r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False)
            return float(np.sum(np.abs(r) ** 2)), np.ones(N, dtype=complex)
        x, c = vp_lm(sk, N, free_idx, x0, iters=iters)
        if c < best[0]:
            al = np.ones(N, dtype=complex); al[free_idx] = x[:nf] + 1j * x[nf:]
            best = (c, al)
    return best

def main():
    m, W, H = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    shape = sys.argv[4].split(','); menu = int(sys.argv[5]); seed = int(sys.argv[6])
    nstarts, steps, cap = int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9])
    oddfirst = len(sys.argv) > 10 and sys.argv[10] == 'odd'
    t0 = time.time()
    for s in range(nstarts):
        sk0, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
        if sk0 is None:
            print(f"start {s}: build failed", flush=True); continue
        G, R = sk0.G, sk0.R
        rng = random.Random(seed * 13 + s)
        tsh = [t for t, (k, sm, d) in enumerate(rec) if k == 'L']
        if not tsh:
            print("no shared move in shape"); return
        tsh = tsh[0]
        Dm = rec[tsh][2]
        cands = beta_candidates(G, R, Dm, rng, nrand=1500)
        if oddfirst:
            odd = [c for c in cands if c[1] % 2 == 1]
            cands = odd if odd else cands
        if not cands:
            print(f"start {s}: no admissible twist", flush=True); continue
        orders = sorted(set(o for (_, o, _) in cands))
        rec = [(k, list(sm), d) for (k, sm, d) in rec]
        rec[tsh] = ('T', rec[tsh][1], Dm)
        tw = {tsh: rng.choice(cands)[0]}
        def build(rc, twd):
            sk = build_T(m, W, H, rc, seed * 100 + s, twd, R)
            if not sk.moves[tsh]['absorbed']:
                return None
            if sum(len(T) for T in sk.gens) > 300:
                return None
            return sk
        sk = build(rec, tw)
        tries = 0
        while sk is None and tries < 30:
            tw = {tsh: rng.choice(cands)[0]}; sk = build(rec, tw); tries += 1
        if sk is None:
            print(f"start {s}: no absorbed twist builds", flush=True); continue
        cur, cural = full_res(sk, cap, rng)
        if cur is None:
            cur = np.inf
        hist = [float("%.3g" % cur)]
        for it in range(steps):
            rec2 = [(k, list(sm), d) for (k, sm, d) in rec]; tw2 = dict(tw)
            op = rng.random()
            if op < 0.15:
                tw2[tsh] = rng.choice(cands)[0]
            elif op < 0.3:
                k_, sm_, d_ = rec2[tsh]
                src = rng.randrange(sk.moves[tsh]['rep'])
                if rng.random() < 0.5 and len(sm_) > 1:
                    sm_.pop(rng.randrange(len(sm_)))
                elif rng.random() < 0.5:
                    sm_[rng.randrange(len(sm_))] = (sm_[0][0], rand_free_b(G, rng, m))
                else:
                    sm_.append((src, rand_free_b(G, rng, m)))
            else:
                r2 = mutate(sk, rec2, rng)
                if r2 is None:
                    continue
                rec2 = r2
            try:
                sk2 = build(rec2, tw2)
            except Exception:
                continue
            if sk2 is None:
                continue
            if [1 + len(mv["images"]) for mv in sk2.moves] != [1 + len(mv["images"]) for mv in sk.moves]:
                continue
            r2, al2 = full_res(sk2, cap, rng)
            if r2 is None:
                continue
            if r2 < cur - 1e-6 or (abs(r2 - cur) < 1e-6 and rng.random() < 0.2):
                if r2 < cur - 1e-6:
                    hist.append(float("%.3g" % r2))
                sk, rec, tw, cur, cural = sk2, rec2, tw2, r2, al2
                pickle.dump({"m": m, "W": W, "H": H, "recipe": rec, "twist": tw, "skseed": seed * 100 + s, "alpha": cural,
                             "res": cur}, open(f"fs_{H}_{seed}_{s}.pkl", "wb"))
            if cur < 1e-10:
                break
        dims = [(nr, nc) for (nr, nc, E) in sk.final] if hasattr(sk, 'final') else None
        line = (f"start {s}: cost {cost} twist orders {orders} chosen order {order_of(G, G.mul(R.Hc[1][0], Dm, np.conj(R.Hc[1][0]).T, np.conj(Dm).T, Dm, tw[tsh], np.conj(Dm).T))} "
                f"history {hist} final closure {dims} ({time.time()-t0:.0f}s)")
        if cur < 1e-10:
            line += "  *** FULL HIT (numerical) ***"
        print(line, flush=True)

if __name__ == "__main__":
    main()
