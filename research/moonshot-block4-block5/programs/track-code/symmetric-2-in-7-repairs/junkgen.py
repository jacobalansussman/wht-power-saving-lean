"""junkgen.py -- mw2-absorb4 (resumed run): the 3W+2 "shared junk absorber" shape on fresh symmetric skeletons.

Base: record_recipe(m, W, H, shape, menu, seed) (mw-sym4's builder; default shape P:inp,P:prev,L:prev,P:prev = 3W+1).
Absorber: ONE more tau-fixed move J = delta * orbit-sum of summands taken from the layer-1 generators (distance 3, not
downstream of the filters), each summand landing in a junk frame of the last orbit (a frame holding a non-kernel term).
Modes: 'single' (one summand), 'cover' (greedy: one summand per distinct junk frame reachable with the same delta, up to
4 summands).  Kinds: S (trivial character), A/B (characters chi, chi^2), T (order-3 twist, odd_twist.py).
Each skeleton: closure + variable-projection LM over all content coefficients (full exact-structure residual).
Usage: junkgen.py m W H shape menu seed0 nbases ntries starts cap mode"""
import sys, pickle, time
from odd_twist import *

def absorb(rec, tw0, m, W, H, skseed, R, ntries, starts, cap, mode, rng, tag):
    base = build_T(m, W, H, rec, skseed, tw0, R)
    G = base.G; R = base.R
    last = base.moves[-1]
    lastgens = [last['rep']] + last['images']
    junk = set()
    for i in lastgens:
        for (l, g) in base.gens[i]:
            if G.dist(g) > 0:
                junk.add((G.frm(g), l))
    l1 = base.moves[0]; l1gens = [l1['rep']] + l1['images']
    if not hasattr(R, '_fd'):
        R._fd = fixed_dirs(R, rng)
    opts = {}
    for k_, dd in enumerate(R._fd):
        for src in l1gens:
            for (l, g) in base.gens[src]:
                for b in G.progress_letters(g, rng, want=40, Dmat=dd, steps=1):
                    g2 = G.intern(G.mul(dd, b, G.MATS[g]))
                    key = (G.frm(g2), l)
                    if key in junk:
                        opts.setdefault(k_, {}).setdefault(key, (src, b))
    items = []
    for k_, byf in opts.items():
        dd = R._fd[k_]
        if mode == 'single':
            for key, sb in byf.items():
                items.append((dd, [sb]))
        else:
            items.append((dd, list(byf.values())[:4]))
    rng.shuffle(items)
    out = []
    for k, (dd, summ) in enumerate(items[:ntries]):
        kinds = [('S', None), ('A', None)] + ([('B', None)] if W >= 3 else [])
        kinds += [('T', bt) for (bt, o, btot) in beta_candidates_odd(G, R, dd, None)[:1]]
        for kind, bt in kinds:
            rec2 = rec + [(kind, list(summ), dd)]
            tw = dict(tw0)
            if kind == 'T':
                tw[len(rec)] = bt
            try:
                sk = build_T(m, W, H, rec2, skseed, tw, R)
            except AssertionError:
                continue
            if sk.moves[-1]['images'] or (kind == 'T' and not sk.moves[-1].get('absorbed')):
                continue
            if not sk.setup_final(cap):
                continue
            prep_fast(sk)
            N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
            best = np.inf
            for st in range(starts):
                if nf == 0:
                    r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False); best = float(np.sum(np.abs(r) ** 2)); break
                x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
                x, c = vp_lm(sk, N, free_idx, x0, iters=100)
                best = min(best, c)
            cost = sum(1 + len(mv['images']) for mv in sk.moves)
            out.append((best, kind, len(summ), cost))
            if best < 1e-12:
                print(f"*** FULL HIT *** {tag} try {k} kind {kind} cost {cost}", flush=True)
                pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": skseed}, open(f"HIT_jg_{tag}_{k}_{kind}.pkl", "wb"))
    return len(junk), len(items), out

def main():
    m, W, H, shape, menu = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5])
    seed0, nb, ntries, starts, cap, mode = (int(sys.argv[6]), int(sys.argv[7]), int(sys.argv[8]), int(sys.argv[9]),
                                            int(sys.argv[10]), sys.argv[11])
    t0 = time.time(); R = None; allbest = []
    for s in range(nb):
        seed = seed0 * 100 + s
        sk0, cost, rec = record_recipe(m, W, H, shape, menu, seed)
        if sk0 is None:
            print(f"base {s}: build failed", flush=True); continue
        R = sk0.R
        # base floor
        bres = np.inf
        if sk0.setup_final(cap):
            prep_fast(sk0)
            N = sk0.nalpha; fi = [a for a in range(N) if a not in sk0.fixed]; nf = len(fi)
            rng0 = random.Random(seed)
            for st in range(2):
                x, c = vp_lm(sk0, N, fi, np.array([rng0.gauss(0, 1) for _ in range(2 * nf)]), iters=100) if nf else (None, np.inf)
                bres = min(bres, c)
        rng = random.Random(seed + 7)
        nj, nit, out = absorb(rec, {}, m, W, H, seed, R, ntries, starts, cap, mode, rng, f"{H}_{seed}")
        out.sort()
        allbest += out[:1]
        print(f"base {s}: cost {cost}, base full floor {bres:.4f}; junk frames {nj}; absorbers {nit}; "
              f"best with absorber {[(round(b, 4), k, ns, c) for (b, k, ns, c) in out[:4]]} ({time.time()-t0:.0f}s)", flush=True)
    allbest.sort()
    print("overall best", allbest[:8])

if __name__ == "__main__":
    main()
