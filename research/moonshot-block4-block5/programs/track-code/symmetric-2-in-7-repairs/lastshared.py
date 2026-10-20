"""lastshared.py -- mw2-absorb4 (resumed run): the near miss plus a SECOND generator at the kernel frame.

Junk lemma (NOTE.md s3): in the near-miss shape the kernel frame of array 0 is reached by one generator only (the last
orbit's g7), its array-1 part is junk that only array 1's private filter phi can cancel, and then phi must be a unit
(to reach K) and a zero divisor (to cancel the junk) at once -- for every filter and every twist, odd orders included.
v-ms-fano's sigma identity needs TWO different pieces reaching the target whose junk parts differ by sigma.  The
cheapest symmetric way to get a second piece at the kernel frame is a second SHARED move from the middle shared
generator to the kernel frame: delta (b u + gamma tau(b u)) (orbit sum, kind S / A, or twisted with an odd-order
gamma, kind T).  Cost 3W + 2 (W = 2: 8 = plain cost, a positive control of the mechanism; W = 6: 6 in 20, above the bar).
Usage: lastshared.py recipe.pkl ntries starts [cap]"""
import sys, pickle, time
from odd_twist import *

d = pickle.load(open(sys.argv[1], 'rb'))
ntries, starts = int(sys.argv[2]), int(sys.argv[3]); cap = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
m, W, H = d['m'], d['W'], d['H']
rec = [(k, list(sm), dd) for (k, sm, dd) in d['recipe']]
tw0 = d.get('twist', {})
base = build_T(m, W, H, rec, d['skseed'], tw0)
G, R = base.G, base.R
sh = [mv['rep'] for mv in base.moves if not mv['images']]
print("shared generators", sh, "terms", [[(l, G.dist(g)) for (l, g) in base.gens[i]] for i in sh])
gs = sh[-1]
rng = random.Random(41)
fd = fixed_dirs(R, rng)
print("fixed directions", len(fd))
lead = [g for (l, g) in base.gens[gs] if l == 0 and G.dist(g) == 1]
opts = []
for dd in fd:
    for g in lead:
        for b in G.progress_letters(g, rng, want=12, Dmat=dd, steps=1):
            if G.dist(G.intern(G.mul(dd, b, G.MATS[g]))) == 0:
                opts.append((dd, b))
print("(direction, letter) options taking the shared generator's array-0 term to K:", len(opts), flush=True)
rng.shuffle(opts)
t0 = time.time(); res = []
for k, (dd, b) in enumerate(opts[:ntries]):
    kinds = [('S', None), ('A', None)]
    oc = beta_candidates_odd(G, R, dd, None)
    kinds += [('T', bt) for (bt, o, btot) in oc[:3]]
    for kind, bt in kinds:
        rec2 = rec + [(kind, [(gs, b)], dd)]
        tw = dict(tw0)
        if kind == 'T':
            tw[len(rec)] = bt
        try:
            sk = build_T(m, W, H, rec2, d['skseed'], tw, R)
        except AssertionError:
            continue
        if sk.moves[-1]['images']:
            print(f"try {k} {kind}: last shared move is not tau-fixed"); continue
        if kind == 'T' and not sk.moves[-1].get('absorbed'):
            continue
        if not sk.setup_final(cap):
            print(f"try {k} {kind}: closure cap"); continue
        prep_fast(sk)
        N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
        best = np.inf
        for st in range(starts):
            if nf == 0:
                r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False); best = float(np.sum(np.abs(r) ** 2)); break
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, c = vp_lm(sk, N, free_idx, x0, iters=120)
            best = min(best, c)
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        cost = sum(1 + len(mv['images']) for mv in sk.moves)
        o = order_of(G, bt) if bt is not None else 0
        res.append((best, k, kind, o))
        print(f"try {k} kind {kind}{' twist order ' + str(o) if bt is not None else ''}: cost {cost} closure {dims} full res^2 {best:.6f} ({time.time()-t0:.0f}s)", flush=True)
        if best < 1e-12:
            print("*** FULL HIT ***", flush=True)
            pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": d['skseed']}, open(f"HIT_last_{k}_{kind}.pkl", "wb"))
res.sort(); print("best", res[:10])
