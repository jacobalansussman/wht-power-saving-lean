"""junkabsorb.py -- mw2-absorb4 (resumed run): escape from the near miss's dilemma by a SHARED JUNK ABSORBER.

Dilemma (NOTE.md s3): with one route into each kernel frame through a filter phi, phi must be a unit (to reach K) and,
when the last orbit's junk can be cancelled only by the other arrays' private filters, a zero divisor.  Escape: let a
shared move J, built from the layer-1 generators (distance 3) and NOT downstream of the filters, put terms into the
junk frames of the last orbit, so that the junk can be cancelled without forcing phi to be a zero divisor.  Cost 3W+2
(W = 2: 8 = plain cost, a positive control of the mechanism; W = 6: 6 in 20).
Recipe: the near-miss limit (or any recipe) + J = delta (b u + tau(b u)) [kind S, A or T], u a layer-1 generator,
delta b u's term landing in a frame that holds a junk term of the last orbit.
Usage: junkabsorb.py recipe.pkl ntries starts [cap] [pairs]"""
import sys, pickle, time
from odd_twist import *

d = pickle.load(open(sys.argv[1], 'rb'))
ntries, starts = int(sys.argv[2]), int(sys.argv[3]); cap = int(sys.argv[4]) if len(sys.argv) > 4 else 4000
m, W, H = d['m'], d['W'], d['H']
rec = [(k, list(sm), dd) for (k, sm, dd) in d['recipe']]
tw0 = d.get('twist', {})
base = build_T(m, W, H, rec, d['skseed'], tw0)
G, R = base.G, base.R
last = base.moves[-1]
lastgens = [last['rep']] + last['images']
# junk frames: frames of the last orbit's terms that are not the kernel frame, with their array
junk = set()
for i in lastgens:
    for (l, g) in base.gens[i]:
        if G.dist(g) > 0:
            junk.add((G.frm(g), l))
print("last orbit generators", lastgens, "junk (frame, array) pairs", len(junk))
layer1 = base.moves[0]
l1gens = [layer1['rep']] + layer1['images']
rng = random.Random(43)
fd = fixed_dirs(R, rng)
opts = []
for dd in fd:
    for src in l1gens:
        for (l, g) in base.gens[src]:
            for b in G.progress_letters(g, rng, want=60, Dmat=dd, steps=1):
                g2 = G.intern(G.mul(dd, b, G.MATS[g]))
                if (G.frm(g2), l) in junk:
                    opts.append((dd, src, b))
print("(direction, source, letter) options landing in a junk frame:", len(opts), flush=True)
rng.shuffle(opts)
pairs = len(sys.argv) > 5 and sys.argv[5] == 'pairs'
if pairs:
    # two summands with the same direction landing in DIFFERENT junk frames: one shared move covering both
    tgt = {}
    for (dd, src, b) in opts:
        for (l, g) in base.gens[src]:
            g2 = G.intern(G.mul(dd, b, G.MATS[g]))
            if (G.frm(g2), l) in junk:
                tgt[(G.key(dd), src, G.key(b))] = (G.frm(g2), l)
    P2 = []
    for i in range(len(opts)):
        for j in range(i + 1, len(opts)):
            (d1_, s1, b1_), (d2_, s2, b2_) = opts[i], opts[j]
            if G.key(d1_) != G.key(d2_):
                continue
            if tgt[(G.key(d1_), s1, G.key(b1_))] == tgt[(G.key(d2_), s2, G.key(b2_))]:
                continue
            P2.append((d1_, [(s1, b1_), (s2, b2_)]))
    rng.shuffle(P2)
    items = P2
    print("two-summand absorbers covering two junk frames:", len(P2), flush=True)
else:
    items = [(dd, [(src, b)]) for (dd, src, b) in opts]
t0 = time.time(); res = []
for k, (dd, summ) in enumerate(items[:ntries]):
    kinds = [('S', None), ('A', None)] + [('T', bt) for (bt, o, btot) in beta_candidates_odd(G, R, dd, None)[:2]]
    for kind, bt in kinds:
        rec2 = rec + [(kind, list(summ), dd)]
        tw = dict(tw0)
        if kind == 'T':
            tw[len(rec)] = bt
        try:
            sk = build_T(m, W, H, rec2, d['skseed'], tw, R)
        except AssertionError:
            continue
        if sk.moves[-1]['images']:
            continue
        if kind == 'T' and not sk.moves[-1].get('absorbed'):
            continue
        if not sk.setup_final(cap):
            print(f"try {k} {kind}: closure cap"); continue
        prep_fast(sk)
        N = sk.nalpha; free_idx = [a for a in range(N) if a not in sk.fixed]; nf = len(free_idx)
        best = np.inf; bal = np.ones(N, dtype=complex)
        for st in range(starts):
            if nf == 0:
                r, _ = vp_eval(sk, np.ones(N, dtype=complex), [], want_J=False); best = float(np.sum(np.abs(r) ** 2)); break
            x0 = np.array([rng.gauss(0, 1) for _ in range(2 * nf)])
            x, c = vp_lm(sk, N, free_idx, x0, iters=120)
            if c < best:
                best = c; bal = np.ones(N, dtype=complex); bal[free_idx] = x[:nf] + 1j * x[nf:]
        if not res or best < min(r_[0] for r_ in res) - 1e-9:
            pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": d['skseed'], "alpha": bal, "res": best},
                        open(f"junkbest_{'pairs' if pairs else 'single'}.pkl", "wb"))
        dims = [(nr, nc) for (nr, nc, E) in sk.final]
        cost = sum(1 + len(mv['images']) for mv in sk.moves)
        res.append((best, k, kind))
        print(f"try {k} kind {kind}: cost {cost} closure {dims} full res^2 {best:.6f} ({time.time()-t0:.0f}s)", flush=True)
        if best < 1e-12:
            print("*** FULL HIT ***", flush=True)
            pickle.dump({"m": m, "W": W, "H": H, "recipe": rec2, "twist": tw, "skseed": d['skseed']}, open(f"HIT_junk_{k}_{kind}.pkl", "wb"))
res.sort(); print("best", res[:10])
