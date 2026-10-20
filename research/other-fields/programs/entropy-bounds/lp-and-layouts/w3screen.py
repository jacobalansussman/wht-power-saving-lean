# w3screen.py -- session 2. Exact screen for three-array coherent units at block m: is the cut value of a layout of
# three opposite pairs equal to 3m?
#  1. COMPATIBLE layouts (some two pairs i, j and ends u of i, w of j with d(u, w') = d(w, u') = m, primes = partners)
#     have cut value 3m by Lemma W3 (PROOF.md section 8): no computation needed.
#  2. Otherwise: pull back along psi(v) = (d(u1,v), d(u2,v), d(u3,v)) (u_i an end of pair i; all 8 choices tried):
#     the cut LP on the image graph H (at most (m+1)^3 vertices) is a lower bound on the cut value of Gamma_m, because
#     psi is non-expanding (adjacent frames have profiles at l-infinity distance <= 1, i.e. adjacent or equal in H) and
#     keeps every pair at distance m (coordinate i goes from 0 to m). Value 3m on H => value 3m on Gamma_m.
#  3. Layouts where every H falls short are written out for a full check on Gamma_m.
import sys, json, random, time, itertools
import numpy as np
import w3cut, w3compat, w3image

from fractions import Fraction
def exact_value(res):
    """exact lower bound from the dual of the final LP: lam >= 0 rationalised, scaled so that every edge load is <= 1
    (checked in exact integers/fractions). Valid for every coherent unit; since moves are integers, M >= ceil(value)."""
    val, lam, rows, rhs = res
    if lam is None: return Fraction(0)
    lam = [Fraction(max(0.0, float(x))).limit_denominator(10 ** 6) for x in lam]
    R = np.array(rows) > 0.5
    E = R.shape[1]; load = [Fraction(0)] * E
    for i, l in enumerate(lam):
        if l == 0: continue
        for e in np.nonzero(R[i])[0].tolist(): load[e] += l
    mx = max(load) if E else Fraction(1)
    sc = Fraction(1) if mx <= 1 else 1 / mx
    return sum(l * int(rhs[i]) for i, l in enumerate(lam)) * sc

THR_FRAC = False
def screen_one(g, P, orients=None, levels=('image3', 'image6', 'full')):
    """returns (exact value, level, orientation) for the first level that certifies M >= 3m (exact value > 3m - 1),
    or (best exact value, None, None)."""
    m = g.m; best = Fraction(0)
    thr = Fraction(3 * m) - Fraction(1, 10 ** 5) if THR_FRAC else 3 * m - 1
    if 'image3' in levels:
        for o in (orients or itertools.product((0, 1), repeat=3)):
            us = [P[i][o[i]] for i in range(3)]
            nH, EH, lab, keys = w3image.image_graph(g, us)
            v = exact_value(w3image.cut_value(nH, EH, [(lab[s], lab[t]) for s, t in P]))
            best = max(best, v)
            if v > thr: return v, 'image3', o
    if 'image6' in levels:
        nH, EH, lab, keys = w3image.image_graph(g, [x for p in P for x in p])
        v = exact_value(w3image.cut_value(nH, EH, [(lab[s], lab[t]) for s, t in P]))
        best = max(best, v)
        if v > thr: return v, 'image6', None
    if 'full' in levels:
        v = exact_value(w3image.cut_value(g.n, g.edges, P))
        best = max(best, v)
        if v > thr: return v, 'full', None
    return best, None, None

if __name__ == '__main__':
    m = int(sys.argv[1]); N = int(sys.argv[2]); seed = int(sys.argv[3]); tag = sys.argv[4] if len(sys.argv) > 4 else ''
    THR_FRAC = 'frac' in tag
    t0 = time.time(); g = w3cut.G(m); rng = random.Random(seed)
    print(f"Gamma_{m}: {g.n} frames ({time.time()-t0:.0f}s)", flush=True)
    ncomp = 0; nimg = 0; fails = []; rec = []
    for tr in range(N):
        P = w3cut.rand_3k2(g, rng)
        c = w3compat.any_compat(g.D, m, P)
        if c:
            ncomp += 1; continue
        v, lev, o = screen_one(g, P)
        if lev is not None:
            nimg += 1; rec.append(dict(P=P, level=lev, orient=o, value=str(v), value_float=float(v)))
            print(f"  certified at level {lev} (orientation {o}): {P} exact value {float(v):.6f} ({time.time()-t0:.0f}s)", flush=True)
        else:
            fails.append(dict(P=P, best=str(v))); print(f"  FAIL {P}: best exact value {float(v):.6f}", flush=True)
        print(f"  {tr+1} layouts: compatible {ncomp}, certified by pull-back {nimg}, open {len(fails)}", flush=True)
        json.dump(dict(m=m, seed=seed, N=tr + 1, compatible=ncomp, image=rec, fails=fails), open(f"w3screen_m{m}_s{seed}{tag}.json", "w"))
    print(f"m={m} seed={seed}: {N} random layouts of three opposite pairs with six distinct frames: compatible {ncomp}, "
          f"image-certified {nimg}, open {len(fails)} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(dict(m=m, seed=seed, N=N, compatible=ncomp, image=rec, fails=fails), open(f"w3screen_m{m}_s{seed}{tag}.json", "w"))
