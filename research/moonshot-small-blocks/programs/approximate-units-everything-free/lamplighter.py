"""lamplighter.py -- xw-border.  The free letters that commute with every fibrewise D1 sub-call (row-1-only letters)
contain the LAMPLIGHTER group Z2 wr Z: shift = a1 (running XOR on row 1), lamp = fX2 (X on row 1 of the first column).
Checks at f = 4..8: (i) the conjugates lamp_k = a1^k fX2 a1^-k are X-flips of row 1 by L^k e0 (L = I + shift);
(ii) lamps commute, square to 1, and lamp_0..lamp_K are independent (all 2^(K+1) products distinct) for K < f;
(iii) a1 has order 2^ceil(log2 f) (unbounded in f); (iv) every lamp_k and a1 commute with D1 for a random fibrewise
sub-call value (so the whole group is kept)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from keptlib import Fib

rng = np.random.default_rng(5)
for f in range(4, 9):
    F = Fib(f); G = F.G; n = F.n
    a1 = G.mono("a1")[f]; a1i = a1.inv(); lamp = G.mono("fX2")[f]
    def conj(g, k):
        h = g
        for _ in range(abs(k)):
            h = (a1i.then(h).then(a1)) if k > 0 else (a1.then(h).then(a1i))   # then = apply self first
        return h
    # order of a1
    o, h = 1, a1
    while not np.array_equal(h.sig, np.arange(n)):
        h = h.then(a1); o += 1
    K = min(f - 1, 5)
    lamps = [conj(lamp, k) for k in range(K + 1)]
    # which row-1 flip vector is lamp_k?  read off from the image of index 0
    vecs = []
    for L_ in lamps:
        x0 = int(L_.sig[0]); vecs.append([(x0 >> (2 * j + 1)) & 1 for j in range(f)])
    # products distinct?
    seen = set()
    for mask in range(1 << (K + 1)):
        g = G.L[f].ident()
        for k in range(K + 1):
            if mask >> k & 1: g = g.then(lamps[k])
        seen.add(g.sig.tobytes())
    comm = all(np.array_equal(lamps[i].then(lamps[j]).sig, lamps[j].then(lamps[i]).sig) for i in range(K + 1) for j in range(K + 1))
    sq = all(np.array_equal(l.then(l).sig, np.arange(n)) for l in lamps)
    Y = F.rand_Y(rng); v = rng.standard_normal((n, 2)) + 1j * rng.standard_normal((n, 2))
    kept = max(np.linalg.norm(F.D(1, G.apply_mono(l, v), Y) - G.apply_mono(l, F.D(1, v, Y))) for l in lamps + [a1])
    print("f=%d: order(a1) = %d; lamp_k flip vectors (row 1, columns 0..f-1) k=0..%d: %s; commute %s, square 1 %s, "
          "%d distinct products of %d; max |[D1(Y), g]| %.1e" % (f, o, K, vecs, comm, sq, len(seen), 1 << (K + 1), kept), flush=True)
