"""cfgprobe.py (b5-twobranch): pure-Python probe of split3's configuration search (no numpy): how often do n1 array-1
branches adjacent to one A1 exist, and how often are all class-family images through the other branches at layer <= LIMY."""
import sys, random, time
from collections import Counter
from sp2 import *
from fastP import fast_find_P as fp

def H_of(c2, c3):
    return mul(inv(c2), mul(D, mul(inv(c3), D)))

seed, ntry, nlet, LIMY = [int(a) for a in sys.argv[1:5]]
rng = random.Random(seed); t0 = time.time()
M0 = image(D, L0)
N1s = [N for N in sorted(neighbors(L0)) if dist(M0, N) == 2]
st = Counter(); imgs = Counter()
for t in range(ntry):
    N1 = rng.choice(N1s)
    brs = [N for N in sorted(neighbors(N1)) if dist(L0, N) == 2 and dist(M0, N) == 3]
    anchors = [X for X in sorted(neighbors(M0)) if dist(X, L0) == 2 and dist(X, N1) == 3]
    P1, Q = rng.sample(anchors, 2)
    if dist(P1, Q) != 2:
        st['dPQ!=2'] += 1; continue
    byA = {}
    for N in brs:
        if dist(P1, N) == 4 and dist(Q, N) in (2, 3):
            for A in neighbors(N):
                if dist(P1, A) == 5 and dist(Q, A) <= 4:
                    byA.setdefault(A, []).append(N)
    good = [(A, Ns) for A, Ns in byA.items() if len(Ns) >= 2]
    st['A1 with >=2 branches: %d' % min(len(good), 3)] += 1
    if not good:
        continue
    A1, Ns = rng.choice(sorted(good, key=lambda x: str(x[0])))
    Na, Nb = rng.sample(Ns, 2)
    c2 = fp(rng, N1, M0)
    nok = 0
    for _ in range(nlet):
        c3a = fp(rng, image(D, image(c2, Na)), M0); c3b = fp(rng, image(D, image(c2, Nb)), M0)
        Ha, Hb = H_of(c2, c3a), H_of(c2, c3b)
        Ya = image(Hb, image(inv(Ha), A1)); Yb = image(Ha, image(inv(Hb), A1))
        imgs[((dist(P1, Ya), dist(Q, Ya)), (dist(P1, Yb), dist(Q, Yb)))] += 1
        if max(dist(P1, Ya), dist(Q, Ya), dist(P1, Yb), dist(Q, Yb)) <= LIMY or Ya == A1 or Yb == A1:
            nok += 1
    st['letter pairs ok'] += nok; st['letter pairs tried'] += nlet
    print("try %d: profiles %s, ok %d/%d (%.0fs)" % (t, [(dist(P1, N), dist(Q, N)) for N in (Na, Nb)], nok, nlet, time.time() - t0), flush=True)
print(dict(st))
print("image layers:", sorted(imgs.items(), key=lambda x: -x[1])[:25])
