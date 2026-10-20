# stretch2.py -- Lemma W (two-class stretch): if the pairs of a layout split into two classes, each class with an
# orientation (low_i, high_i) such that d(low_i, high_j) = m for all i != j in the class, then the cut value is W*m and
# every coherent unit with that layout has M >= W*m. This program measures how often random layouts qualify.
import sys, random, itertools, json
import numpy as np
import w3cut
def stretchable(D, m, P, idx):
    """one 1-Lipschitz function stretches all pairs in idx iff some orientation has d(low_i, high_j) = m for i != j."""
    for o in itertools.product((0, 1), repeat=len(idx)):
        lo = [P[i][o[k]] for k, i in enumerate(idx)]; hi = [P[i][1 - o[k]] for k, i in enumerate(idx)]
        if all(D[lo[a], hi[b]] == m for a in range(len(idx)) for b in range(len(idx)) if a != b): return o
    return None
def two_class(D, m, P):
    W = len(P)
    for mask in range(1 << (W - 1)):          # pair 0 always in class A
        A = [i for i in range(W) if not (mask >> i & 1)]; B = [i for i in range(W) if mask >> i & 1]
        oa = stretchable(D, m, P, A)
        if oa is None: continue
        ob = stretchable(D, m, P, B) if B else ()
        if ob is not None: return (A, oa, B, ob)
    return None
def rand_layout(g, rng, W, distinct=True):
    while True:
        P = []; used = set()
        for _ in range(W):
            s = rng.randrange(g.n); t = rng.choice(g.opp[s]); P.append((s, t)); used |= {s, t}
        if not distinct or len(used) == 2 * W: return P
if __name__ == '__main__':
    m = int(sys.argv[1]); W = int(sys.argv[2]); N = int(sys.argv[3]); seed = int(sys.argv[4])
    g = w3cut.G(m); rng = random.Random(seed)
    hit = sum(1 for _ in range(N) if two_class(g.D, m, rand_layout(g, rng, W)))
    print(f"m={m} W={W}: {hit}/{N} random layouts (distinct frames) are two-class stretchable")
