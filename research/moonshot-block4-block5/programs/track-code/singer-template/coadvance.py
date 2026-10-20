"""coadvance.py (mw2-singer4): can ONE move advance two arrays' terms that a shared depth-2 move parked over N<x>?
Terms g_i = D1 rho_i D1 (both through N<x> = the frame of D1, then on the interval at depth 2, frames N<x,v_i>).
A move D1 b advances g along the interval iff the neighbour b^-1 N1 of L0 (N1 = D1^-1 L0) is mapped by g^-1 to a
forward neighbour of g's frame; so the 30 neighbour classes of b decide everything (exhaustive over b).
Prints, over sampled pairs (rho_1, rho_2): how many pairs admit a common advancing b, split by same frame / different
frames, and (Singer) for pairs related by eta-conjugation.  usage: coadvance.py samples seed"""
import sys, random
import numpy as np
from s4lib import S4, fast_frame_key

def main():
    n, seed = int(sys.argv[1]), int(sys.argv[2])
    S = S4(4, seed=seed); G = S.G; rng = random.Random(seed)
    # 30 neighbour classes of L0: representatives b with distinct frames of D1 b
    nb = {}
    for _ in range(4000):
        b = S.rand_free()
        k = fast_frame_key(G.mul(S.D1, b))
        nb.setdefault(k, b)
        if len(nb) >= 30:
            break
    reps = list(nb.values())
    # depth-2 interval terms through N<x>
    terms = []
    for _ in range(3000):
        r = S.rand_free(); g = G.intern(G.mul(S.D1, r, S.D1))
        if S.depth(g) == 2 and S.dK(g) == 2:
            terms.append(g)
        if len(terms) >= 200:
            break
    def adv(b, g):
        h = G.intern(G.mul(S.D1, b, G.MATS[g]))
        return S.depth(h) == 3 and S.dK(h) == 1
    advsets = {g: frozenset(i for i, b in enumerate(reps) if adv(b, g)) for g in terms}
    stats = {'same': [0, 0], 'diff': [0, 0]}
    for _ in range(n):
        g1, g2 = rng.sample(terms, 2)
        kind = 'same' if S.fk(g1) == S.fk(g2) else 'diff'
        stats[kind][0] += 1
        stats[kind][1] += bool(advsets[g1] & advsets[g2])
    sizes = sorted(set(len(v) for v in advsets.values()))
    print(f'neighbour classes of L0 found: {len(reps)}; depth-2 interval terms through N<x>: {len(terms)}; '
          f'advancing classes per term: {sizes}')
    print(f'pairs (total, with a common advancing move): same frame {stats["same"]}, different frames {stats["diff"]}')
    # Singer pairs: g and eta^k g eta^-k (if still through N<x> up to a left free factor, that is the HWc situation)
    sing = [0, 0]
    for g in terms[:60]:
        for k in range(1, 7):
            h = S.mulid(S.etap_id[k], g, S.inv(S.etap_id[k]))
            if S.depth(h) != 2 or S.dK(h) != 2:
                continue
            hs = frozenset(i for i, b in enumerate(reps) if adv(b, h))
            sing[0] += 1; sing[1] += bool(advsets[g] & hs)
    print(f'Singer-conjugate pairs (g, eta^k g eta^-k): total {sing[0]}, with a common advancing move {sing[1]}')

if __name__ == '__main__':
    main()
