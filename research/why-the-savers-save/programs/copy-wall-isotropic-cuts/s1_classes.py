"""s1_classes.py -- for every t.i. cut g (dim 1..4): |C_g|, omega(C_g), octads inside, octads met, and a greedy
partition of C_g into orthogonal cliques (k cliques => |T| <= k * omega(T) for every T inside C_g)."""
from base import *
from collections import Counter
import time
t0 = time.time()
O = octads(); print('octads in P112:', len(O))
Omask = [sum(1 << i for i in o) for o in O]
TI = ti_subspaces(4)
for d in TI: print('t.i. subspaces of dim', d, ':', len(TI[d]))

def clique_partition(C):
    """greedy partition of C into mutually orthogonal sets; try several orders, keep the fewest parts"""
    import random
    best = None
    for seed in range(60):
        rnd = random.Random(seed)
        rem = list(C); rnd.shuffle(rem); rem = set(rem)
        parts = []
        while rem:
            # largest clique in rem (exact via omega-like search returning the set)
            bestset = [[]]
            def r(cur, cand):
                if len(cur) + len(cand) <= len(bestset[0]): return
                if not cand:
                    if len(cur) > len(bestset[0]): bestset[0] = cur
                    return
                T = rnd.choice(sorted(cand))
                r(cur + [T], cand & orth[T]); r(cur, cand - {T})
            r([], set(rem))
            parts.append(bestset[0]); rem -= set(bestset[0])
        if best is None or len(parts) < len(best): best = parts
    return best

rows = Counter(); examples = {}
for d in (1, 2, 3, 4):
    for key, basis in TI[d].items():
        C = C_of(basis)
        Cm = sum(1 << i for i in C)
        inside = sum(1 for m in Omask if m & Cm == m)
        met = sum(1 for m in Omask if m & Cm)
        w = omega(C) if C else 0
        sig = (d, len(C), w, inside, met)
        rows[sig] += 1
        examples.setdefault(sig, basis)
print('dim | |C_g| | omega | octads inside | octads met | cost 9-2d | met/cost | count | clique-partition parts (best of 60 greedy)')
for sig in sorted(rows):
    d, c, w, ins, met = sig
    parts = clique_partition(C_of(examples[sig])) if c else []
    print('%d | %3d | %d | %4d | %4d | %d | %6.1f | %4d | %d parts sizes %s' % (d, c, w, ins, met, 9 - 2 * d,
          met / (9 - 2 * d), rows[sig], len(parts), sorted((len(p) for p in parts), reverse=True)))
print('time %.1f s' % (time.time() - t0))
