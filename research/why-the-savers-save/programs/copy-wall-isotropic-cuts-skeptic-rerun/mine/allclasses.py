"""check the clique-partition certificate on EVERY isotropic h (not one example per signature)"""
import sys, random
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-isotropic-cuts'))
from base import *
from collections import Counter
iso = [h for h in range(1, 512) if not dot9(h, h)]
print('isotropic nonzero h:', len(iso))
def maxclique(rem, rnd):
    best = [[]]
    def r(cur, cand):
        if len(cur) + len(cand) <= len(best[0]): return
        if not cand:
            if len(cur) > len(best[0]): best[0] = cur
            return
        T = rnd.choice(sorted(cand))
        r(cur + [T], cand & orth[T]); r(cur, cand - {T})
    r([], set(rem)); return best[0]
def partition(C, tries=30):
    bestp = None
    for seed in range(tries):
        rnd = random.Random(seed); rem = set(C); parts = []
        while rem:
            q = maxclique(rem, rnd); parts.append(q); rem -= set(q)
        if bestp is None or len(parts) < len(bestp): bestp = parts
    return bestp
worst = 0; hist = Counter()
for h in iso:
    C = C_of([h])
    parts = partition(C)
    # verify certificate
    assert sorted(x for p in parts for x in p) == sorted(C)
    for p in parts:
        for a in p:
            for b in p:
                if a != b: assert b in orth[a]
    hist[(len(C), omega(C), len(parts))] += 1
    worst = max(worst, len(parts))
for k in sorted(hist): print('|C_h| %d omega %d parts %d : %d classes' % (k + (hist[k],)))
print('max parts over all 255 isotropic h:', worst)
# every port lies in some C_h? and h with C_h empty?
print('min |C_h|', min(len(C_of([h])) for h in iso))
