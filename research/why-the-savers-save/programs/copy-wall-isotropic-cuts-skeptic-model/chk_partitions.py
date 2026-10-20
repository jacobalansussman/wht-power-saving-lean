"""v-g-break-1-model: check, for EVERY isotropic line class h (not one representative per signature),
a partition of C_h into orthogonal cliques with few parts; and recount the K0=7 exact kills."""
import sys
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-isotropic-cuts'))
from base import *
import random
from collections import Counter
iso = [h for h in range(1, 512) if not dot9(h, h)]
def maxclique(rem, rnd):
    best = [[]]
    def r(cur, cand):
        if len(cur) + len(cand) <= len(best[0]): return
        if not cand:
            best[0] = cur; return
        T = rnd.choice(sorted(cand))
        r(cur + [T], cand & orth[T]); r(cur, cand - {T})
    r([], set(rem)); return best[0]
worst = 0; dist = Counter(); kills = Counter()
for h in iso:
    C = C_of([h])
    bestp = None
    for seed in range(15):
        rnd = random.Random(seed); rem = set(C); parts = []
        while rem:
            q = maxclique(rem, rnd); parts.append(q); rem -= set(q)
        if bestp is None or len(parts) < bestp: bestp = len(parts)
    # verify cliques are orthogonal (implicitly by construction) and count
    dist[(len(C), bestp)] += 1; worst = max(worst, bestp)
    comp = set(range(N)) - set(C)
    kills[(len(C), omega(comp))] += 1
print('classes', len(iso)); print('(|C_h|, parts) distribution:', sorted(dist.items())); print('worst parts', worst)
print('(|C_h|, omega(complement)):', sorted(kills.items()))
print('exact K0=7 kills (omega(comp) >= 8):', sum(c for (s, w), c in kills.items() if w >= 8))
