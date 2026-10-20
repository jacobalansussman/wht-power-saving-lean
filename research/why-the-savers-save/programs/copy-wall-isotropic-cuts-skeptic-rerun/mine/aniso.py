"""1-dim cuts <h> for ANISOTROPIC h too (h.h = 1): |C_h|, omega, clique-partition parts.  Theorem C gives cost 7 to a
total read at an 8-dim span w-perp containing h and served inside C_h; M1 omits these.  Does the K >= 9 counting
(K_h >= omega(T_h) >= |T_h| / parts) still kill K = 9 if they are allowed?"""
import sys, random
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../copy-wall-isotropic-cuts'))
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from base import *
from collections import Counter
import importlib.util
spec = importlib.util.spec_from_file_location('ac', __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), 'allclasses.py'))
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
def partition(C, tries=20):
    bestp = None
    for seed in range(tries):
        rnd = random.Random(seed); rem = set(C); parts = []
        while rem:
            q = maxclique(rem, rnd); parts.append(q); rem -= set(q)
        if bestp is None or len(parts) < len(bestp): bestp = parts
    return bestp
an = [h for h in range(1, 512) if dot9(h, h)]
labset = set(UI)
hist = Counter(); worst = 0
for h in an:
    C = C_of([h]); parts = partition(C) if C else []
    kind = 'port label' if h in lab else ('removed-octad label' if h in UI[:8] else 'other')
    w = omega(C) if C else 0
    hist[(kind, len(C), w, len(parts))] += 1
    worst = max(worst, len(parts))
for k in sorted(hist): print('%-20s |C_h| %3d omega %d parts %2d : %d' % (k + (hist[k],)))
print('anisotropic h:', len(an), ' max parts:', worst, ' => pure 1-dim-cut designs need K >= ceil(112/%d) = %d' % (worst, -(-112 // worst)))
