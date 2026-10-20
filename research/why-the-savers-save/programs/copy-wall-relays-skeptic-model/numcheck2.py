import sys, itertools
sys.dont_write_bytecode = True
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../copy-wall-relays'))
from collections import Counter
import hub_e8_skeleton as W
from patterns import hamming_e8_hub
Lg, Gg = hamming_e8_hub()
GW = W.GI
P = [39, 52, 57, 71, 76, 78, 103, 109, 114]
def prof(G, P):
    live = [i for i in range(120) if i not in P]
    c = Counter(sum(1 for p in P if G[i][p] != 0) for i in live)
    d = Counter(sum(1 for q in P if q != p and G[p][q] != 0) for p in P)
    return sorted(c.items()), sorted(d.items())
print('word P in word (#370 e8.py) numbering :', prof(GW, P))
print('same indices in hamming_e8 numbering  :', prof(Gg, P))
# map word -> geo
gi = {}
for i, r in enumerate(Lg):
    gi[tuple(r)] = i; gi[tuple(-x for x in r)] = i
w2g = [gi.get(tuple(int(x) for x in r)) for r in W.LINES]
print('unmapped', sum(x is None for x in w2g))
Pg = sorted(w2g[p] for p in P)
print('word P mapped to geo numbering:', Pg, prof(Gg, Pg))
