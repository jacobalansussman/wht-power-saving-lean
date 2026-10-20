"""aut2.py (t-km): random automorphisms of P112 in GL(9,2) (label set and label dot parity preserved), by randomised
backtracking (first solution per seed), then random products to find elements of order 7 and 3 and their
port permutations.  Saves gens.npz (list of L as images of e_0..e_8, and the port permutations).
usage: python -I aut2.py NSEEDS"""
import sys, time, random
sys.path.insert(0, __import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '.'))
from model import *
import numpy as np
from collections import Counter
LS = set(LAB)
def _spanset(bs):
    s = {0}
    for b in bs: s |= {x ^ b for x in s}
    return s
basis = [LAB[0]]
for _k in range(8):
    _best = None
    for u in LAB:
        if u in _spanset(basis): continue
        n = len(_spanset(basis + [u]) & LS)
        if _best is None or n > _best[0]: _best = (n, u)
    basis.append(_best[1])
span = {0: 0}
for i, b in enumerate(basis):
    span.update({x ^ b: c | (1 << i) for x, c in list(span.items())})
coord = {u: span[u] for u in LAB}
dot = lambda a, b: par(a & b)
G0 = [[dot(a, b) for b in basis] for a in basis]
bylev = [[] for _ in range(10)]
for u in LAB: bylev[coord[u].bit_length()].append(coord[u])
def img(c, ims):
    x = 0; i = 0
    while c:
        if c & 1: x ^= ims[i]
        c >>= 1; i += 1
    return x
def first(rng):
    order = list(LAB)
    def bt(ims):
        k = len(ims)
        if k == 9: return ims
        cands = order[:]; rng.shuffle(cands)
        for cand in cands:
            if any(dot(cand, ims[j]) != G0[k][j] for j in range(k)): continue
            ims2 = ims + [cand]
            if all(img(c, ims2) in LS for c in bylev[k + 1]):
                r = bt(ims2)
                if r: return r
        return None
    return bt([])
def to_std(ims): return [img(span[1 << i], ims) for i in range(9)]
def apply(L, x): return img(x, L)
PI = {u: i for i, u in enumerate(LAB)}
def perm_of(L): return np.array([PI[apply(L, u)] for u in LAB])
def compose(A, B): return [apply(A, apply(B, 1 << i)) for i in range(9)]   # A o B
def order(L):
    k = 1; M = L
    ident = [1 << i for i in range(9)]
    while M != ident: M = compose(L, M); k += 1
    return k
def power(L, n):
    M = [1 << i for i in range(9)]
    for _ in range(n): M = compose(L, M)
    return M
rng = random.Random(1)
gens = []
t0 = time.time()
for s in range(int(sys.argv[1])):
    ims = first(rng); L = to_std(ims)
    assert all(apply(L, u) in LS for u in LAB)
    gens.append(L)
print('random automorphisms', len(gens), 'orders', [order(L) for L in gens], '%.1fs' % (time.time() - t0), flush=True)
found = {}
for t in range(4000):
    a = rng.choice(gens); b = rng.choice(gens); c = rng.choice(gens)
    M = compose(a, compose(b, c))
    o = order(M)
    for p in (7, 3, 9, 5):
        if o % p == 0 and p not in found:
            found[p] = power(M, o // p)
    if 7 in found and 3 in found: break
print('found prime orders:', sorted(found), flush=True)
np.savez('gens.npz', gens=np.array(gens), **{'g%d' % p: np.array(L) for p, L in found.items()},
         **{'perm%d' % p: perm_of(L) for p, L in found.items()})
for p, L in found.items():
    pm = perm_of(L)
    cyc = Counter()
    seen = set()
    for i in range(len(pm)):
        if i in seen: continue
        j = i; n = 0
        while j not in seen: seen.add(j); j = pm[j]; n += 1
        cyc[n] += 1
    print('element of order', p, 'port cycle type', dict(cyc))
