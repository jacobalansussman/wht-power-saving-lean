"""bothfwd.py -- mw2-absorb4: can the LAST move go forward (to the K frame) on BOTH arrays' parts of the shared move?
For a limit skeleton, take the shared generator's array-0 term h0 and array-1 term h1 (distance 1 each) and sample free
letters b (engine and Wd) with dist(D1 b h0) = 0; count those with dist(D1 b h1) = 0 too.  Usage: bothfwd.py lim.pkl tsh nsamp"""
import sys, pickle
from twist import *
d = pickle.load(open(sys.argv[1], 'rb')); tsh = int(sys.argv[2]); ns = int(sys.argv[3])
sk = build_T(d['m'], d['W'], d['H'], d['recipe'], d['skseed'], d.get('twist', {}))
G = sk.G; n = G.n
gsh = sk.moves[tsh]['rep']
T = sk.gens[gsh]
h0 = [g for (l, g) in T if l == 0]; h1 = [g for (l, g) in T if l == 1]
print("shared generator terms: array0 dists", [G.dist(g) for g in h0], "array1 dists", [G.dist(g) for g in h1])
rng = random.Random(1)
one = both = 0; dists = {}
for it in range(ns):
    if it % 2: b = rand_free_b(G, rng, G.m)
    else:
        sig = list(range(n)); rng.shuffle(sig); b = G.mono(sig, [rng.randrange(4) for _ in range(n)])
    g0 = G.intern(G.mul(G.D1, b, G.MATS[h0[0]]))
    if G.dist(g0) != 0: continue
    one += 1
    g1 = G.intern(G.mul(G.D1, b, G.MATS[h1[0]]))
    dd = G.dist(g1); dists[dd] = dists.get(dd, 0) + 1
    if dd == 0: both += 1
print(f"samples {ns}: forward on array 0: {one}; distance reached on array 1 by those: {dists}; forward on both: {both}")
# also the progress letters the builder would offer
pl = G.progress_letters(h0[0], rng, want=40)
dd = {}
for b in pl:
    g1 = G.intern(G.mul(G.D1, b, G.MATS[h1[0]])); x = G.dist(g1); dd[x] = dd.get(x, 0) + 1
print("builder progress letters for array 0:", len(pl), "-> array-1 distance histogram", dd)
