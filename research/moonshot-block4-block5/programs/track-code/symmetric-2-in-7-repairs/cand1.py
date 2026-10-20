"""cand1.py -- inspect the relaxation-passing 2-in-7 skeleton (shape P:inp,P:prev,L:prev,P:prev, menu 2, seed 66, s 1)"""
import sys
from relax_search import *
m, W, H, shape, menu, seed, s = 4, 2, 'Z2swap', 'P:inp,P:prev,L:prev,P:prev'.split(','), 2, 66, 1
sk, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
print("cost", cost, "orbits", [1 + len(mv["images"]) for mv in sk.moves], "terms per generator", [len(T) for T in sk.gens])
print("summands per move:", [(k, len(sm)) for (k, sm, d) in rec])
print("relaxation residual (60 LM iterations, 4 starts):", relax_res(sk, iters=60, starts=4))
G = sk.G
tok = []; gd = []
for i, T in enumerate(sk.gens):
    d = {}
    for j, tt in enumerate(T):
        d[tt] = len(tok) + 1; tok.append((i, j))
    gd.append(d)
for cap in (6000, 20000, 60000):
    cl = G.closure(gd, (0, G.KID), cap)
    if cl is None:
        print("closure cap", cap, "exceeded"); continue
    rows, cols, ent = cl
    print("closure: rows", len(rows), "cols", len(cols), "entries", len(ent)); break
