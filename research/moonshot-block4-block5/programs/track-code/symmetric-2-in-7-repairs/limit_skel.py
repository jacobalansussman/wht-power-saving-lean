"""limit_skel.py -- mw2-absorb4: the projective limit of mw-sym4's border 2-in-7.
In the seed-66/1 skeleton the relaxation residual is |Ft-1|^2 + |t|^2 (1+|a1|^2)(1+|a3|^2), F = (1+a1)(1+a3)(1+a5):
it tends to 0 as a5 -> infinity, i.e. as the FIRST summand of the shared (L) move gets coefficient 0 relative to the
second.  This script builds that limit skeleton (L move with the second summand only), checks the relaxation, and
measures the full closure.  Usage: limit_skel.py seed start [keep_index] [cap]"""
import sys, pickle, time
from relax_search import *
from relax_sym import relax_matrix

m, W, H, shape, menu = 4, 2, 'Z2swap', 'P:inp,P:prev,L:prev,P:prev'.split(','), 2
seed, s = int(sys.argv[1]), int(sys.argv[2])
keep = int(sys.argv[3]) if len(sys.argv) > 3 else 1
cap = int(sys.argv[4]) if len(sys.argv) > 4 else 60000
sk0, cost, rec = record_recipe(m, W, H, shape, menu, seed * 100 + s)
rec2 = [(k, list(sm), d) for (k, sm, d) in rec]
k2, sm2, d2 = rec2[2]
rec2[2] = (k2, [sm2[keep]], d2)
sk = build_from_recipe(m, W, H, rec2, seed * 100 + s)
G = sk.G
print("terms per generator", [len(T) for T in sk.gens], "orbits", [1 + len(mv['images']) for mv in sk.moves], "fixed", sk.fixed, "nalpha", sk.nalpha)
A, rowinfo, Krow, free_idx = relax_matrix(sk)
names = [f"a{a}" for a in free_idx]
for r, (dd, l, g) in enumerate(rowinfo):
    ents = [(i, A[(r, i)]) for i in range(len(sk.gens)) if (r, i) in A]
    print(f"row {r:2d} dist {dd} arr {l}: " + " | ".join(f"g{i}: {p.__repr__(names)}" for i, p in ents))
pickle.dump({"recipe": rec2, "skseed": seed * 100 + s}, open(f"limit_{seed}_{s}_k{keep}.pkl", "wb"))
t0 = time.time()
tok = []; gd = []
for i, T in enumerate(sk.gens):
    d = {}
    for j, tt in enumerate(T):
        d[tt] = len(tok) + 1; tok.append((i, j))
    gd.append(d)
for l in range(W):
    cl = G.closure(gd, (l, G.KID), cap)
    if cl is None:
        print(f"array {l}: closure cap {cap} exceeded ({time.time()-t0:.0f}s)"); continue
    rows, cols, ent = cl
    print(f"array {l}: closure rows {len(rows)} cols {len(cols)} entries {len(ent)} ({time.time()-t0:.0f}s)", flush=True)
