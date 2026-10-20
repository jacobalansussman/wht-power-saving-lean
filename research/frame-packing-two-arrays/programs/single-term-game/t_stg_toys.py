#!/usr/bin/env python3
"""t_stg_toys.py: exhaustive single-term game (stg.py) on toy frame games, every paid double coset with m >= 2."""
import sys, time
import numpy as np
from toylab import *
from sweep import GAMES
from stg import *

def paid_classes(G, P):
    Pset = set(P)
    seen = -np.ones(G.n, dtype=np.int32); dcs = []
    for g in range(G.n):
        if seen[g] < 0:
            k = len(dcs); dcs.append(g)
            for p in P:
                for q in P: seen[G.m(G.m(p, g), q)] = k
    out = []; used = set()
    for k, h in enumerate(dcs):
        if h in Pset or k in used: continue
        ki = int(seen[int(G.inv[h])]); used |= {k, ki}
        out.append([h] if ki == k else [h, dcs[ki]])
    return out

def main():
    name = sys.argv[1]; D = int(sys.argv[2]); cap = int(sys.argv[3]) if len(sys.argv) > 3 else 400000
    only = int(sys.argv[4]) if len(sys.argv) > 4 else None
    G, P, _ = GAMES[name]()
    for ci, hr in enumerate(paid_classes(G, P)):
        if only is not None and ci != only: continue
        try:
            g = from_toylab(G, P, hr, name="%s/c%d" % (name, ci))
        except AssertionError as e:
            continue
        if min(g.lv) < 0 or g.m < 2: continue
        print("== %s class %d: frames %d, degree %s, m %d, twists %d, free %d, target orbit %d, level-m frames %d" % (
            name, ci, g.N, sorted({len(x) for x in g.nb}), g.m, g.T, len(g.free), len(g.targets),
            sum(1 for x in g.lv if x == g.m)), flush=True)
        bfs(g, D if D > 0 else 2 * g.m, report=lambda s: print(s, flush=True), max_states=cap)


if __name__ == "__main__":
    main()
