#!/usr/bin/env python3
"""t_stg_2m1.py: for every paid class with m = 3 of the named toy games, decide EXACTLY whether a single-term
two-input 2 in 2m - 1 = 2 in 5 exists (bfs3: dedup to 3 moves, every 4-move state checked for (I), exact
one-move goal test at 5)."""
import sys, time
from stg import *
from sweep import GAMES
from t_stg_toys import paid_classes
for name in sys.argv[1:]:
    G, P, _ = GAMES[name]()
    for ci, hr in enumerate(paid_classes(G, P)):
        try: g = from_toylab(G, P, hr, name="%s/c%d" % (name, ci))
        except AssertionError: continue
        if min(g.lv) < 0 or g.m != 3: continue
        print("== %s class %d: frames %d, degree %s, m %d, twists %d, free %d" % (name, ci, g.N, sorted({len(x) for x in g.nb}), g.m, g.T, len(g.free)), flush=True)
        bfs3(g, 2 * g.m - 1, report=lambda s: print(s, flush=True))
