#!/usr/bin/env python3
"""t_stg_sub.py: (I) and the goal in engine SUB-GAMES (cube B_m = apartment with Weyl twists; Hamming H(k,3) =
products of block-1 frame games), and toy games, with bfs2 (last layer without dedup)."""
import sys
from stg import *
from sweep import GAMES
from t_stg_toys import paid_classes
spec, D = sys.argv[1], int(sys.argv[2]); cap = int(sys.argv[3]) if len(sys.argv) > 3 else 2000000
if spec.startswith("H"):
    k = int(spec[1]); g = hamming_game(k, wreath=spec.endswith("w"))
else:
    nm, ci = spec.split(":")
    if nm == "B5":
        from sweep import hyperoct; G, P, _ = hyperoct(5)
    else: G, P, _ = GAMES[nm]()
    g = from_toylab(G, P, paid_classes(G, P)[int(ci)], name=spec)
print("== %s: frames %d, degrees %s, m %d, twists %d, free %d, target orbit %d" % (
    g.name, g.N, sorted({len(x) for x in g.nb}), g.m, g.T, len(g.free), len(g.targets)), flush=True)
sp = tuple(int(x) for x in sys.argv[4].split('/')) if len(sys.argv) > 4 else None
if len(sys.argv) > 5 and sys.argv[5] == 'two':
    bfs3(g, D, report=lambda s: print(s, flush=True), split=sp)
else:
    bfs2(g, D, report=lambda s: print(s, flush=True), max_states=cap, split=sp)
