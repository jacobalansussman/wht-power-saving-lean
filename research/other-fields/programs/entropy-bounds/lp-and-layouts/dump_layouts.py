# dump_layouts.py -- write the frames (Lagrangian subspaces of F_2^6, as bases of (a|b) bit vectors, a = low 3 bits)
# of the block-3 layouts where the cut LP is 13 or 14, with the distance and adjacency facts that define them.
import json, dpg
verts, idx, adj, dist = dpg.build(3)
lays = {'cut13': [(0, 95), (0, 80), (80, 95), (12, 80), (84, 60)], 'cut14': [(0, 128), (0, 124), (124, 128), (8, 124), (76, 44)]}
def vec(v): return f"({v & 7:03b}|{v >> 3:03b})"
out = {}
for name, D in lays.items():
    fr = sorted({x for p in D for x in p})
    out[name] = dict(sessions=D, frames={x: [vec(b) for b in dpg.basis_of(verts[x])] for x in fr},
                     distances={f"{x}-{y}": dist[x][y] for x in fr for y in fr if x < y})
    print(name, "sessions", D)
    for x in fr: print("   frame", x, "basis", out[name]['frames'][x])
    print("   distances", out[name]['distances'])
json.dump(out, open('layouts_block3.json', 'w'), indent=1)
