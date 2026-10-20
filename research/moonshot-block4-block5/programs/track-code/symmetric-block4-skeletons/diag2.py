"""diag2.py -- the structure of a failing skeleton: array 0's closure columns (which generator) and the rows they hit."""
import sys
from symsolve import *
m, W, H, shape, menu, seed, s = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4].split(','), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7])
sk, cost = build_skeleton(m, W, H, shape, menu, seed * 100 + s)
sk.setup_final(6000); prep_fast(sk)
G = sk.G
names = {}
for i in range(W): names[i] = f"x{i}"
for t, mv in enumerate(sk.moves):
    names[mv["rep"]] = f"m{t}"
    for k, im in enumerate(mv["images"]): names[im] = f"m{t}'{k+1}"
tok = []; gd = []
for i, T in enumerate(sk.gens):
    d = {}
    for j, tt in enumerate(T):
        d[tt] = len(tok) + 1; tok.append((i, j))
    gd.append(d)
rows, cols, ent = G.closure(gd, (0, G.KID), 6000)
inv = {r: t for t, r in rows.items()}
print("shape", shape, "cost", cost, "orbits", [1 + len(mv["images"]) for mv in sk.moves], "terms per generator", [len(T) for T in sk.gens])
for (kind, i, b), c in sorted(cols.items(), key=lambda x: x[1]):
    hit = [(r, inv[r][0], G.dist(inv[r][1])) for (r, cc, co) in ent if cc == c]
    print(f" column {c}: generator {names.get(i, i)} -> rows (row, input, dist-to-K): {sorted(set(hit))}")
print(" rows:", [(r, inv[r][0], G.dist(inv[r][1])) for r in range(len(rows))])
# which generators hold a term in the same FRAME as each junk row (but not free-equivalent)?
for r in range(1, len(rows)):
    l, gid = inv[r]
    fk = G.frm(gid)
    same = [(names.get(i, i), j) for i, T in enumerate(sk.gens) for j, (ll, g) in enumerate(T) if ll == l and G.frm(g) == fk]
    print(f" row {r} (input {l}, dist {G.dist(gid)}): generator terms in the same frame: {same[:6]}")
