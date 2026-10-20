"""p4_where.py WHICH TSING [maxU]: which targeted candidate of probe4.py is the file TSING (the t.sing a running
probe last wrote) for?  Rebuilds probe4.py's targeted candidate list for WHICH exactly as probe4.py does, writes the
first-order system of each candidate in Singular syntax and compares it with the first ideal in TSING.
Prints the index of the matching candidate (= how many candidates were finished before it).  For one remote run (b5:0),
whose copy of probe4.py printed no progress lines (16 Oct)."""
import sys, re
from f2q import Field, Config, popc
from gb2 import supports_in, orders, Ring, system, pstr
from ztest import links, cliques

which, tsing = sys.argv[1], sys.argv[2]; maxU = int(sys.argv[3]) if len(sys.argv) > 3 else 7
F = Field(5); k = int(which.split(":")[1])
line = [l for l in open("survivors_b5w4.txt") if l.startswith("SURVIVES")][k]
ks = [int(v) for v in re.search(r"\(([0-9, ]+)\)", line).group(1).split(",")]
moves = [int(v) for v in re.search(r"dirs \[([0-9, ]+)\]", line).group(1).split(",")]
cfg = Config(F, [F.target(k) for k in ks], moves)
tight, wit = links(cfg); cl = [y for y in cliques(cfg, wit) if len(y) >= 2]
sc = {}
for Q in range(1, len(cfg.Dmask)):
    d = cfg.Dmask[Q]; r = popc(cfg.pmask[Q]); sd = popc(d)
    if sd > maxU: continue
    for y in cl:
        ym = sum(1 << i for i in y)
        if d & ym != ym: continue
        u = d
        for a in y:
            for b in y: u |= wit[(min(a, b), max(a, b))][2]
        if popc(u) > maxU: continue
        key = (sd + r - 2 * len(y), popc(u))
        if u not in sc or key < sc[u]: sc[u] = key
cands = sorted(sc, key=lambda u: (sc[u], u))
txt = open(tsing).read()
ids = [g.replace("\n", "") for g in re.findall(r"I = (.*?);\nG = std", txt, re.S)]
print("candidates", len(cands), "ideals in file", len(ids), flush=True)
# Different U of one size can give the same text (variables are numbered by position in U), so match BOTH ideals
# of the two-order filter (the second order is random.Random(u).choice(ords), as in probe4.py) and list all matches.
import random
hits = []
for ci, u in enumerate(cands):
    U = [i for i in range(cfg.M) if u >> i & 1]
    sup = supports_in(cfg, U)
    if len(sup) * cfg.W ** 2 != ids[0].count(",") + 1: continue   # system() gives W^2 polys per support
    R = Ring(cfg.W, U); ords = list(orders(cfg, U))
    if ",".join(pstr(p, "sing") for p in system(cfg, R, ords[0], sup)) != ids[0]: continue
    first = [ords[0], random.Random(u).choice(ords)]
    ss = [",".join(pstr(p, "sing") for p in system(cfg, R, o, sup)) for o in first]
    if ss[:len(ids)] == ids:
        hits.append(ci); print("MATCH candidate index", ci, "U", U, flush=True)
print("matches", hits, flush=True)
