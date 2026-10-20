#!/usr/bin/env python3
"""vsym.py (v-p6-refute): test stg.canon (Lemma 4.1) on states of my BFS: for random symmetries sigma (free right
multiplication on each input independently, and the input swap), canon(sigma s) == canon(s), sigma commutes with the
move set ({canon(succ)} equal), and sigma preserves L1 + L2 and the goal."""
import sys, random
import os as _os; sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "..", "single-term-game"))
import stg, vstg
from vcross import build
def apply(s, G, st, a, b, sw):
    A1, A2, L = st; N = G.N
    if sw:
        A1, A2 = A2, A1; L = tuple(sorted((s.inv[t], sum(1 << s.tp[t][x] for x in range(N) if (xs >> x) & 1)) for (t, xs) in L))
    pa, pb = s.tp[a], s.tp[b]
    B1 = sum(1 << pa[x] for x in range(N) if (A1 >> x) & 1); B2 = sum(1 << pb[x] for x in range(N) if (A2 >> x) & 1)
    LL = tuple(sorted((s.comp[b][s.comp[t][s.inv[a]]], sum(1 << pa[x] for x in range(N) if (xs >> x) & 1)) for (t, xs) in L))
    return (B1, B2, LL)
name, D, ntest = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
s, G = build(name); g = s.g; rng = random.Random(3)
free = [i for i, p in enumerate(s.tp) if p[g.raw] == g.raw]
st0 = s.close(1 << g.raw, 1 << g.raw, ()); seen = {st0}; layer = [st0]; allst = [st0]
for M in range(1, D + 1):
    new = []
    for st in layer:
        for s2 in s.moves(st):
            c = s.close(*s2)
            if c not in seen: seen.add(c); new.append(c)
    layer = new; allst += new
bad = 0; tm = G.tmask
lv = lambda A: max(g.lv[x] for x in range(g.N) if (A >> x) & 1)
for k in range(ntest):
    st = rng.choice(allst); a, b, sw = rng.choice(free), rng.choice(free), rng.randrange(2)
    st2 = apply(s, G, st, a, b, sw)
    if s.close(*st2) != st2: bad += 1; print("image not closed"); continue
    if stg.canon(G, *st) != stg.canon(G, *st2): bad += 1; print("canon differs", st, st2)
    if lv(st[0]) + lv(st[1]) != lv(st2[0]) + lv(st2[1]): bad += 1
    if bool((st[0] & tm) and (st[1] & tm)) != bool((st2[0] & tm) and (st2[1] & tm)): bad += 1
    S1 = {stg.canon(G, *s.close(*x)) for x in s.moves(st)}; S2 = {stg.canon(G, *s.close(*x)) for x in s.moves(st2)}
    if S1 != S2: bad += 1; print("successor classes differ")
print(name, "depth", D, "states", len(allst), "symmetry tests", ntest, "failures", bad)
