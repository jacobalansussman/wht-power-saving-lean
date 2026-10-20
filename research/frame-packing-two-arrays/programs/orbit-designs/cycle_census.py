#!/usr/bin/env python3
"""cycle_census.py (pk2-orbits): how often does a CYCLE orbit exist?  For trunk designs (2 in 2m - 1) after one
landing, take every junk row of the delivering states and every landing site (register slot, later output, source
walker); a cycle sigma exists at the site iff a SECOND junk row of the same input sits on a frame X2 != X that the same
source walker can also reach at the same output (same distance), since then sigma = chart(Y2 -> Y) * stab works.
Usage: cycle_census.py m trials seed"""
import sys, os, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "heights")); sys.path.insert(0, os.path.join(HERE, "..", "frame-game")); sys.path.insert(0, HERE)
from fpa import *
import design as D, orb
m, trials, seed = map(int, sys.argv[1:4])
S = Sp(m); rng = random.Random(seed)
tot_sites = cyc_sites = designs = des_with = 0
for tr in range(trials):
    P = D.Prog(S, m, rng, dup=0)
    if not P.ok: continue
    res = P.evaluate()
    orph = P.orphan_list(res)
    if orph: P.try_repair(rng.choice(orph))
    st = P.run()
    dl = orb.deliver_idx(P, 2)
    rows = orb.junk_rows(S, st, dl, 2)
    O = orb.Orbit(S, P, 2, rng, ["stab"])
    designs += 1; got = False
    for (l, j) in rows:
        X = S.pos(j)
        for (reg, k, t, src, g, w) in O.sites(l, X):
            tot_sites += 1
            lw = S.level(w)
            ok = any(r[0] == l and S.pos(r[1]) != X and
                     S.rank_rows(list(S.pos(g)) + list(S.pos(r[1]))) - S.m == lw for r in rows)
            if ok: cyc_sites += 1; got = True
    des_with += got
print(f"m={m}: {designs} designs after one landing; landing sites on delivering-state junk: {tot_sites}; "
      f"sites admitting a cycle orbit: {cyc_sites}; designs with at least one: {des_with}")
