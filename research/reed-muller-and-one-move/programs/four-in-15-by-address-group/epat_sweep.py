#!/usr/bin/env python3
"""epat_sweep.py (w8-rmexact): deletion + E-pattern test (epat.py) for 4 in 15 with the moves spread evenly
(all 15 directions once, so every direction gets exactly n = |H| virtual moves, for any H), over EVERY cyclic
subgroup of GL(4,2) of the given orders.  Exact integers."""
import sys, itertools
sys.path.insert(0, '.')
from forms import setup, pair
orders = [int(v) for v in sys.argv[1].split(',')]
W, M = 4, 15
m = 4; D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']
I = [1 << i for i in range(m)]; Ienc = D['enc_sym'](I)
subs = {}
for g in GL:
    grp = [I]; h = g
    while h != I: grp.append(h); h = mm(g, h)
    if len(grp) in orders: subs[frozenset(tuple(x) for x in grp)] = grp
def aug(r, adj, match, seen):
    for c in adj[r]:
        if c in seen: continue
        seen.add(c)
        if c not in match or aug(match[c], adj, match, seen): match[c] = r; return True
    return False
from collections import Counter
res = Counter()
for key, gr in sorted(subs.items(), key=lambda kv: sorted(kv[0])):
    n = len(gr); T = [D['enc_sym'](mm(h, tr(h))) for h in gr]
    orth = all(t == Ienc for t in T)
    nodd = [sum(pair(F, S) for S in T) for F in range(nf)]
    C = n
    sl = [C * len(supp[F]) - W * nodd[F] for F in range(nf)]
    if min(sl[1:]) < 0: res[(n, orth, 'deletion')] += 1; continue
    tight = [F for F in range(1, nf) if sl[F] == 0 and nodd[F] > 0]
    zero = set()
    for F in tight:
        for x in supp[F]:
            for y in supp[F]: zero.add((x, y))
    verdict = 'pass'
    if tight:
        for F in range(1, nf):
            P = sorted(supp[F]); need = sl[F]
            if need == 0: continue
            # direction-level b-matching with capacity C per direction = matching on C copies
            rows = [(x, a) for x in P for a in range(C)]
            adj = {r: [(y, b) for y in P for b in range(C) if (r[0], y) not in zero] for r in rows}
            match = {}
            mmv = sum(1 for r in rows if aug(r, adj, match, set()))
            if mmv < need: verdict = 'E-pattern'; break
    res[(n, orth, verdict, len(tight) > 0)] += 1
for k in sorted(res, key=str): print(k, res[k])
