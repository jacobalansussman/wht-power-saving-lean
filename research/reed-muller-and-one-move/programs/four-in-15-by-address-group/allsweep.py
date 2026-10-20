#!/usr/bin/env python3
"""allsweep.py (w8-rmexact): for EVERY cyclic subgroup H of GL(4,2) of the given orders, 4 in 15 with the moves spread
evenly (all 15 directions once: every direction gets n = |H| virtual moves), run the three exact tests in turn:
deletion, E-pattern matching (epat.py), rectangle partitions on 4-point supports (rcover.py).  Prints a table of
verdicts by (order, orthogonal, verdict)."""
import sys, itertools
sys.path.insert(0, '.')
from forms import setup, pair
orders = [int(v) for v in sys.argv[1].split(',')]
part = sys.argv[2] if len(sys.argv) > 2 else '0/1'; pi_, pn_ = map(int, part.split('/'))
W, M = 4, 15
m = 4; D = setup(m); GL = D['GL']; supp = D['supp']; nf = D['nf']; dirs = D['dirs']; mm, tr = D['matmul'], D['transpose']
I = [1 << i for i in range(m)]; Ienc = D['enc_sym'](I)
S = [frozenset(supp[F]) for F in range(nf)]
subs = {}
for g in GL:
    grp = [I]; h = g
    while h != I: grp.append(h); h = mm(g, h)
    if len(grp) in orders: subs[frozenset(tuple(x) for x in grp)] = grp
subs = [subs[k] for k in sorted(subs, key=lambda k: sorted(k))][pi_::pn_]
def aug(r, adj, match, seen):
    for c in adj[r]:
        if c in seen: continue
        seen.add(c)
        if c not in match or aug(match[c], adj, match, seen): match[c] = r; return True
    return False
def verdict(gr):
    n = len(gr); T = [D['enc_sym'](mm(h, tr(h))) for h in gr]
    nodd = [sum(pair(F, s) for s in T) for F in range(nf)]
    C = {x: n for x in dirs}
    sl = [n * len(S[F]) - W * nodd[F] for F in range(nf)]
    if min(sl[1:]) < 0: return 'deletion', 0
    tight = [F for F in range(1, nf) if sl[F] == 0 and nodd[F] > 0]
    zero = set()
    for F in tight:
        for x in S[F]:
            for y in S[F]: zero.add((x, y))
    if tight:
        for F in range(1, nf):
            P = sorted(S[F]); need = sl[F]
            if need == 0: continue
            rows = [(x, a) for x in P for a in range(n)]
            adj = {r: [(y, b) for y in P for b in range(n) if (r[0], y) not in zero] for r in rows}
            match = {}
            if sum(1 for r in rows if aug(r, adj, match, set())) < need: return 'E-pattern', len(tight)
    msl = {}
    def ms(U):
        if U not in msl: msl[U] = min(sl[G] for G in range(1, nf) if U <= S[G])
        return msl[U]
    for F in range(1, nf):
        P = sorted(S[F])
        if len(P) > 4 or sl[F] == 0: continue
        cells = [(x, y) for x in P for y in P if (x, y) not in zero]
        if not cells: return 'rect-cover', len(tight)
        ci = {c: i for i, c in enumerate(cells)}; full = (1 << len(cells)) - 1
        best = {}
        subsets = [frozenset(c) for r in range(1, len(P) + 1) for c in itertools.combinations(P, r)]
        for R in subsets:
            for Cc in subsets:
                mk = 0
                for x in R:
                    for y in Cc:
                        if (x, y) in ci: mk |= 1 << ci[(x, y)]
                if not mk: continue
                u = min(n * len(R), n * len(Cc), ms(R | Cc))
                if mk not in best or u < best[mk]: best[mk] = u
        rects = list(best.items()); INF = 10 ** 9; dp = [INF] * (full + 1); dp[0] = 0
        for st in range(full):
            if dp[st] == INF: continue
            low = (~st) & (st + 1)
            for mk, u in rects:
                if mk & low and not (mk & st):
                    ns = st | mk
                    if dp[st] + u < dp[ns]: dp[ns] = dp[st] + u
        if dp[full] < sl[F]: return 'rect-cover', len(tight)
    return 'pass', len(tight)
from collections import Counter
tab = Counter()
for gr in subs:
    T = [D['enc_sym'](mm(h, tr(h))) for h in gr]
    orth = all(t == Ienc for t in T)
    v, nt = verdict(gr)
    tab[(len(gr), orth, v, nt > 0)] += 1
    print("%s order %d gen %s tight forms %d" % (v.upper(), len(gr), gr[1], nt), flush=True)
for k in sorted(tab, key=str): print("order %2d orthogonal %-5s verdict %-10s has tight forms %-5s : %d" % (k[0], k[1], k[2], k[3], tab[k]))
