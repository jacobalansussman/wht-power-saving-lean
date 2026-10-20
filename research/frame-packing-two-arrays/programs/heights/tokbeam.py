#!/usr/bin/env python3
"""tokbeam.py (p6-prove, 15 Oct 2026): beam search in the permissive single-term two-input module game
(same rules as tokgame2.py) at blocks 3 and 4, looking for a state with max level(S1) + max level(S2) > M.
Twists: random Sp(2m,2) elements with sigma(x0) = y0, built as c_y0 o s o c_x0^-1 (s random in Stab(L0)), plus a
few structured ones (the swap Kbar-type element mapping L0 <-> LT).  A hit is only a CANDIDATE: the permissive
closure unlocks every cycle, which the exact module may not.  Usage: tokbeam.py m beam width seed
"""
import sys, random
from collections import deque
from frames import Frames

def main(m, beam, width, seed):
    rng = random.Random(seed)
    F = Frames(m); N = F.N; lv = F.dist[F.L0]
    # c[x]: an Sp element with c[x](L0) = x (BFS over transvections)
    gens = [F.tv(v) for v in F.vecs]
    c = {F.L0: F.ident()}; q = deque([F.L0])
    while q:
        x = q.popleft()
        for g in gens:
            h = F.mul(g, c[x]); y = F.apply_pt(h, F.L0)
            if y not in c: c[y] = h; q.append(y)
    pcache = {}
    def perm(g):
        if g not in pcache: pcache[g] = F.perm(g)
        return pcache[g]
    def invp(p):
        r = [0] * N
        for i, j in enumerate(p): r[j] = i
        return tuple(r)
    def compp(a, b): return tuple(a[b[i]] for i in range(N))
    stabgens = [g for g in gens if F.apply_pt(g, F.L0) == F.L0]
    def rand_stab():
        g = F.ident()
        for _ in range(12): g = F.mul(rng.choice(stabgens), g)
        return g
    cinv = {x: invp(perm(c[x])) for x in c}
    def rand_twist(x0, y0):
        return compp(perm(c[y0]), compp(perm(rand_stab()), cinv[x0]))
    # a swap-type element: LT = c[LT](L0); kswap maps L0 -> LT and LT -> ? (we just use c[LT] composed with stab)
    def closure(S1, S2, lins):
        while True:
            par = {}
            def f(a):
                while par.setdefault(a, a) != a: a = par[a]
                return a
            E = []
            for (s, xs) in lins:
                for x in xs:
                    u, w = (1, x), (2, s[x]); E.append((u, w)); par[f(u)] = f(w)
            comps = {}
            for v in list(par): comps.setdefault(f(v), [set(), 0])[0].add(v)
            for (u, w) in E: comps[f(u)][1] += 1
            changed = False
            for V, ne in comps.values():
                if ne >= len(V) or any((x in S1) if sd == 1 else (x in S2) for (sd, x) in V):
                    for (sd, x) in V:
                        if sd == 1 and x not in S1: S1 = S1 | {x}; changed = True
                        if sd == 2 and x not in S2: S2 = S2 | {x}; changed = True
            if not changed: return S1, S2
    def score(st):
        S1, S2, lins, M = st
        L1 = max(lv[x] for x in S1); L2 = max(lv[x] for x in S2)
        pot = max([lv[x] + lv[s[x]] for (s, xs) in lins for x in xs] + [0])
        return (L1 + L2 - M) + 0.25 * pot - 0.5 * (L1 + L2 - M < 0) + rng.random() * 0.01
    start = (frozenset([F.L0]), frozenset([F.L0]), (), 0)
    states = [start]; best = -99
    for depth in range(1, 2 * m + 1):
        cand = []
        for (S1, S2, lins, M) in states:
            opts = []
            for x in S1:
                for x2 in F.nb[x]:
                    if x2 not in S1: opts.append(("s1", x2))
            for x in S2:
                for x2 in F.nb[x]:
                    if x2 not in S2: opts.append(("s2", x2))
            for li, (s, xs) in enumerate(lins):
                for x in xs:
                    for x2 in F.nb[x]:
                        if x2 not in xs: opts.append(("ln", li, x2))
            rng.shuffle(opts); opts = opts[:width]
            for _ in range(width // 2):
                x0 = rng.choice(sorted(S1)); y0 = rng.choice(sorted(S2))
                opts.append(("new", rand_twist(x0, y0), rng.choice(F.nb[x0])))
            for o in opts:
                if o[0] == "s1": n = (S1 | {o[1]}, S2, lins)
                elif o[0] == "s2": n = (S1, S2 | {o[1]}, lins)
                elif o[0] == "ln":
                    s, xs = lins[o[1]]; nl = list(lins); nl[o[1]] = (s, xs + (o[2],)); n = (S1, S2, tuple(nl))
                else: n = (S1, S2, lins + ((o[1], (o[2],)),))
                a, b = closure(n[0], n[1], n[2])
                cand.append((a, b, n[2], M + 1))
        cand.sort(key=score, reverse=True)
        seen = set(); states = []
        for st in cand:
            k = (st[0], st[1], tuple(sorted((s[:3], xs) for (s, xs) in st[2])))
            if k in seen: continue
            seen.add(k); states.append(st)
            if len(states) >= beam: break
        top = max(max(lv[x] for x in s[0]) + max(lv[x] for x in s[1]) - s[3] for s in states)
        best = max(best, top)
        print("m=%d depth=%d beam=%d best (L1+L2-M) in beam=%d" % (m, depth, len(states), top), flush=True)
    print("RESULT m=%d seed=%d: max over search of L1+L2-M = %d (%s)" % (m, seed, best, "violation candidate" if best > 0 else "no violation"))

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
