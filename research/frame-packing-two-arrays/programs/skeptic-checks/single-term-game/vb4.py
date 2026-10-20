#!/usr/bin/env python3
"""vb4.py (v-p6-refute): my own exact decider for a single-term two-input 2 in D (permissive Lemma G closure), used to
re-check p6-refute's Proposition 4.4 (block-4 cube, no 2 in 7) independently of stg.py.
  - moves and closure: vstg.S (unpruned; checked against the explicit module in vmod.py);
  - canonical dedup (my code) under free right multiplication on each input and the swap, to depth D-2;
  - depth D-1: every successor closed and checked for (I) and the goal;
  - depth D: for every depth-(D-1) state whose two inputs both own a frame (available or locked) at level >= m-1,
    an exact goal-in-one test written from scratch: enumerate the END(S) of every possible new generator (single
    loop or pair edge); the goal holds iff each input still lacking the target has its target vertex in the merged
    component {ends} + their locked components, and that merged component is unlocked (a loop, an available end,
    or both ends already in one locked component = a new cycle).
Usage: vb4.py GAME D [split_layer part nparts]"""
import sys, time
import vstg, os
CHECK = os.environ.get("VB4_CHECK") == "1"

def main():
    name, D = sys.argv[1], int(sys.argv[2])
    split = tuple(int(x) for x in sys.argv[3:6]) if len(sys.argv) > 5 else None
    s = vstg.S(vstg.make(name)); g = s.g; N = g.N; tp = s.tp
    free = [i for i, p in enumerate(tp) if p[g.raw] == g.raw]
    if os.environ.get("VB4_LOW") == "1":
        low = next(x for x in range(N) if g.lv[x] == g.m - 1)
        g.target = {tp[a][low] for a in free}; g.m = g.m - 1; s.tmask = sum(1 << x for x in g.target)
    tmask = s.tmask; m = g.m
    bitsof = lambda M: [x for x in range(N) if (M >> x) & 1]
    def img(p, M):
        r = 0
        for x in range(N):
            if (M >> x) & 1: r |= 1 << p[x]
        return r
    def canon(st):
        A1, A2, L = st; best = None
        for sw in (0, 1):
            if sw: B1, B2, LB = A2, A1, [(s.inv[t], img(tp[t], xs)) for (t, xs) in L]
            else: B1, B2, LB = A1, A2, list(L)
            i1 = [(img(tp[a], B1), a) for a in free]; m1 = min(i1)[0]
            i2 = [(img(tp[b], B2), b) for b in free]; m2 = min(i2)[0]
            if best is not None and (m1, m2) > best[:2]: continue
            for a in [a for v, a in i1 if v == m1]:
                ia = s.inv[a]
                for b in [b for v, b in i2 if v == m2]:
                    key = (m1, m2, tuple(sorted((s.comp[b][s.comp[t][ia]], img(tp[a], xs)) for (t, xs) in LB)))
                    if best is None or key < best: best = key
        return best
    lvmax = lambda A: max(g.lv[x] for x in bitsof(A))
    # created-pair end table: Y[(x0, y0, x2)] = {t(x2) : t(x0) = y0}
    Y = {}
    for (x0, y0), ts in s.cre.items():
        for x2 in s.step[x0]:
            Y[(x0, y0, x2)] = sorted({tp[t][x2] for t in ts})
    tgt = sorted(g.target)
    def goal_in_one(st):
        A1, A2, L = st
        need = [l for l, A in ((0, A1), (1, A2)) if not (A & tmask)]
        if not need: return True
        par = {}
        def f(a):
            while par.get(a, a) != a: a = par[a]
            return a
        adj = {}
        for (t, xs) in L:
            for x in bitsof(xs):
                u, w = x, N + tp[t][x]
                ru, rw = f(u), f(w)
                if ru != rw: par[ru] = rw
                adj.setdefault(u, []).append((w, t)); adj.setdefault(w, []).append((u, t))
        comp = {v: f(v) for v in adj}
        # target vertices still needed: (l, T) for T in the target orbit; goal needs, for each needed l, SOME T
        needv = {l: [(T if l == 0 else N + T) for T in tgt] for l in need}
        def ok(ends, unlocked):
            if not unlocked: return False
            cs = {comp[e] for e in ends if e in comp}
            for l in need:
                if not any((v in ends) or (v in comp and comp[v] in cs) for v in needv[l]): return False
            return True
        av = lambda v: ((A1 >> v) & 1) if v < N else ((A2 >> (v - N)) & 1)
        # singles (loops)
        for x in bitsof(A1):
            for x2 in s.step[x]:
                if ok([x2], True): return True
        for y in bitsof(A2):
            for y2 in s.step[y]:
                if ok([N + y2], True): return True
        def edge(u, w):
            unl = av(u) or av(w) or (u in comp and w in comp and comp[u] == comp[w])
            return ok([u, w], unl)
        for (t, xs) in L:
            for x in bitsof(xs):
                for x2 in s.step[x]:
                    if (xs >> x2) & 1: continue           # the same copy again: no new generator
                    if edge(x2, N + tp[t][x2]): return True
        for x0 in bitsof(A1):
            for y0 in bitsof(A2):
                if (x0, y0) not in s.cre: continue
                for x2 in s.step[x0]:
                    for y2 in Y[(x0, y0, x2)]:
                        if edge(x2, N + y2): return True
        for u0 in [v for v in adj if v < N]:
            seen = {u0: (0, 0)}; stck = [u0]
            while stck:
                v = stck.pop(); tv, d = seen[v]
                for (nb, t) in adj[v]:
                    if nb in seen: continue
                    nt = s.comp[t][tv] if v < N else s.comp[s.inv[t]][tv]
                    seen[nb] = (nt, d + 1); stck.append(nb)
                    if nb >= N and d + 1 >= 3:
                        for x2 in s.step[u0]:
                            if edge(x2, N + tp[nt][x2]): return True
        return False
    print("== %s: frames %d, m %d, twists %d, free %d, target orbit %d, D %d, split %s" % (
        name, N, m, len(tp), len(free), len(tgt), D, split), flush=True)
    st0 = s.close(1 << g.raw, 1 << g.raw, ())
    layer = {canon(st0): st0}; allkeys = set(layer)
    for M in range(1, D - 1):
        t0 = time.time(); new = {}; mx = -1; ngoal = 0
        items = sorted(layer.items())
        if split and M == split[0] + 1:
            items = [it for i, it in enumerate(items) if i % split[2] == split[1]]
            print("   split: part %d/%d of layer %d (%d states)" % (split[1], split[2], split[0], len(items)), flush=True)
        for key, st in items:
            for s2 in s.moves(st):
                c = s.close(*s2); k = canon(c)
                if k in allkeys or k in new: continue
                new[k] = c; mx = max(mx, lvmax(c[0]) + lvmax(c[1]))
                if (c[0] & tmask) and (c[1] & tmask): ngoal += 1
        allkeys |= set(new); layer = new
        print("%s M=%d new canonical states %d max(L1+L2) %d (I) %s goal %d [%.0fs]" % (
            name, M, len(new), mx, "holds" if mx <= M else "VIOLATED", ngoal, time.time() - t0), flush=True)
    t0 = time.time(); n1 = nbad = ng1 = ntest = ng = 0; ex = None; nmis = [0]
    for key, st in layer.items():
        loc = set()
        for s2 in s.moves(st):
            c = s.close(*s2)
            if c in loc or c == st: continue
            loc.add(c); n1 += 1
            a, b = lvmax(c[0]), lvmax(c[1])
            if a + b > D - 1: nbad += 1
            if (c[0] & tmask) and (c[1] & tmask): ng1 += 1
            Q0 = a; Q1 = b
            for (t, xs) in c[2]:
                for x in bitsof(xs):
                    Q0 = max(Q0, g.lv[x]); Q1 = max(Q1, g.lv[tp[t][x]])
            if Q0 < m - 1 or Q1 < m - 1: continue
            ntest += 1
            if CHECK:
                bf = any(((cc[0] & tmask) and (cc[1] & tmask)) for cc in (s.close(*x) for x in s.moves(c)))
                if bf != goal_in_one(c): print("GOAL-IN-ONE MISMATCH", c, bf, flush=True); nmis[0] += 1
            if goal_in_one(c):
                ng += 1
                if ex is None: ex = c
    print("%s M=%d successors %d (I) %s goal %d | M=%d: tested %d, goal-in-one %d [%.0fs]" % (
        name, D - 1, n1, "holds" if not nbad else "VIOLATED %d" % nbad, ng1, D, ntest, ng, time.time() - t0), flush=True)
    if CHECK: print("   brute-force check of goal_in_one on every tested state: %d mismatches" % nmis[0])
    if ex: print("   example depth-%d state with a goal in one: %s" % (D - 1, ex,), flush=True)

if __name__ == "__main__":
    main()
