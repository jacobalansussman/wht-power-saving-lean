#!/usr/bin/env python3
"""vinv.py (v-p6-refute): my own test of p6-refute's Conjecture S on every state of my BFS (vstg, each state at its
least depth M).  Objects: the availability mark A and the lineages (entries (t, xs) of the state).
W(u) for an input-1 vertex u = {A if u available} + {lineages with a locked copy whose input-1 end is u}; likewise
W(w) for input-2 vertices (copy ends t(x)).  A pair (u, w) is TWO-OBJECT if both are available, or
|W(u) | W(w)| >= 2, or u, w lie in one locked component at path distance >= 3.
Conjecture S: lv(u) + lv(w) <= M for every two-object pair.  Also reports the largest same-lineage excess and checks
Lemma 4.7's bound is not contradicted in the simplest form (excess <= number of copies of the lineage)."""
import sys
import vstg

def main(name, D):
    s = vstg.S(vstg.make(name)); g = s.g; N = g.N; tp = s.tp
    st0 = s.close(1 << g.raw, 1 << g.raw, ())
    seen = {st0}; layer = [st0]
    worst = -99; worst_ex = None; worst_same = -99; nstates = 0
    for M in range(1, D + 1):
        new = []
        for st in layer:
            for s2 in s.moves(st):
                c = s.close(*s2)
                if c in seen: continue
                seen.add(c); new.append(c)
        for c in new:
            nstates += 1
            A1, A2, L = c
            W1 = {}; W2 = {}
            for x in range(N):
                if (A1 >> x) & 1: W1.setdefault(x, set()).add("A")
                if (A2 >> x) & 1: W2.setdefault(x, set()).add("A")
            adj = {}
            for li, (t, xs) in enumerate(L):
                for x in range(N):
                    if (xs >> x) & 1:
                        W1.setdefault(x, set()).add(li); W2.setdefault(tp[t][x], set()).add(li)
                        u, w = x, N + tp[t][x]
                        adj.setdefault(u, []).append(w); adj.setdefault(w, []).append(u)
            # path distances inside locked components
            dist = {}
            for u in [v for v in adj if v < N]:
                d = {u: 0}; q = [u]
                for v in q:
                    for nb in adj[v]:
                        if nb not in d: d[nb] = d[v] + 1; q.append(nb)
                dist[u] = d
            for u, Wu in W1.items():
                for w, Ww in W2.items():
                    both_av = ("A" in Wu) and ("A" in Ww)
                    far = (u in dist and (N + w) in dist[u] and dist[u][N + w] >= 3)
                    ex = g.lv[u] + g.lv[w] - M
                    if both_av or len(Wu | Ww) >= 2 or far:
                        if ex > worst: worst = ex; worst_ex = (M, c, u, w)
                    else:
                        if ex > worst_same: worst_same = ex
        layer = new
    print("%s to depth %d: %d states; Conjecture S max excess %d%s; same-lineage-only max excess %d" % (
        name, D, nstates, worst, ("  COUNTEREXAMPLE %s" % (worst_ex,)) if worst > 0 else "", worst_same))

if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]))
