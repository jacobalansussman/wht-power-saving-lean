#!/usr/bin/env python3
"""tokgame2.py (p6-prove, 15 Oct 2026): the single-term two-input MODULE game with the full (permissive) closure,
exhaustive at block 2 (Lemma M, Lemmas U and G of PROOF.md section 4).

Quotient picture (positions = Lagrangians of Gamma_2, twists in Sp(4,2) as permutations of the 15 points):
  singles: S1, S2 = positions where input 1 / input 2 is available alone;
  lineages: (sigma, set of side-1 positions x of its moved copies); a copy is the edge (1,x)-(2,sigma x).
Closure (permissive, safe for lower bounds): build the graph on vertices (1,x), (2,y) with one edge per (lineage,
copy) (copies of one lineage never form a cycle by themselves); a component is UNLOCKED if it contains an available
vertex or a cycle (|edges| >= |vertices|); every vertex of an unlocked component becomes available.  Repeat.
Moves (cost 1): move a single; move a copy of a lineage (same lineage); move a created pair (x0 in S1, sigma x0 in S2,
new lineage); move a COMPOSITE: in a locked (tree) component, the unique path from (1,u) to (2,w) gives twist
sigma_k o sigma_{k-1}^-1 o ... o sigma_1 (alternating), a new lineage unless the path is one edge.
Check at every state: (I) max level(S1) + max level(S2) <= M.  Goal: LT in S1 and in S2.
"""
import sys
from frames import Frames

def main(m, D):
    F = Frames(m); N = F.N
    G = F.all_elements(); gi = {g: k for k, g in enumerate(G)}
    inv = [gi[tuple(sorted(range(N), key=lambda i: g[i]))] for g in G]
    comp = lambda a, b: gi[tuple(G[a][G[b][i]] for i in range(N))]          # a after b
    stab = [k for k, g in enumerate(G) if g[F.L0] == F.L0]
    lv = F.dist[F.L0]
    conj = {h: [comp(h, comp(s, inv[h])) for s in range(len(G))] for h in stab}
    def components(S1, S2, lins):
        par = {}
        def f(a):
            while par.setdefault(a, a) != a: a = par[a]
            return a
        edges = []
        for (s, xs) in lins:
            for x in xs:
                u, w = (1, x), (2, G[s][x]); edges.append((u, w, s, x))
                par[f(u)] = f(w)
        comps = {}
        for v in list(par): comps.setdefault(f(v), {"V": set(), "E": []})["V"].add(v)
        for e in edges: comps[f(e[0])]["E"].append(e)
        return comps
    def closure(S1, S2, lins):
        while True:
            changed = False
            for c in components(S1, S2, lins).values():
                V = c["V"]
                avail = any((S1 >> x) & 1 if side == 1 else (S2 >> x) & 1 for (side, x) in V)
                if avail or len(c["E"]) >= len(V):
                    for (side, x) in V:
                        if side == 1 and not (S1 >> x) & 1: S1 |= 1 << x; changed = True
                        if side == 2 and not (S2 >> x) & 1: S2 |= 1 << x; changed = True
            if not changed: return S1, S2
    def canon(S1, S2, lins):
        best = None
        for h in stab:
            g = G[h]; cj = conj[h]
            for sw in (0, 1):
                if not sw:
                    a = sum(1 << g[x] for x in range(N) if (S1 >> x) & 1)
                    b = sum(1 << g[x] for x in range(N) if (S2 >> x) & 1)
                    L = tuple(sorted((cj[s], tuple(sorted(g[x] for x in xs))) for (s, xs) in lins))
                else:
                    a = sum(1 << g[x] for x in range(N) if (S2 >> x) & 1)
                    b = sum(1 << g[x] for x in range(N) if (S1 >> x) & 1)
                    L = tuple(sorted((cj[inv[s]], tuple(sorted(g[G[s][x]] for x in xs))) for (s, xs) in lins))
                key = (a, b, L)
                if best is None or key < best: best = key
        return best
    def maxlv(S): return max(lv[x] for x in range(N) if (S >> x) & 1)
    def composites(S1, S2, lins):
        """(twist, u) for every side-1 vertex u and side-2 vertex w joined by a path of >= 3 edges in a locked
        tree component"""
        out = set()
        for c in components(S1, S2, lins).values():
            V = c["V"]
            if any((S1 >> x) & 1 if side == 1 else (S2 >> x) & 1 for (side, x) in V) or len(c["E"]) >= len(V):
                continue
            adj = {}
            for (u, w, s, x) in c["E"]:
                adj.setdefault(u, []).append((w, s)); adj.setdefault(w, []).append((u, s))
            for start in [v for v in V if v[0] == 1]:
                # DFS carrying the composite map from side-1 position of start to the current vertex's position
                stack = [(start, None, F.ident and None, 0)]
                seen = {start}
                st = [(start, None)]   # (vertex, twist so far as a permutation index mapping start pos -> vertex pos)
                tw = {start: None}
                while st:
                    v, t = st.pop()
                    for (nb, s) in adj.get(v, []):
                        if nb in tw: continue
                        if v[0] == 1:   # side-1 -> side-2 along edge with twist s: compose s after t
                            nt = s if t is None else comp(s, t)
                        else:           # side-2 -> side-1: compose s^-1 after t
                            nt = comp(inv[s], t)
                        tw[nb] = nt; st.append((nb, nt))
                        if nb[0] == 2:
                            # path length >= 3 means t was not None at a side-1 vertex other than start
                            if v != start: out.add((nt, start[1]))
        return out
    S1, S2 = closure(1 << F.L0, 1 << F.L0, ())
    layer = {canon(S1, S2, ()): (S1, S2, ())}
    for M in range(1, D + 1):
        new = {}
        for (S1, S2, lins) in layer.values():
            cand = []
            for x in range(N):
                if (S1 >> x) & 1:
                    for x2 in F.nb[x]:
                        if not (S1 >> x2) & 1: cand.append((S1 | 1 << x2, S2, lins))
                if (S2 >> x) & 1:
                    for x2 in F.nb[x]:
                        if not (S2 >> x2) & 1: cand.append((S1, S2 | 1 << x2, lins))
            for li, (s, xs) in enumerate(lins):
                for x in xs:
                    for x2 in F.nb[x]:
                        if x2 not in xs:
                            nl = list(lins); nl[li] = (s, tuple(sorted(xs + (x2,)))); cand.append((S1, S2, tuple(nl)))
            for x0 in range(N):
                if not (S1 >> x0) & 1: continue
                for s in range(len(G)):
                    if (S2 >> G[s][x0]) & 1:
                        for x2 in F.nb[x0]: cand.append((S1, S2, lins + ((s, (x2,)),)))
            for (t, u) in composites(S1, S2, lins):
                for x2 in F.nb[u]: cand.append((S1, S2, lins + ((t, (x2,)),)))
            for (a, b, L) in cand:
                a, b = closure(a, b, L)
                key = canon(a, b, L)
                if key not in new: new[key] = (a, b, L)
        layer = new
        mx = max(maxlv(a) + maxlv(b) for (a, b, L) in layer.values())
        g = any(((a >> F.LT) & 1) and ((b >> F.LT) & 1) for (a, b, L) in layer.values())
        bad = [(a, b, L) for (a, b, L) in layer.values() if maxlv(a) + maxlv(b) > M]
        print("m=%d M=%d states=%d  max(L1+L2)=%d  (I) %s  goal reached: %s" %
              (m, M, len(layer), mx, "holds" if mx <= M else "VIOLATED (%d states)" % len(bad), g), flush=True)
        if bad:
            a, b, L = bad[0]
            print("   example: S1 levels", sorted(lv[x] for x in range(N) if (a >> x) & 1),
                  "S2 levels", sorted(lv[x] for x in range(N) if (b >> x) & 1), "lineages", L, flush=True)

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]))
