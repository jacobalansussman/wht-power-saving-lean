#!/usr/bin/env python3
"""tokgame.py (p6-prove, 15 Oct 2026): exhaustive search of the TOKEN GAME at blocks 1 and 2.

The token game (PROOF.md section 3) is the engine game for TWO inputs, ANY number of scratch arrays, and
single-term contents (every entry of every register is 0 or one term c.g at every moment).  State:
  S1, S2  = positions (Lagrangians) where input 1 / input 2 is held ALONE in some register (copies are free);
  pair nodes (x, sigma, lineage): a register holding (a g1 x1 + b g2 x2), input 1 at position x, input 2 at
            position sigma(x); sigma in Sp(2m,2) is fixed for the lineage; the coefficient vector (a,b) is fixed
            per lineage (generic, so different lineages are independent).
Free closure: at a location (x, sigma) the objects {lineages there} + {e1 if x in S1} + {e2 if sigma x in S2}
  span C^2 as soon as there are two of them; then x joins S1 and sigma x joins S2 (extraction / rank 2).
Paid steps (cost 1): move a single token to a neighbour; move a pair node to (x', sigma) with x' ~ x (old copy
  kept); create a lineage at (x0, sigma) with x0 in S1, sigma x0 in S2 and move it once.
Checks, for every state reachable with <= D paid moves (up to the symmetry Stab(L0) x input swap):
  (I)  max level in S1 + max level in S2 <= M      (level = distance from L0)
  goal: LT in S1 and LT in S2 needs M >= 2m.
"""
import sys, itertools
from frames import Frames

def main(m, D):
    F = Frames(m); N = F.N
    G = F.all_elements(); gi = {g: k for k, g in enumerate(G)}
    inv = [gi[tuple(sorted(range(N), key=lambda i: g[i]))] for g in G]
    comp = lambda a, b: gi[tuple(G[a][G[b][i]] for i in range(N))]          # a after b
    stab = [k for k, g in enumerate(G) if g[F.L0] == F.L0]
    lv = F.dist[F.L0]
    conj = {}
    for h in stab:
        hi = inv[h]
        conj[h] = [comp(h, comp(s, hi)) for s in range(len(G))]
    def closure(S1, S2, lins):
        changed = True
        while changed:
            changed = False
            loc = {}
            for li, L in enumerate(lins):
                for (x, s) in L: loc.setdefault((x, s), set()).add(li)
            for (x, s), ls in loc.items():
                y = G[s][x]
                k = len(ls) + ((S1 >> x) & 1) + ((S2 >> y) & 1)
                if k >= 2 and (not (S1 >> x) & 1 or not (S2 >> y) & 1):
                    S1 |= 1 << x; S2 |= 1 << y; changed = True
        return S1, S2
    def canon(S1, S2, lins):
        best = None
        for h in stab:
            g = G[h]; cj = conj[h]
            for sw in (0, 1):
                if not sw:
                    a = sum(1 << g[x] for x in range(N) if (S1 >> x) & 1)
                    b = sum(1 << g[x] for x in range(N) if (S2 >> x) & 1)
                    L = tuple(sorted(tuple(sorted((g[x], cj[s]) for (x, s) in Ln)) for Ln in lins))
                else:
                    a = sum(1 << g[x] for x in range(N) if (S2 >> x) & 1)
                    b = sum(1 << g[x] for x in range(N) if (S1 >> x) & 1)
                    L = tuple(sorted(tuple(sorted((g[G[s][x]], cj[inv[s]]) for (x, s) in Ln)) for Ln in lins))
                key = (a, b, L)
                if best is None or key < best: best = key
        return best
    def maxlv(S):
        return max(lv[x] for x in range(N) if (S >> x) & 1)
    S1, S2 = closure(1 << F.L0, 1 << F.L0, ())
    layer = {canon(S1, S2, ()): (S1, S2, ())}
    worst = []; goal_at = None
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
            for li, L in enumerate(lins):
                s = next(iter(L))[1]
                for (x, _) in L:
                    for x2 in F.nb[x]:
                        if (x2, s) not in L:
                            nl = list(lins); nl[li] = tuple(sorted(L + ((x2, s),)))
                            cand.append((S1, S2, tuple(nl)))
            for x0 in range(N):
                if not (S1 >> x0) & 1: continue
                for s in range(len(G)):
                    if (S2 >> G[s][x0]) & 1:
                        for x2 in F.nb[x0]:
                            cand.append((S1, S2, lins + (((x0, s), (x2, s)),)))
            for (a, b, L) in cand:
                a, b = closure(a, b, L)
                key = canon(a, b, L)
                if key not in new: new[key] = (a, b, L)
        layer = new
        mx = max(maxlv(a) + maxlv(b) for (a, b, L) in layer.values())
        g = any(((a >> F.LT) & 1) and ((b >> F.LT) & 1) for (a, b, L) in layer.values())
        # goal up to right multiplication of each column by a free letter: both inputs at level m
        g2 = any(maxlv(a) == m and maxlv(b) == m for (a, b, L) in layer.values())
        print("m=%d M=%d states=%d  max(L1+L2)=%d  (I) %s   goal LT/LT reached: %s   both at level m: %s"
              % (m, M, len(layer), mx, "holds" if mx <= M else "VIOLATED", g, g2), flush=True)
        if mx > M: worst.append(M)
    return worst

if __name__ == "__main__":
    m = int(sys.argv[1]); D = int(sys.argv[2])
    main(m, D)
