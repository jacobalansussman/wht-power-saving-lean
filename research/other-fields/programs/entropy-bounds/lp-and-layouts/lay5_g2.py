# lay5_g2.py -- the Gamma_2 layout with 5 sessions where cuts fail (cut LP 9 < 10 = routing; Theorem 2W: truth 10).
# Contract the 9 non-endpoint frames into one cell and see whether the contraction keeps the routing value and the
# cut gap, and how many LP elements it needs.
import itertools, collections, json, sys
import numpy as np
from scipy.optimize import linprog
import dpg
from layout_cut2 import cut_lp
verts, idx, adj, dist = dpg.build(2)
n = len(verts)
D = [(0, 7), (3, 1), (4, 10), (0, 10), (10, 7)]
arcs = [(u, v) for u in range(n) for v in adj[u]]
print("full-graph directed cut LP:", cut_lp(n, arcs, D)[0], " routing", 2 * len(D))
ends = sorted({x for p in D for x in p})
print("endpoints", ends, " edges among endpoints:", [(a, b) for a in ends for b in ends if a < b and b in adj[a]])
def quotient(part):
    k = max(part) + 1
    qa = sorted({(part[u], part[v]) for (u, v) in arcs if part[u] != part[v]})
    sess = [(part[s], part[t]) for (s, t) in D]
    nb = collections.defaultdict(set)
    for (a, b) in qa: nb[a].add(b)
    def qd(s, t):
        fr = {s}; seen = {s}; d = 0
        while t not in fr:
            fr = {y for x in fr for y in nb[x] if y not in seen}; seen |= fr; d += 1
        return d
    R = sum(qd(s, t) for (s, t) in sess)
    rows = []; rhs = []
    for q in range(1, 1 << k):
        Q = {t for t in range(k) if q >> t & 1}
        need = sum(1 for (s, t) in sess if s not in Q and t in Q)
        if need: rows.append([-1.0 if (a not in Q and b in Q) else 0.0 for (a, b) in qa]); rhs.append(-need)
    cut = linprog(np.ones(len(qa)), A_ub=np.array(rows), b_ub=np.array(rhs), bounds=(0, None), method='highs').fun
    return k, qa, sess, R, cut
# one cell for all non-endpoints
part = [ends.index(v) if v in ends else len(ends) for v in range(n)]
k, qa, sess, R, cut = quotient(part)
print(f"7-cell contraction: arcs {len(qa)}, sessions {len(sess)}, ground {len(qa)+len(sess)}, routing {R}, cut LP {cut}")
# try splitting the rest into 2 cells to see if fewer arcs / still a gap
best = None
others = [v for v in range(n) if v not in ends]
for mask in range(1, 1 << (len(others) - 1)):
    part2 = part[:]
    for i, v in enumerate(others):
        if mask >> i & 1: part2[v] = len(ends) + 1
    k2, qa2, s2, R2, c2 = quotient(part2)
    if R2 == R and c2 < R - 1e-9:
        if best is None or len(qa2) < best[0]: best = (len(qa2), c2, part2)
print("best 8-cell contraction keeping routing 10 with a cut gap:", (best[0], best[1]) if best else None)
json.dump(dict(D=D, part7=part, part8=best[2] if best else None), open('lay5_g2.json', 'w'))
