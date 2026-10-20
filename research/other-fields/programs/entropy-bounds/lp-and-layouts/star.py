# star.py -- the K_{2,3} "star" contraction of Gamma_2 found by quot3 (agglomerative k = 5):
# A = {v}, X_i = the two other points of each of the 3 lines through v (GQ(2,2) picture), B = the 8 frames opposite v.
# Builds the quotient instance from Gamma_2 itself, checks its structure, and computes its routing value and cut LP.
import itertools, collections
import numpy as np
from scipy.optimize import linprog
import dpg

verts, idx, adj, dist = dpg.build(2)
n = len(verts); v = 0
nb = adj[v]
# lines through v: maximal cliques {v, x, y}; pair neighbours that are adjacent to each other
pairs = [(x, y) for x in nb for y in nb if x < y and y in adj[x]]
assert len(pairs) == 3
part = {v: 'A'}
for i, (x, y) in enumerate(pairs): part[x] = part[y] = f'X{i+1}'
for u in range(n):
    if u not in part: part[u] = 'B'
assert sum(1 for u in range(n) if part[u] == 'B') == 8 and all(dist[v][u] == 2 for u in range(n) if part[u] == 'B')
edges = [(i, j) for i in range(n) for j in adj[i] if i < j]
cnt = collections.Counter(tuple(sorted((part[i], part[j]))) for (i, j) in edges)
print("edges between parts:", dict(cnt))
opp = [(s, t) for s in range(n) for t in range(n) if dist[s][t] == 2]
dem = collections.Counter((part[s], part[t]) for (s, t) in opp if part[s] != part[t])
print("session groups (ordered part pairs, counts):", dict(dem), " internal:", sum(1 for (s, t) in opp if part[s] == part[t]))
P = ['A', 'X1', 'X2', 'X3', 'B']
qadj = {p: set() for p in P}
for (a, b) in cnt:
    if a != b: qadj[a].add(b); qadj[b].add(a)
def qd(a, b):
    if a == b: return 0
    if b in qadj[a]: return 1
    return 2
R = sum(c * qd(a, b) for (a, b), c in dem.items())
cross = sum(c for (a, b), c in cnt.items() if a != b)
print(f"routing value on the quotient R = {R}; crossing edges {cross} of {len(edges)}; rho = {(R/len(opp))/(cross/len(edges)):.5f}")
# quotient arcs: one arc group per ordered adjacent part pair; cut LP: min sum c s.t. for every Q, c(in Q) >= demand into Q
arcs = [(a, b) for a in P for b in P if a != b and b in qadj[a]]
rows = []; rhs = []
for q in range(1, 1 << 5):
    Q = {P[t] for t in range(5) if q >> t & 1}
    need = sum(c for (a, b), c in dem.items() if a not in Q and b in Q)
    if need: rows.append([-1.0 if (a not in Q and b in Q) else 0.0 for (a, b) in arcs]); rhs.append(-need)
res = linprog(np.ones(len(arcs)), A_ub=np.array(rows), b_ub=np.array(rhs), bounds=(0, None), method='highs')
print(f"cut LP on the quotient = {res.fun}  -> per session bound {res.fun/len(opp)*len(edges)/cross:.5f}  (routing gives {R}, {R/len(opp)*len(edges)/cross:.5f})")
print("arc groups:", len(arcs), " session groups:", len(dem), " ground set:", len(arcs) + len(dem))
