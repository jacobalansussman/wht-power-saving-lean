#!/usr/bin/env python3
# skeptic_collapsed_edges.py (the skeptic's program) : (1) re-derive the reading's forced pairs and the 7-against-8 count with separate loops;
# (2) test the case the first reading did not treat: maps of the Braverman-He physical graph G (rebuilt from the lists in its appendix) into a frame graph at block size 3 or 4 that send some edges to EQUAL frames (those edges cost no move).
# (3) the exponent arithmetic by the engine's inequality.
import re, sys, collections, math
import json
Q = sys.argv[1]                          # the JSON file written by extract_lists.py
R = set(json.load(open(Q))["deleted_edges"])
C = {102, 163, 164, 341, 386, 406, 476, 513, 862}; J = {103, 349, 470, 517, 603, 606, 674, 716, 861}; B = R - C - J
add = lambda x, y: ((x[0] + y[0]) % 3, (x[1] + y[1]) % 3)
mul = lambda x, y: ((x[0] * y[0] - x[1] * y[1]) % 3, (x[0] * y[1] + x[1] * y[0]) % 3)
Z0, ONE = (0, 0), (1, 0); eta = [(r % 3, r // 3) for r in range(9)]
T = [(Z0, Z0, ONE)] + [(Z0, ONE, eta[r]) for r in range(9)] + [(ONE, eta[r], eta[s]) for r in range(9) for s in range(9)]
def dot(x, y):
    s = Z0
    for a, b in zip(x, y): s = add(s, mul(a, b))
    return s
nbr = [[l for l in range(91) if dot(T[p], T[l]) == Z0] for p in range(91)]
edges = [(i, 91 + nbr[i][r]) for i in range(91) for r in range(10)]
adj = collections.defaultdict(set)
for q, (a, b) in enumerate(edges):
    if q not in R: adj[a].add(b); adj[b].add(a)
nE = sum(len(adj[v]) for v in range(182)) // 2
def bfs(s):
    d = {s: 0}; fr = [s]
    while fr:
        nx = []
        for u in fr:
            for w in adj[u]:
                if w not in d: d[w] = d[u] + 1; nx.append(w)
        fr = nx
    return d
D = [bfs(v) for v in range(182)]
sess = [edges[q] for q in sorted(B)]
print("edges", nE, "sessions", len(sess), "session distances", sorted(collections.Counter(D[s][t] for s, t in sess).items()))
# forced[(x,y)] : x<y lie, in this order or the other, on one shortest route of a session
forced = set()
for s, t in sess:
    for x in range(182):
        for y in range(x + 1, 182):
            k = D[x][y]
            if D[s][x] + k + D[y][t] == 5 or D[s][y] + k + D[x][t] == 5: forced.add((x, y))
F = lambda a, b: (min(a, b), max(a, b)) in forced
tot = collections.Counter(D[x][y] for x in range(182) for y in range(x + 1, 182))
fk = collections.Counter(D[x][y] for x, y in forced)
print("pairs by distance", sorted(tot.items())); print("forced pairs by distance", sorted(fk.items()))
# (1) the 7-against-8 count
viol = collections.Counter()
for x, y in forced:
    if D[x][y] != 3: continue
    for u, v in ((x, y), (y, x)):
        Z = [z for z in adj[u] if D[z][v] == 2]
        ok = all(F(u, z) and F(z, v) for z in Z) and all(F(a, b) for a in Z for b in Z if a < b)
        if len(Z) >= 8: viol[(len(Z), ok)] += 1
print("block size 5, no edge collapsed: forced pairs 3 apart with 8 or 9 first steps (count, all needed pairs forced):", sorted(viol.items()), "total", sum(viol.values()))
# (2a) block size 4: a session route of 5 edges goes to a walk of at most 5 steps between frames 4 apart, so a forced
#      pair at G-distance k lands on frames at least k-1 apart.  Two vertices 2 apart can then not share a frame.
all2 = fk[2] == tot[2]
print("block size 4: every pair 2 apart is forced:", all2, "=> the collapsed edges form a matching, at most 91; moves at least", nE - 91, "against 157*4 =", 628)
# (2b) block size 3: a forced pair at G-distance k lands on frames at least k-2 apart: forced pairs 3, 4, 5 apart never share
#      a frame.  Collapsed edges = edges inside the classes of vertices that share a frame.  A cycle of collapsed edges has
#      every three consecutive edges a-b-c-d with (a,d) NOT forced.  Peel the graph of two-edge paths.
nodes = [(a, b, c) for b in range(182) for a in adj[b] for c in adj[b] if a != c]
out = {n: [(n[1], n[2], d) for d in adj[n[2]] if d != n[1] and not F(n[0], d)] for n in nodes}
n3 = sum(1 for n in nodes for d in adj[n[2]] if d != n[1]) // 2
free3 = sum(len(v) for v in out.values()) // 2
print("block size 3: paths of 3 edges", n3, "of which with unforced ends (could be collapsed whole)", free3)
alive = set(nodes); changed = True
while changed:
    changed = False
    for n in list(alive):
        if not any(o in alive for o in out[n]): alive.discard(n); changed = True
cyc_edges = {(min(a, b), max(a, b)) for (a, b, c) in alive} | {(min(b, c), max(b, c)) for (a, b, c) in alive}
print("  two-edge paths that survive the peeling (a collapsed cycle needs them):", len(alive), "; edges on them:", len(cyc_edges))
# crude bound for the number of collapsed edges: a forest on 182 vertices plus every edge that could lie on a collapsed cycle
bound = 181 + len(cyc_edges)
print("  collapsed edges at most 181 + %d = %d ; a saving at block size 3 needs at least %d ; moves at least %d against 157*3 = 471" % (len(cyc_edges), bound, nE - 470, nE - bound))
# (3) exponent arithmetic, engine inequality in the limit of large s: sum over blocks (k/m)^z < W
def x_rank1(moves, W, m): return 1 - math.log(moves / W) / math.log(m)
print("x for 735 rank-1 blocks, W = 157, m = 5: %.4f ; certificate 744/166: %.4f" % (x_rank1(735, 157, 5), x_rank1(744, 166, 5)))
for m in (3, 4, 5): print("  m = %d: x = 0.1 needs moves per array below %.4f (saved fraction %.4f)" % (m, m ** 0.9, 1 - m ** 0.9 / m))
rho = 50 / 785
print("whole-block bound rho/((1-rho) Lbar): Lbar = ln 5: %.4f ; Lbar needed for 0.1: %.3f, i.e. mean block rank about %.2f of 5" % (rho / (1 - rho) / math.log(5), rho / (1 - rho) / 0.1, 5 / math.exp(rho / (1 - rho) / 0.1)))
