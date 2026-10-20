#!/usr/bin/env python3
# placement_test.py : can the 182-vertex physical graph G of arXiv:2610.10108v1 (as rebuilt from the lists in its appendix) be mapped
# into a frame graph X_m (dual polar graph of Sp(2m,2)) so that every physical edge goes to a pair of adjacent frames and
# the two ends of every session go to frames at distance 5?   (Sessions are at G-distance exactly 5: rebuild_and_field.out.)
# Necessary condition tested here [derived in the section README]:
#   such a map sends every G-path of length 5 between the ends of a session to a shortest path of X_m, so all pairs on it keep
#   their distance ("forced pairs").  In X_m, for frames U, V at distance 3 exactly 2^3 - 1 = 7 neighbours of U are at
#   distance 2 from V (frames.out: (2^2 - 1) = 3 common neighbours at distance 2; the general count is the textbook c_3).
#   So: a forced pair (x, y) at distance 3 with 8 or more neighbours z of x at G-distance 2 from y whose images are forced
#   to be pairwise different (each pair z, z' itself a forced pair at distance 2) rules the map out.
import re, sys, collections, itertools
import json
Q = sys.argv[1]                          # the JSON file written by extract_lists.py
R = set(json.load(open(Q))["deleted_edges"])
C = {102, 163, 164, 341, 386, 406, 476, 513, 862}; J = {103, 349, 470, 517, 603, 606, 674, 716, 861}; B = R - C - J
p3 = 3
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
sessions = [edges[q] for q in sorted(B)]
forced = collections.defaultdict(set)            # forced[(x,y)] = set of session indices that force the pair (unordered key x<y)
for si, (s, t) in enumerate(sessions):
    assert D[s][t] == 5
    on = [v for v in range(182) if D[s][v] + D[v][t] == 5]
    for x in on:
        for y in on:
            if x < y and (D[s][x] + D[x][y] + D[y][t] == 5 or D[s][y] + D[y][x] + D[x][t] == 5):
                forced[(x, y)].add(si)
fk = collections.Counter(D[x][y] for (x, y) in forced)
tot = collections.Counter(D[x][y] for x in range(182) for y in range(x + 1, 182))
print("unordered pairs of G by distance:", sorted(tot.items()))
print("of these, forced to keep their distance under any such map:", sorted(fk.items()))
isf = lambda a, b: (min(a, b), max(a, b)) in forced
viol = 0; n3 = 0; best = collections.Counter(); sizeZ = collections.Counter(); examples = []
for (x, y) in list(forced):
    if D[x][y] != 3: continue
    for (u, v) in ((x, y), (y, x)):
        n3 += 1
        Z = [z for z in adj[u] if D[z][v] == 2]
        sizeZ[len(Z)] += 1
        mx = 0
        for k in range(len(Z), 0, -1):
            if any(all(isf(a, b) for a, b in itertools.combinations(S, 2)) for S in itertools.combinations(Z, k)):
                mx = k; break
        best[mx] += 1
        if mx >= 8:
            viol += 1
            if len(examples) < 3: examples.append((u, v, sorted(Z)))
print("forced pairs at distance 3, each taken from both ends:", n3)
print("  number of first steps |Z| (size: count):", sorted(sizeZ.items()))
print("  largest set of first steps with pairwise different images forced (size: count):", sorted(best.items()))
print("  cases with 8 or more, against the 7 a frame graph allows:", viol)
print("  examples (x, y, first steps; vertex ids: point i = i, line j = 91 + j):", examples)
# how many sessions are touched: a session is 'hit' if it forces the pair (x,y) of some violating case
bad_sessions = set()
for (x, y) in forced:
    if D[x][y] != 3: continue
    for (u, v) in ((x, y), (y, x)):
        Z = [z for z in adj[u] if D[z][v] == 2]
        if len(Z) >= 8 and any(all(isf(a, b) for a, b in itertools.combinations(S, 2)) for S in itertools.combinations(Z, 8)):
            bad_sessions |= forced[(x, y)]
print("sessions that force at least one violating pair:", len(bad_sessions), "of", len(sessions), " (the code keeps a gain only while fewer than 10 sessions are dropped: 735 < 5*(157-k))")
