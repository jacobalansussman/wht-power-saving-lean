#!/usr/bin/env python3
# first_steps_and_drops.py : (a) in the frame graphs at block sizes 3 and 4, for frames U, V at distance 3, count the neighbours of U at
# distance 2 from V (the number the embedding test compares with);  (b) greedy: drop the 9 sessions involved in most violations
# and count the violations left (9 is the most the code can lose and still beat 5 moves per session).
import sys, re, collections, itertools
import os
HERE = os.path.dirname(os.path.abspath(__file__)) + "/"
src = open(HERE + "frame_graphs.py").read().split("for m in (2, 3, 4):")[0]
exec(src)
for m in (3, 4):
    order, adj = build(m); N = len(order)
    def bfs(s):
        d = {s: 0}; fr = [s]
        while fr:
            nx = []
            for u in fr:
                for w in adj[u]:
                    if w not in d: d[w] = d[u] + 1; nx.append(w)
            fr = nx
        return d
    d0 = bfs(0); cnt = collections.Counter(); c4 = collections.Counter()
    for v in range(N):
        if d0[v] == 3:
            dv = bfs(v); cnt[sum(1 for z in adj[0] if dv[z] == 2)] += 1
        if d0[v] == 4:
            dv = bfs(v); c4[sum(1 for z in adj[0] if dv[z] == 3)] += 1
    print(f"block size {m}: neighbours of U at distance 2 from V, over all V at distance 3 from U (count: how many V): {sorted(cnt.items())}"
          + (f"; at distance 4: neighbours at distance 3: {sorted(c4.items())}" if c4 else ""))
    sys.stdout.flush()
# ---- (b)
src2 = open(HERE + "placement_test.py").read().split("sessions = [edges[q] for q in sorted(B)]")[0]
sys.argv = [sys.argv[0], sys.argv[1]]
exec(src2.split("import re, sys, collections, itertools")[1])
sessions = [edges[q] for q in sorted(B)]
def violations(alive):
    forced = collections.defaultdict(set)
    for si in alive:
        s, t = sessions[si]
        on = [v for v in range(182) if D[s][v] + D[v][t] == 5]
        for x in on:
            for y in on:
                if x < y and (D[s][x] + D[x][y] + D[y][t] == 5 or D[s][y] + D[y][x] + D[x][t] == 5): forced[(x, y)].add(si)
    isf = lambda a, b: (min(a, b), max(a, b)) in forced
    out = []
    for (x, y) in forced:
        if D[x][y] != 3: continue
        for (u, v) in ((x, y), (y, x)):
            Z = [z for z in adj[u] if D[z][v] == 2]
            if len(Z) < 8: continue
            for S in itertools.combinations(Z, 8):
                if all(isf(a, b) for a, b in itertools.combinations(S, 2)):
                    sup = set(forced[(x, y)])
                    out.append(sup); break
    return out
alive = set(range(len(sessions)))
for step in range(10):
    v = violations(alive)
    print(f"sessions kept {len(alive)}: violating cases {len(v)}")
    if not v or step == 9: break
    c = collections.Counter(si for sup in v for si in sup)
    alive.discard(c.most_common(1)[0][0])
