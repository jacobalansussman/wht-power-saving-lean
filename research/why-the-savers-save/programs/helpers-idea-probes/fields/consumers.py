#!/usr/bin/env python3
"""g-fields: who consumes each packet plane P of the U8 word, and how many copies the consumption geometry needs.
For each rank-3 plane P: its 3 visitors; for each visitor holding P afterwards, the reads of it (p += c*q) after its visit,
with reader class, reader frame dim, time. Then a Dilworth-style chain cover: consumption events (time, gate frame) that
one register could serve in a nested chain (t1 < t2 and f1 <= f2), computed exactly (bipartite matching).
Usage: python3 -I consumers.py  (needs timeline.json)"""
import json
from collections import Counter, defaultdict
T = json.load(open('timeline.json')); ev = {int(q): E for q, E in T['ev'].items()}; cls = {int(q): k for q, k in T['cls'].items()}
dims = T['dims']; Bf = [tuple(b) for b in T['frames']]; start = T['start']
def inside(a, b):
    for x in Bf[a]:
        y = x
        for z in Bf[b]: y = min(y, y ^ z)
        if y: return False
    return True
# replay frames per register per event (same climb rule)
cur = {q: start[q] for q in cls}
visit3 = defaultdict(list)    # plane-basis -> [(reg, t)]
held = {}                      # reg -> (plane, t of visit)
for q in cls:
    if dims[start[q]] == 3: visit3[Bf[start[q]]].append((q, -1))
allev = sorted((t, q, ph, role, p, f) for q, E in ev.items() for (t, ph, role, p, f) in E)
for t, q, ph, role, p, f in allev:
    if f is None or cur[q] == f: continue
    if dims[f] <= dims[cur[q]] or not inside(cur[q], f): continue
    cur[q] = f
    if dims[f] == 3: visit3[Bf[f]].append((q, t))
planes = sorted(visit3)
print('planes', len(planes), 'visitor counts', Counter(len(v) for v in visit3.values()))
res = []
copies_needed = Counter(); consumers_hist = Counter(); vis_cls = Counter()
for P in planes:
    vs = visit3[P]; vis_cls[tuple(sorted(cls[q] for q, _ in vs))] += 1
    cons = []
    for q, tv in vs:
        for (t, ph, role, p, f) in ev[q]:
            if t > tv and role == 'r' and cls[p] in ('y',):
                cons.append((t, f, q, p))
    consumers_hist[(len(cons), len({c[3] for c in cons}))] += 1
    # chain cover of consumption events by frames (time order + inclusion)
    E = sorted(cons)
    n = len(E); adj = [[j for j in range(n) if E[i][0] < E[j][0] and E[i][1] is not None and E[j][1] is not None and inside(E[i][1], E[j][1])] for i in range(n)]
    match = [-1]*n
    def aug(i, seen):
        for j in adj[i]:
            if j in seen: continue
            seen.add(j)
            if match[j] < 0 or aug(match[j], seen): match[j] = i; return True
        return False
    m = sum(1 for i in range(n) if aug(i, set()))
    copies_needed[n - m] += 1
    # static (time-free) cover: inclusion only
    adj2 = [[j for j in range(n) if i != j and E[i][1] is not None and E[j][1] is not None and inside(E[i][1], E[j][1]) and (E[i][1] != E[j][1] or i < j)] for i in range(n)]
    match = [-1]*n
    def aug2(i, seen):
        for j in adj2[i]:
            if j in seen: continue
            seen.add(j)
            if match[j] < 0 or aug2(match[j], seen): match[j] = i; return True
        return False
    m2 = sum(1 for i in range(n) if aug2(i, set()))
    res.append(dict(plane=list(P), visitors=[(q, cls[q], tv) for q, tv in vs], consumers=[(t, dims[f] if f is not None else None, cls[q], p) for t, f, q, p in E], chain_cover_timed=n - m, chain_cover_static=n - m2))
print('visitor classes per plane', dict(vis_cls))
print('(target reads of plane holders after visit, distinct targets) histogram', dict(consumers_hist))
print('timed chain cover (copies the read geometry needs)', dict(copies_needed))
print('static chain cover', dict(Counter(r['chain_cover_static'] for r in res)))
for r in res[:3]: print(json.dumps(r))
json.dump(res, open('consumers.json', 'w'))
