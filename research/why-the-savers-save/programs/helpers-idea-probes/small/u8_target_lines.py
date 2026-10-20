#!/usr/bin/env python3
"""g-small: is each U8 target's existing rank-2 frame a 2-label line of a packet plane P inside its cap u_t^perp?
(where a 'pre-merge half read' could be taken without a new target climb), and the max matching targets -> planes
through such lines (each target has one rank-2 frame, so it can take at most one half read).
Usage: python3 -I u8_target_lines.py CANONICAL-WORD.json.gz"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
dot = lambda a, b: bin(a & b).count('1') & 1
B = [basis(f) for f in F]; dims = [len(b) for b in B]
def sub(a, b): return all(inspan(x, B[b]) for x in B[a])
cur = list(c['start']); chain = defaultdict(list)
for r in range(n): chain[r].append(cur[r])
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not sub(cur[r], f): return
    cur[r] = f; chain[r].append(f)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
planes = sorted({tuple(B[f]) for r in range(n) for f in chain[r] if dims[f] == 3})
hits = Counter(); match_edges = defaultdict(list)
for t in range(v):
    y = v + t
    L = [f for f in chain[y] if dims[f] == 2]
    if not L: hits['no rank-2 frame'] += 1; continue
    Lb = B[L[0]]
    labs = [s for s in range(v) if inspan(ports[s], Lb)]
    inP = [i for i, P in enumerate(planes) if all(inspan(x, list(P)) for x in Lb) and all(not dot(x, ports[t]) for x in P)]
    hits[(len(labs), len(inP) > 0)] += 1
    match_edges[t] = inP
print('targets by (#port labels on their rank-2 frame, frame lies in a packet plane inside the cap):', dict(hits))
# matching targets -> distinct planes
m = {}
def aug(t, seen):
    for p in match_edges[t]:
        if p in seen: continue
        seen.add(p)
        if p not in m or aug(m[p], seen): m[p] = t; return True
    return False
print('max matching targets->planes via their existing rank-2 line:', sum(aug(t, set()) for t in match_edges))
