#!/usr/bin/env python3
"""g-break: the 56 non-plane dim-1 entering helpers (M copies): gates per phase, and whether each can reach a plane
that holds one of the 56 cold 0->4 scratch helpers (the role they would take).  Usage: python3 -I mcopy_check.py WORD"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
def inspan(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x == 0
B = [basis(f) for f in F]; dims = [len(b) for b in B]
ORDER = ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']
cur = list(c['start']); P_of = {}; path = defaultdict(list)
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f; path[r].append(f)
    if dims[f] == 3: P_of[r] = tuple(B[f])
gates = defaultdict(Counter)
for ph in ORDER:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1]); gates[e[1]][ph] += 1
        if e[0] == 'add': climb(e[2], e[-1]); gates[e[2]][ph] += 1
for r in range(n): climb(r, c['final'][r])
mc = [r for r in range(2*v, n) if dims[c['start'][r]] == 1 and r not in P_of]
scr = [r for r in range(2*v, n) if dims[c['start'][r]] == 0 and dims[c['final'][r]] == 4]
scrP = {P_of[r] for r in scr if r in P_of}
print('M copies', len(mc), ' gates by phase (total):', dict(sum((gates[r] for r in mc), Counter())))
print('M copies with any gate in A:', sum(1 for r in mc if gates[r]['A']))
print('cold 0->4 scratch', len(scr), 'on distinct planes', len(scrP), ' their gates:', dict(sum((gates[r] for r in scr), Counter())))
ok = defaultdict(list)
for r in mc:
    st = B[c['start'][r]]
    for P in scrP:
        if all(inspan(x, list(P)) for x in st) and any(all(inspan(x, B[f]) for x in P) for f in path[r]):
            ok[r].append(P)
match = {}
def aug(u, seen):
    for P in ok[u]:
        if P in seen: continue
        seen.add(P)
        if P not in match or aug(match[P], seen): match[P] = u; return True
    return False
print('M copies reaching a scratch plane on their own chain:', sum(1 for r in mc if ok[r]), ' perfect matching size:', sum(1 for r in mc if aug(r, set())))
print('M-copy dim paths:', Counter(tuple([1] + [dims[f] for f in path[r]]) for r in mc).most_common(5))
