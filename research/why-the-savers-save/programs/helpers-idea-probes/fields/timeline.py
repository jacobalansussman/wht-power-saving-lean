#!/usr/bin/env python3
"""g-fields: per-register gate timeline of the U8 word (time index over the replay order used by u8_census.py).
Usage: python3 -I timeline.py WORD.json.gz  -> timeline.json (per register: list of [t, phase, role(w/r), partner, frame])"""
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
B = [basis(f) for f in F]; dims = [len(b) for b in B]
cen = set(c['centers'])
def cls(q):
    if q in cen: return 'C'
    if q < v: return 'x'
    if q < 2*v: return 'y'
    h = q - 2*v
    return 'fac' if h < 224 else ('face' if h < 392 else 'M')
ev = defaultdict(list); t = 0
order = ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']
for ph in order:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): ev[r].append([t, ph, 'cen', None, f]); t += 1
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            ev[e[1]].append([t, ph, 'w', e[2], e[-1]]); ev[e[2]].append([t, ph, 'r', e[1], e[-1]])
        else:
            ev[e[1]].append([t, ph, 's', None, e[-1]])
        t += 1
json.dump(dict(T=t, cls={q: cls(q) for q in range(n)}, ev=ev, dims=dims, frames=B, start=c['start'], final=c['final'], ports=c['ports']), open('timeline.json', 'w'))
# summary: per class, phases touched, first/last time
for k in ['x', 'y', 'C', 'fac', 'face', 'M']:
    regs = [q for q in range(n) if cls(q) == k]
    ph = Counter(p for q in regs for (_, p, *_r) in ev[q])
    fl = [(ev[q][0][0], ev[q][-1][0]) for q in regs if ev[q]]
    print(k, len(regs), 'phases', dict(ph), 'first t range', min(a for a, b in fl), max(a for a, b in fl), 'last t range', min(b for a, b in fl), max(b for a, b in fl))
print('T', t, 'phase boundaries:')
t = 0
for ph in order:
    L = 8 if ph == 'centres' else len(c.get(ph) or [])
    print(' ', ph, t, t + L); t += L
