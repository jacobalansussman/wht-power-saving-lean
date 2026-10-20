"""Infer the hub's climb rule: replay a canonical word (no copies) under several rules and compare with local_profile."""
import sys, json, gzip
from collections import Counter
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2 * v + R; dims = [len(f) for f in F]
def inside(a, b):
    basis = list(b)
    def red(x):
        for y in sorted(basis, reverse=True):
            x = min(x, x ^ y)
        return x
    # reduce with echelon basis
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    for x in a:
        for z in bb: x = min(x, x ^ z)
        if x: return False
    return True
def run(phases, reads_climb, cen=True):
    cur = list(c['start']); H = Counter(); bad = 0
    def climb(r, f):
        nonlocal bad
        if f is None or cur[r] == f: return
        if dims[f] <= dims[cur[r]] or not inside(F[cur[r]], F[f]): bad += 1; return
        H[dims[f] - dims[cur[r]]] += 1; cur[r] = f
    for ph in phases:
        if ph == 'centres':
            for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
            continue
        for e in c.get(ph) or []:
            climb(e[1], e[-1])
            if reads_climb and e[0] == 'add': climb(e[2], e[-1])
    for r in range(n): climb(r, c['final'][r])
    return dict(sorted(H.items())), bad
prof = {int(k): v for k, v in c['local_profile'].items()}
print('hub profile', dict(sorted(prof.items())), 'cst', sum(dims[f] for f in c['center_frames']))
for name, ph in [('A,cen,P,cut,comp,B', ['A', 'centres', 'precondition', 'cut', 'compensation', 'B']),
                 ('A,cen,P,cut,comp,B,cleanup,tail', ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'cleanup', 'tail']),
                 ('A,cen,P,cut,comp,B,tail,cleanup', ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup'])]:
    for rc in (True, False):
        h, bad = run(ph, rc)
        # centres counted as climbs of their registers; hub profile counts centre blocks separately? show both
        print(f'{name:36s} reads_climb={rc!s:5s} hist={h} rejected_steps={bad} match={h == prof}')
