#!/usr/bin/env python3
"""v-an2-word: are the complex Coxeter-Todd lines' products re-phasable to powers of i?  Rephasing-invariant test:
the triple products <a,b><b,c><c,a> over mutually non-orthogonal triples.  If one has a phase outside {1, i, -1, -i},
no choice of line representatives (any complex phases) makes every product a power of i times a positive real.
Own construction (x_i congruent mod theta, sum = 0 mod 3).  usage: python3 -I vct12_phase.py"""
import itertools, cmath, math
from collections import Counter
def enorm(z): a, b = z; return a*a - a*b + b*b
def emul(z, u): a, b = z; c, d = u; return (a*c - b*d, a*d + b*c - b*d)
def econj(z): a, b = z; return (a - b, -b)
def eadd(z, u): return (z[0]+u[0], z[1]+u[1])
TH = (1, 2)
def div(z, d):
    p = emul(z, econj(d)); n = enorm(d); return p[0] % n == 0 and p[1] % n == 0
small = [(a, b) for a in range(-4, 5) for b in range(-4, 5) if enorm((a, b)) <= 6]
byn = {}
for z in small: byn.setdefault(enorm(z), []).append(z)
V = []
def rec(pre, bud):
    if len(pre) == 6:
        if bud == 0:
            s = (0, 0)
            for z in pre: s = eadd(s, z)
            if div(s, (3, 0)) and all(div((pre[i][0]-pre[0][0], pre[i][1]-pre[0][1]), TH) for i in range(1, 6)):
                V.append(tuple(pre))
        return
    for k, zs in byn.items():
        if k <= bud:
            for z in zs: rec(pre + [z], bud - k)
rec([], 6)
units = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, -1), (1, 1)]
seen, C = set(), []
for x in V:
    if x in seen: continue
    for u in units: seen.add(tuple(emul(u, z) for z in x))
    C.append(x)
w = cmath.exp(2j * math.pi / 3)
def cx(z): return z[0] + z[1] * w
def herm(x, y): return sum(cx(a) * cx(b).conjugate() for a, b in zip(x, y))
n = len(C)
Hm = [[herm(C[i], C[j]) for j in range(n)] for i in range(n)]
ph = Counter()
for i, j, k in itertools.combinations(range(n), 3):
    if abs(Hm[i][j]) > 1e-9 and abs(Hm[j][k]) > 1e-9 and abs(Hm[k][i]) > 1e-9:
        t = Hm[i][j] * Hm[j][k] * Hm[k][i]
        ph[round(math.degrees(cmath.phase(t))) % 360] += 1
print('complex lines', n, '; phases (degrees) of triple products over mutually non-orthogonal triples:', dict(sorted(ph.items())))
bad = {d: c for d, c in ph.items() if d % 90}
print('triples with a phase that is not a multiple of 90 degrees:', sum(bad.values()), '->',
      'NOT re-phasable to powers of i' if bad else 're-phasable test inconclusive')
