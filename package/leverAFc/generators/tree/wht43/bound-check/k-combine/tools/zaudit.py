"""k-lower: exact rational audit (Shea's rule, toolkit/verification/scalar_audit.py) for the zero checks
GS.zCheck: per batch, sources = gl (content 1), replay A, scatter, B keeping positive / negative accumulations apart;
units = lcm of the denominators per register class; sw = bit length of the largest digit + 1.  x registers are
never written, so ux = us (an x source then decodes to 1 like a slot source).  Then the EXACT packed mirror of the
kernel term (zplan.zCheck) on every batch at the chosen SPar, and at sw - 1 (margin).
usage: zaudit.py <cert.gcert.json.gz> <out.json>"""
import sys, json, gzip, os
from fractions import Fraction as Q
from math import lcm
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import zmirror as zm
from zplan import plan, zCheck

def adds(g):
    assert g[0] == 'out'
    return [(t, g[2], Q(a, b)) for t, a, b in g[3]]

def replay(c, gl, rows):
    v, nr = c['v'], 2 * c['v'] + c['R']
    cls = lambda r: 'x' if r < v else 'y' if r < 2 * v else 's'
    P, N = [dict() for _ in range(nr)], [dict() for _ in range(nr)]
    den, mx = {k: 1 for k in 'xys'}, {k: Q(1) for k in 'xys'}
    for r in gl: P[r][r] = Q(1)
    def add(t, s, a):
        pos, neg = (P[s], N[s]) if a >= 0 else (N[s], P[s])
        for dst, src in ((P[t], pos), (N[t], neg)):
            for key, value in src.items():
                z = dst.get(key, Q(0)) + abs(a) * value
                if z:
                    dst[key] = z
                    k = cls(t)
                    den[k] = lcm(den[k], z.denominator)
                    mx[k] = max(mx[k], z)
    for g in c['A']:
        for t, s, a in adds(g): add(t, s, a)
    for t, row in enumerate(c['scat']['table']):
        for k, a, b in row: add(v + t, c['ret'][k][1], Q(a, b))
    for g in c['B']:
        for t, s, a in adds(g): add(t, s, a)
    g = set(gl)
    bad = []
    for r in rows:
        if r in g: continue
        net = {k: z for k in P[r].keys() | N[r].keys() if (z := P[r].get(k, Q(0)) - N[r].get(k, Q(0)))}
        if net: bad.append(r)
    return den, mx, bad

c0 = json.loads(gzip.open(sys.argv[1], 'rt').read())
pl = plan(c0)
jobs = [('bit', b, val, gl, pl['slots']) for b, val, gl in pl['bits']] + [('chunk', i, None, gl, pl['el']) for i, gl in enumerate(pl['chunks'])]
den, mx, bads = {k: 1 for k in 'xys'}, {k: Q(1) for k in 'xys'}, []
for kind, i, val, gl, rows in jobs:
    d1, m1, bad = replay(c0, gl, rows)
    for k in 'xys':
        den[k] = lcm(den[k], d1[k]); mx[k] = max(mx[k], m1[k])
    if bad: bads.append([kind, i, val, bad[:5]])
U = lcm(den['s'], den['x'])
uy = den['y']
largest = max(mx['s'] * U, mx['y'] * uy, mx['x'] * U)
assert largest.denominator == 1
sw = largest.numerator.bit_length() + 1
d = (2 * c0['v'] + c0['R'] - 1).bit_length()
spar = [d, sw, U, U, uy]
c = zm.lean_data(c0)
res = []
for kind, i, val, gl, rows in jobs:
    r1 = zm.evaluate(zCheck, c, tuple(spar), gl, rows)
    r0 = zm.evaluate(zCheck, c, (d, sw - 1, U, U, uy), gl, rows)
    res.append(dict(kind=kind, i=i, val=val, size=len(gl), value=r1['value'], why=r1.get('why'), adds=r1.get('adds'),
                    seconds=r1['seconds'], margin_value=r0['value']))
out = dict(rational_zero_failures=bads, denominators=den, maximum={k: str(x) for k, x in mx.items()}, spar=spar,
           all_true=all(r['value'] for r in res), margin_all_false=all(not r['margin_value'] for r in res),
           batches=res, plan=dict(el=pl['el'], nb=pl['nb'], chunks=pl['chunks'], bits=[[b, v, g] for b, v, g in pl['bits']]))
json.dump(out, open(sys.argv[2], 'w'))
print(json.dumps({k: out[k] for k in ('rational_zero_failures', 'denominators', 'maximum', 'spar', 'all_true', 'margin_all_false')}, indent=1))
print([ (r['kind'], r['i'], r['val'], r['size'], r['value'], r['margin_value'], r['seconds']) for r in res])
