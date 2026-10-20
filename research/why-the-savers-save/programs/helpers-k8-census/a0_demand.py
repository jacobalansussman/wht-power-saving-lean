#!/usr/bin/env python3
"""g-search a0: the U8 demand and who serves it. Replays the CANONICAL word over Q (hub demand order:
precondition, A, B, cleanup; centres left out, as in u8_nonplane.py), records what every target READS (content of the
read register at read time, as a source set), and the label geometry of the demand.
Usage: python3 -I a0_demand.py WORD out.json"""
import sys, json, gzip
from fractions import Fraction as Fr
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
dot = lambda a, b: bin(a & b).count('1') & 1
start = c['start']; F = c['frames']
val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: Fr(1)}
reads = []   # (phase, target, reg, frozenset sources, coef)
treads = Counter(); hh = Counter()
for ph in ['precondition', 'A', 'B', 'cleanup']:
    for e in c.get(ph) or []:
        if e[0] == 'add':
            _, a, b, num, den = e[:5]; k = Fr(num, den)
            if cls(a) == 'y':
                reads.append((ph, a - v, b, cls(b), sorted(val[b]), str(k)))
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            _, a, num, den = e[:4]; k = Fr(num, den); val[a] = {key: x*k for key, x in val[a].items()}
dem = {t: {s: str(x) for s, x in val[v+t].items()} for t in range(v)}
nnz = sum(len(d) for d in dem.values())
print('target nnz', nnz)
# geometry of demand: dot(u_s,u_t)
geo = Counter((dot(ports[s], ports[t])) for t in range(v) for s in dem[t])
print('demand by dot(u_s,u_t):', geo, ' totals of dot0/dot1 pairs:', Counter(dot(ports[s], ports[t]) for t in range(v) for s in range(v)))
print('coef values in demand', Counter(x for d in dem.values() for x in d.values()))
print('demand sizes per target', Counter(len(d) for d in dem.values()))
print('is s==t in demand?', sum(1 for t in range(v) if t in dem[t]))
# target reads by class of register read and size of its content
print('target reads: class/size', Counter((r[3], len(r[4])) for r in reads))
print('target reads by phase/class', Counter((r[0], r[3]) for r in reads))
json.dump(dict(demand=dem, reads=reads, ports=ports), open(sys.argv[2], 'w'))
