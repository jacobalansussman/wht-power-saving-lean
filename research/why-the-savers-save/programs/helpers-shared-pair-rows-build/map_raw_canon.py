#!/usr/bin/env python3
"""pm-build: register and frame map between the raw U8 FULL-WORD (audit_local's input) and the native CANONICAL-WORD."""
import sys, json, gzip
from collections import Counter
sys.path.insert(0, '../helpers-half-reads')
raw = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
can = json.loads(gzip.decompress(open(sys.argv[2], 'rb').read()))
def basis(vs):
    piv = {}
    for a in vs:
        while a:
            bit = a.bit_length()-1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(piv[k] for k in sorted(piv, reverse=True))
CF = [basis(f) for f in can['frames']]
rmap = {}; bad = Counter(); fmap = Counter()
for ph in ('A', 'B', 'cut', 'precondition'):
    a, b = raw[ph], can[ph]
    print(ph, len(a), len(b))
    for x, y in zip(a, b):
        if x[0] != y[0]: bad[ph] += 1; continue
        k = 3 if x[0] == 'add' else 2
        for r, s in zip(x[1:k], y[1:k]):
            if rmap.setdefault(r, s) != s: bad['reg'] += 1
        if x[-1] is not None and y[-1] is not None and ph in ('A','B'):
            fmap[(basis(x[-1]), CF[y[-1]])] += 1
print('bad', bad, 'mapped regs', len(rmap))
v = 112
nonid = {r: s for r, s in rmap.items() if r != s}
print('non-identity mapped', len(nonid), Counter(('x' if r < v else 'y' if r < 2*v else 'h') for r in nonid))
# frame transform guess: tau43
tau = lambda x: x ^ (43 if bin(x & 43).count('1') % 2 else 0)
ok = sum(1 for (a, b), k in fmap.items() if basis([tau(z) for z in a]) == b); print('frame pairs', len(fmap), 'tau43 ok', ok)
ok2 = sum(1 for (a, b), k in fmap.items() if a == b); print('identity ok', ok2)
# ports
P = raw['ports']; Q = can['ports']
print('ports map consistent:', all(Q[rmap[s]] == tau(P[s]) or Q[rmap[s]] == P[s] for s in range(v) if s in rmap))
print('sample', [(s, P[s], rmap.get(s), Q[rmap[s]] if s in rmap else None, tau(P[s])) for s in range(5)])
json.dump({str(k): x for k, x in sorted(rmap.items())}, open('raw2canon_regs.json', 'w'))
