#!/usr/bin/env python3
"""v-pair-market-rules: trace how U8 serves one block's targets (FULL-WORD, literal frames as in the hub checker).
Values: symbolic sources over Q, helpers start empty (dirt ignored), phases A, cut, B in order (cut replayed for values
as the checker does: main = A + cut + B).  For each target of the block: every gate in B touching it, with partner,
partner class, partner content (sources, coefficient), gate frame dim, target frame dim before.
Usage: python3 -I trace_block.py FULL-WORD.json.gz BLOCKS.json DEMAND.json B"""
import sys, json, gzip
from fractions import Fraction as F
w = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']; dem = json.load(open(sys.argv[3]))['demand']; bi = int(sys.argv[4])
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
v = w['v']; R = w['R']; n = 2*v + R
FR = [basis(f) for f in w['frames']]
start = [FR[f] for f in w['start']]
cur = list(start); val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: F(1)}
cls = lambda r: 'x%d' % r if r < v else ('y%d' % (r-v) if r < 2*v else 'h%d' % r)
S8 = set(blocks[bi]); T8 = [v+t for t in blocks[bi]]
log = {r: [] for r in T8}; partners = set()
def fmt(d): return {(cls(k) if k < v else 'D'): str(x) for k, x in d.items()}
for ph in ('A', 'cut', 'compensation', 'B'):
    if ph == 'B':
        # centre / target levels as the checker sets them
        for q, U in zip(w['centers'], w['center_frames']):
            if len(FR[U]) > len(cur[q]): cur[q] = FR[U]
    for i, e in enumerate(w[ph]):
        regs = e[1:3] if e[0] == 'add' else e[1:2]
        U = FR[e[-1]] if e[-1] is not None else None
        if ph == 'B':
            for r in regs:
                if r in log:
                    other = [x for x in regs if x != r]
                    o = other[0] if other else None
                    log[r].append(dict(i=i, role='reads' if r == e[1] else 'is read', op=e[0], k=str(F(e[3], e[4])) if e[0]=='add' else str(F(e[2], e[3])),
                                       partner=cls(o) if o is not None else None, pstart=len(start[o]) if o is not None else None,
                                       pcont=fmt(val[o]) if o is not None else None, pframe=len(cur[o]) if o is not None else None,
                                       gate=len(U), tframe=len(cur[r])))
                    if o is not None: partners.add(o)
        if ph in ('cut', 'compensation'): continue
        if e[0] == 'add':
            a, b = e[1], e[2]; k = F(e[3], e[4])
            for key, x in list(val[b].items()):
                y = val[a].get(key, 0) + k*x
                if y: val[a][key] = y
                else: val[a].pop(key, None)
        else:
            a = e[1]; k = F(e[2], e[3]); val[a] = {key: x*k for key, x in val[a].items()}
        if ph not in ('cut',) and U is not None:
            for r in regs:
                if len(U) > len(cur[r]): cur[r] = U
for r in T8:
    t = r - v
    own = {int(s): x for s, x in dem[str(t)].items() if int(s) in S8}
    print('TARGET', cls(r), 'own-block need', own)
    for g in log[r]:
        if g['op'] == 'add': print('   ', g['i'], g['role'], g['k'], g['partner'], 'pstart', g['pstart'], 'pframe', g['pframe'], 'gate', g['gate'], 'tframe', g['tframe'], 'own:', {k: x for k, x in (g['pcont'] or {}).items() if k != 'D' and int(k[1:]) in S8 and k[0]=='x'}, 'n', len(g['pcont'] or {}))
    print('    ... total B gates', len(log[r]), ' gate dims', [g['gate'] for g in log[r]])
