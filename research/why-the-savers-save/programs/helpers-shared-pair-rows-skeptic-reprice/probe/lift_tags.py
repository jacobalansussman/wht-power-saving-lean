#!/usr/bin/env python3
"""v-pair-market-reprice: which register classes carry a reader's own-plane directions into another target's frame
when every target sits at P_own(t) from its rank-2 stop (pair-market design) and U8's other gates stay as they are.
Propagation is run with and without the 112 entering (M) copies and the target-target gates, so the face side
(plane helpers + sources, which the model keeps at 28 per block with q = 4) is tested alone.
Usage: python3 -I lift_tags.py WORD DEMAND BLOCKS"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
d = json.load(open(sys.argv[2])); U8 = json.load(open(sys.argv[3]))['blocks']
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']; h = c['h']
blk = {s: i for i, x in enumerate(U8) for s in x}
def B_(gens):
    piv = {}
    for a in gens:
        while a:
            bit = a.bit_length() - 1
            if bit in piv: a ^= piv[bit]
            else: piv[bit] = a; break
    for k in sorted(piv):
        for j in piv:
            if j > k and piv[j] >> k & 1: piv[j] ^= piv[k]
    return tuple(sorted(piv.values(), reverse=True))
def cont(a, b): return B_(list(a) + list(b)) == B_(b)
perp = {t: B_([x for x in range(1, 1 << h) if bin(x & ports[t]).count('1') % 2 == 0]) for t in range(v)}
dem = {int(t): {int(s): x for s, x in dd.items()} for t, dd in d['demand'].items()}
sub_readers = defaultdict(list)
for t in range(v):
    per = defaultdict(list)
    for s in dem[t]: per[blk[s]].append(s)
    for b, ss in per.items(): sub_readers[tuple(sorted(ss))].append(t)
own = {t: B_([ports[s] for s in ss]) for ss, ts in sub_readers.items() if len(ts) == 4 for t in ts}
entering = {r for r in range(2*v, n) if F[c['start'][r]]}
print('entering helpers (M copies):', len(entering))
# pairwise: is reader s's own plane inside target t's u^perp?
ok = sum(cont(own[s], perp[t]) for s in range(v) for t in range(v) if s != t)
print('ordered target pairs (s,t) with P_own(s) inside u_t^perp:', ok, 'of', v*(v-1))
def run(skipM, skipYY):
    tags = [set() for _ in range(n)]; arrived = {}
    r2 = set(); cur = [B_(F[f]) for f in c['start']]
    for ph in ['B', 'tail']:
        for gi, e in enumerate(c.get(ph) or []):
            f = e[-1]
            if f is None: continue
            regs = list(e[1:3] if e[0] == 'add' else e[1:2])
            for r in regs:
                if v <= r < 2*v and len(F[f]) >= 2 and r not in r2: r2.add(r); tags[r].add(r - v)
            if len(regs) == 2:
                a, b = regs
                if skipM and (a in entering or b in entering): continue
                if skipYY and v <= a < 2*v and v <= b < 2*v: continue
                for x, y in ((a, b), (b, a)):
                    for tg in tags[y] - tags[x]: arrived[(x, tg)] = (ph, gi, y)
                u = tags[a] | tags[b]; tags[a] = set(u); tags[b] = set(u)
    bad = {}
    for t in range(v):
        bt = [s for s in tags[v + t] if not cont(own[s], perp[t])]
        if bt: bad[t] = bt
    cls = Counter()
    for t, bt in bad.items():
        for s in bt:
            ph, gi, y = arrived[(v + t, s)]
            cls['target' if v <= y < 2*v else 'M copy' if y in entering else 'source' if y < v else 'plane/cold helper'] += 1
    return bad, cls
for skipM, skipYY in ((False, False), (True, False), (True, True)):
    bad, cls = run(skipM, skipYY)
    print(f'skip M copies={skipM} skip target-target={skipYY}: targets outside u^perp {len(bad)}, bad tags {sum(map(len, bad.values()))}, last hop {dict(cls)}')
