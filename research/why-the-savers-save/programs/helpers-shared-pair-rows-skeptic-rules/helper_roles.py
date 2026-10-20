#!/usr/bin/env python3
"""v-pair-market-rules: role of each of U8's 504 helpers (CANONICAL word; census climb rule; values over Q with
symbolic sources, helpers start empty = dirt ignored).  For each helper: start dim, final dim, visits a rank-3 frame?,
read by a target at that target's rank-2 stop (own-block pair row)?, number of distinct targets reading it in B,
max source-content size.  Also: target<->target gates in B (the 'borrowed output pairs').
Usage: python3 -I helper_roles.py CANONICAL-WORD.json.gz BLOCKS.json"""
import sys, json, gzip
from fractions import Fraction as F
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']; blk = {s: i for i, x in enumerate(blocks) for s in x}
v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
def inside(a, b):
    for x in a:
        for z in b: x = min(x, x ^ z)
        if x: return False
    return True
FR = [basis(f) for f in c['frames']]; dims = [len(b) for b in FR]
cur = list(c['start']); visits = defaultdict(set); val = [dict() for _ in range(n)]
for s in range(v): val[s] = {s: F(1)}
for r in range(n): visits[r].add(dims[cur[r]])
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(FR[cur[r]], FR[f]): return
    cur[r] = f; visits[r].add(dims[f])
rank2_partner = Counter(); readers = defaultdict(set); maxc = Counter(); tt = []
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add':
            a, b = e[1], e[2]
            if ph == 'B' and v <= a < 2*v:
                if b >= 2*v: readers[b].add(a)
                if dims[cur[a]] < 2 and e[-1] is not None and dims[e[-1]] == 2: rank2_partner[b] += 1
                if v <= b < 2*v: tt.append((a-v, b-v, dims[e[-1]]))
            climb(a, e[-1]); climb(b, e[-1])
            if ph not in ('cut', 'compensation'):
                k = F(e[3], e[4])
                for key, x in list(val[b].items()):
                    y = val[a].get(key, 0) + k*x
                    if y: val[a][key] = y
                    else: val[a].pop(key, None)
                for r in (a, b): maxc[r] = max(maxc[r], len(val[r]))
        else:
            climb(e[1], e[-1])
            if ph not in ('cut', 'compensation'):
                k = F(e[2], e[3]); val[e[1]] = {key: x*k for key, x in val[e[1]].items()}
for r in range(n): climb(r, c['final'][r])
rows = []
for h in range(2*v, n):
    rows.append(dict(h=h, start=dims[c['start'][h]], final=dims[cur[h]], plane=3 in visits[h],
                     rank2=rank2_partner[h] > 0, nreaders=len(readers[h]), maxc=maxc[h]))
cnt = Counter((r['plane'], r['rank2'], '%d->%d' % (r['start'], r['final'])) for r in rows)
print('helpers by (visits rank 3, partner at a rank-2 stop, start->final dim):')
for k, x in sorted(cnt.items(), key=lambda kv: -kv[1]): print('   ', k, x)
np_ = [r for r in rows if not r['plane']]
print('non-plane helpers:', len(np_), 'of which rank-2-stop pair rows:', sum(r['rank2'] for r in np_))
print('non-plane, not rank-2 partners, by start->final, #target readers, max content:',
      Counter(('%d->%d' % (r['start'], r['final']), r['nreaders'], r['maxc']) for r in np_ if not r['rank2']).most_common(12))
print('rank-2-stop partners by class:', Counter('x' if b < v else ('y' if b < 2*v else 'h') for b in rank2_partner))
print('target<->target gates in B:', len(tt), ' frame dims', Counter(d for _, _, d in tt),
      ' same U8 block:', sum(blk[a] == blk[b] for a, b, _ in tt))
json.dump(dict(rows=rows, tt=tt), open('helper_roles.json', 'w'))
