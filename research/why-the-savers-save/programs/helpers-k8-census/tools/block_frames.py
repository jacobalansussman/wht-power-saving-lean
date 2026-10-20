#!/usr/bin/env python3
"""k8-pairs: per-block helper count in the visit language of pm.py: a register 'serves' block b if it visits a
rank-2 (line) or rank-3 (plane) frame lying inside the span of b's 8 labels (census climb rule).  Lines/planes of
a frame inside the union of TWO block spans but not one are reported as cross (leaver cells).
Usage: python3 -I block_frames.py WORD U8_BLOCKS out.json"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
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
spans = [basis([ports[s] for s in x]) for x in blocks]
def home(f):
    return tuple(i for i, sp in enumerate(spans) if all(inspan(x, sp) for x in B[f]))
cur = list(c['start']); vis = defaultdict(list)
for r in range(n): vis[r].append(cur[r])
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not all(inspan(x, B[f]) for x in B[cur[r]]): return
    cur[r] = f; vis[r].append(f)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
for r in range(n): climb(r, c['final'][r])
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
per = defaultdict(Counter); kinds = Counter(); cross = Counter(); multi = Counter()
hb = {}
for r in range(n):
    homes = []
    for f in vis[r]:
        if dims[f] in (2, 3):
            hm = home(f)
            homes.append((dims[f], hm))
    bs = {h for d, hm in homes for h in hm}
    if cls(r) == 'h':
        hb[r] = sorted(bs)
        path = tuple('L' if d == 2 else 'P' for d, hm in homes if hm)
        xpath = tuple('L' if d == 2 else 'P' for d, hm in homes if not hm)
        if len(bs) == 1: per[min(bs)][''.join(path) or '-'] += 1
        elif len(bs) > 1: multi[len(bs)] += 1
        if xpath: cross[''.join(xpath)] += 1
        kinds[(''.join(path) or '-', ''.join(xpath) or '-')] += 1
tot = {b: sum(per[b].values()) for b in range(len(blocks))}
print('helpers per block (in-block line/plane visits, single block):', tot, 'sum', sum(tot.values()))
print('per block kinds:', {b: dict(per[b]) for b in range(len(blocks))})
print('helpers with in-block frames in >1 block:', dict(multi))
print('helper kinds (in-block path, cross-frame path):', dict(kinds))
# cross frames: rank 2/3 frames with no single home block
xs = Counter()
for r in range(n):
    for f in vis[r]:
        if dims[f] in (2, 3) and not home(f):
            S = tuple(s for s in range(v) if inspan(ports[s], B[f]))
            xs[(dims[f], len(S), cls(r))] += 1
print('visits to rank-2/3 frames outside every single block span (dim, #labels inside, class):', dict(xs))
json.dump(dict(per_block=tot, kinds={b: dict(per[b]) for b in range(len(blocks))}, helper_blocks=hb), open(sys.argv[3], 'w'), indent=0)
