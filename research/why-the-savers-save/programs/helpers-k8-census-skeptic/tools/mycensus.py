#!/usr/bin/env python3
"""v-k8-0 independent census: for every helper, the set of frames it is ever at (start, final, every event's
frame); classify by in-block rank-2/3 frames (frame inside ONE block span of u8_blocks.json).
Also: target<->target gates (both registers are targets y), and which own-block pair rows targets read.
Usage: python3 -I mycensus.py WORD.json.gz u8_blocks.json out.json"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
blocks = json.load(open(sys.argv[2]))['blocks']
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def red(x, bb):
    for z in bb: x = min(x, x ^ z)
    return x
def basis(b):
    bb = []
    for y in b:
        y = red(y, bb)
        if y: bb.append(y); bb.sort(reverse=True)
    return bb
B = [basis(f) for f in F]; dims = [len(b) for b in B]
spans = [basis([ports[s] for s in x]) for x in blocks]
def home(f): return [i for i, sp in enumerate(spans) if all(red(x, sp) == 0 for x in B[f])]
at = defaultdict(set)
for r in range(n): at[r].add(c['start'][r]); at[r].add(c['final'][r])
for r, f in zip(c['centers'], c['center_frames']): at[r].add(f)
tt = Counter(); evs = 0
for ph in ['A', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    for e in c.get(ph) or []:
        f = e[-1]
        regs = [e[1]] + ([e[2]] if e[0] == 'add' else [])
        if f is not None:
            for r in regs: at[r].add(f)
        if e[0] == 'add' and all(v <= r < 2*v for r in regs): tt[(ph, dims[f] if f is not None else None)] += 1
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
perblock = Counter(); noblock = 0; multi = 0; kinds = Counter(); mlines = Counter()
for r in range(2*v, n):
    hs = set(); rk = set()
    for f in at[r]:
        if dims[f] in (2, 3):
            for b in home(f): hs.add(b); rk.add(dims[f])
    vis3 = any(dims[f] == 3 for f in at[r])
    if not hs: noblock += 1; kinds[('none', vis3, dims[c['start'][r]], dims[c['final'][r]])] += 1
    elif len(hs) == 1: perblock[hs.pop()] += 1
    else: multi += 1
out = dict(helpers=R, per_block={b: perblock[b] for b in range(len(blocks))}, sum_inblock=sum(perblock.values()),
           no_inblock_frame=noblock, multi_block=multi, noblock_kinds={str(k): x for k, x in kinds.items()},
           target_target_gates={str(k): x for k, x in tt.items()})
print(json.dumps(out, indent=1)); json.dump(out, open(sys.argv[3], 'w'), indent=1)
