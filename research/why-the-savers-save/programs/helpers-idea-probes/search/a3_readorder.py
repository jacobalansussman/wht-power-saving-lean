#!/usr/bin/env python3
"""g-search a3: in the U8 climb-rule replay, the frame dim of each target BEFORE each of its reads in phase B,
the dim of the read register, and the dim after (= join). Also: which helpers are read by >1 target, and how the
helper's frame moves between its readers (does each read raise the copy?)."""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return bb
B = [basis(f) for f in F]; dims = [len(b) for b in B]
def inside(a, b):
    for x in B[a]:
        y = x
        for z in B[b]: y = min(y, y ^ z)
        if y: return False
    return True
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start'])
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(cur[r], f): return
    cur[r] = f
first = {}; seq = Counter(); helper_readers = defaultdict(list); raise_on_read = Counter()
tdim_at_B = None
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'B': tdim_at_B = Counter(dims[cur[t]] for t in range(v, 2*v))
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        continue
    for e in c.get(ph) or []:
        if e[0] == 'add' and cls(e[1]) == 'y' and ph == 'B':
            a, b = e[1], e[2]
            before = dims[cur[a]]; hb = dims[cur[b]]
            climb(a, e[-1]); climb(b, e[-1])
            seq[(before, hb, dims[cur[a]], cls(b))] += 1
            helper_readers[b].append(a)
            if dims[cur[b]] > hb: raise_on_read[(cls(b), hb, dims[cur[b]])] += 1
            continue
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1])
print('target frame dims at start of B:', tdim_at_B)
print('B target reads (tdim before, read-reg dim, tdim after, class): top 30')
for k, x in seq.most_common(30): print('  ', k, x)
print('read raises the read register frame (class, from, to):', raise_on_read.most_common(20))
print('distinct target readers per read register (in B):', Counter(len(set(x)) for x in helper_readers.values()))
