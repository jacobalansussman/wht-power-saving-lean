#!/usr/bin/env python3
"""g-small: for every packet plane of the U8 word, its visitors and how targets read them (times, frames, signs).
Usage: python3 -I u8_planes_readers.py CANONICAL-WORD.json.gz"""
import sys, json, gzip
from collections import Counter, defaultdict
c = json.loads(gzip.decompress(open(sys.argv[1], 'rb').read()))
F = c['frames']; v = c['v']; R = c['R']; n = 2*v + R; ports = c['ports']
def basis(b):
    bb = []
    for y in b:
        for z in bb: y = min(y, y ^ z)
        if y: bb.append(y)
    return sorted(bb, reverse=True)
B = [basis(f) for f in F]; dims = [len(b) for b in B]
def inside(a, b):
    for x in B[a]:
        y = x
        for z in B[b]: y = min(y, y ^ z)
        if y: return False
    return True
cls = lambda r: 'x' if r < v else ('y' if r < 2*v else 'h')
cur = list(c['start']); T = 0
events = []  # (time, phase, kind, a, b, frame_after_a, frame_after_b)
plane_of = {}  # register -> (plane basis tuple, time entered)
hist = defaultdict(list)  # register -> list of (time, frame)
for r in range(n): hist[r].append((-1, cur[r]))
def climb(r, f):
    if f is None or cur[r] == f: return
    if dims[f] <= dims[cur[r]] or not inside(cur[r], f): return
    cur[r] = f; hist[r].append((T, f))
    if dims[f] == 3 and r not in plane_of: plane_of[r] = (tuple(B[f]), T)
for r in range(n):
    if dims[cur[r]] == 3: plane_of[r] = (tuple(B[cur[r]]), -1)
for ph in ['A', 'centres', 'precondition', 'cut', 'compensation', 'B', 'tail', 'cleanup']:
    if ph == 'centres':
        for r, f in zip(c['centers'], c['center_frames']): climb(r, f)
        T += 1; continue
    for e in c.get(ph) or []:
        T += 1
        climb(e[1], e[-1])
        if e[0] == 'add': climb(e[2], e[-1]); events.append((T, ph, e[1], e[2], e[3], e[4], cur[e[1]], cur[e[2]]))
planes = defaultdict(list)
for r, (P, t) in plane_of.items(): planes[P].append((r, t))
print('planes', len(planes))
# for each plane visitor: who reads it (b == r) and who writes it, by class and phase; and the reader frames
rd = Counter(); rdph = Counter(); nreads = Counter(); tgt_reads_per_visitor = Counter()
readers_per_plane = Counter()
for P, vis in planes.items():
    tg = set()
    for r, t in vis:
        rs = [(T_, ph, a) for (T_, ph, a, b, *_ ) in events if b == r and T_ > t]
        tr = [a for _, _, a in rs if cls(a) == 'y']
        tgt_reads_per_visitor[(cls(r), len(set(tr)))] += 1
        tg |= set(tr)
        for _, ph, a in rs: rdph[(cls(r), ph, cls(a))] += 1
    readers_per_plane[len(tg)] += 1
print('visitor class, #distinct targets reading it after its plane visit:', sorted(tgt_reads_per_visitor.items()))
print('distinct targets reading a plane (any visitor):', sorted(readers_per_plane.items()))
print('reads of plane visitors by (visitor cls, phase, reader cls):', sorted(rdph.items()))
# multiplicity: how many times does each target read the same plane visitor (differential reading)?
mult = Counter()
for P, vis in planes.items():
    for r, t in vis:
        cnt = Counter(a for (T_, ph, a, b, *_ ) in events if b == r and cls(a) == 'y')
        for a, k in cnt.items(): mult[(cls(r), k)] += 1
print('(visitor cls, times one target reads that visitor):', sorted(mult.items()))
# one example plane in detail
P0 = sorted(planes)[0]
print('example plane', P0, 'visitors', planes[P0])
for r, t in planes[P0]:
    print(' visitor', r, cls(r), 'frames', [(tt, dims[f]) for tt, f in hist[r]])
    print('   reads by', [(T_, ph, a, cls(a), dims[fa]) for (T_, ph, a, b, num, den, fa, fb) in events if b == r][:40])
    print('   written by', [(T_, ph, b, cls(b)) for (T_, ph, a, b, *_ ) in events if a == r][:40])
