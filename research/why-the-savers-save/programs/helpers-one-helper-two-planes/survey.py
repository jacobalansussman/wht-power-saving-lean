#!/usr/bin/env python3
"""d-aboveplane step 1a: exact survey of the 168 packet planes of the U8 FULL-WORD (#942's unchanged local word).
For each rank-3 plane P: its block, its 4 source labels, its 3 visitors (arrival/leave times, frames before/after),
the clean content each visitor holds at arrival and at leave (source support), and every later read of a visitor by a
target (time, phase, gate frame) while it still holds a P-spanning row. Also the per-visitor 'S_read' = intersection
of u_t^perp over the targets that read it while it holds packet content (main phases only; compensation reads happen at
frame 0 on initial dirt and are excluded, see notes).
Usage: python3 -I survey.py FULL-WORD.json.gz survey.json"""
import sys, json
from collections import Counter, defaultdict
from fractions import Fraction as F
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
from lib import *

w = load(sys.argv[1]); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
paid, chains, off, segs = frame_replay(w)
assert dict(paid) == {int(k): a for k, a in w['report']['local_profile'].items()}, paid
ctr = set(w['centers'])
def cls(q):
    if q < v: return 'x'
    if q < 2 * v: return 'y'
    hh = q - 2 * v
    return 'fac' if hh < 224 else ('face' if hh < 392 else 'M')
# planes and visitors
vis = defaultdict(list)
for q in range(n):
    ch = chains[q]
    for i, (t, seg, U) in enumerate(ch):
        if len(U) == 3:
            tl = ch[i + 1][0] if i + 1 < len(ch) else None
            nxt = ch[i + 1][2] if i + 1 < len(ch) else None
            prev = ch[i - 1][2] if i > 0 else None
            vis[U].append(dict(q=q, cls=cls(q), t_in=t, t_out=tl, prev=list(prev) if prev is not None else None,
                               next=list(nxt) if nxt is not None else None, start=list(start[q]), end=list(end[q])))
planes = sorted(vis)
lab = {p: s for s, p in enumerate(ports)}
# blocks: sources of each plane; union-find
srcs = {P: [s for s in range(v) if contained((ports[s],), P)] for P in planes}
par = list(range(v))
def find(a):
    while par[a] != a: a = par[a]
    return a
for P in planes:
    for s in srcs[P][1:]: par[find(s)] = find(srcs[P][0])
blk = {P: find(srcs[P][0]) for P in planes}
bid = {b: i for i, b in enumerate(sorted(set(blk.values())))}
# value replay with snapshots of all plane visitors
watch = {d['q'] for P in planes for d in vis[P]}
segl = []
for name, g in segs: segl += [(name, e) for e in g]
rows = [{k: F(1)} for k in range(n)]
hist = defaultdict(list)   # q -> [(t, row)] after each write to q
reads = defaultdict(list)  # q -> [(t, seg, reader, frame, coeff)]
for t, (name, e) in enumerate(segl):
    if e[0] == 'add' and e[2] in watch: reads[e[2]].append((t, name, e[1], e[-1], [e[3], e[4]]))
    apply(e, rows)
    if e[1] in watch: hist[e[1]].append((t, dict(rows[e[1]])))
def row_at(q, t):
    r = {q: F(1)}
    for tt, rr in hist[q]:
        if tt <= t: r = rr
        else: break
    return r
def clean(r): return {k: x for k, x in r.items() if k < v}
def spans(r, P):
    S = [ports[k] for k in r if k < v]
    return len(S) > 0 and basis(S) == P
out = []
for P in planes:
    rec = dict(plane=list(P), block=bid[blk[P]], sources=srcs[P], visitors=[])
    for d in vis[P]:
        q = d['q']; tin = d['t_in']; tout = d['t_out']
        ain = clean(row_at(q, tin)); aout = clean(row_at(q, tout - 1 if tout is not None else 10**9))
        # later target reads while q's clean row spans P (packet-holding reads)
        rd = []
        for (t, seg, rr, fr, co) in reads[q]:
            if t < max(tin, 0): continue
            cr = clean(row_at(q, t - 1))
            if rr >= v and rr < 2 * v and seg != 'compensation':
                rd.append(dict(t=t, seg=seg, reader=rr, frame=list(basis(fr)) if fr else fr, holds_packet=spans(cr, P),
                               support=sorted(cr)))
        d.update(clean_in={str(k): str(x) for k, x in sorted(ain.items())}, clean_out={str(k): str(x) for k, x in sorted(aout.items())},
                 out_spans_P=spans(aout, P), target_reads=rd)
        rec['visitors'].append(d)
    out.append(rec)
json.dump(dict(planes=out, seg_offsets=off), open(sys.argv[2], 'w'))
print('planes', len(planes), 'blocks', len(bid), 'visitor classes', Counter(tuple(sorted(d['cls'] for d in r['visitors'])) for r in out))
print('visitors whose clean row spans P when leaving P', Counter(sum(d['out_spans_P'] for d in r['visitors']) for r in out))
print('target reads of packet-holding visitors by phase', Counter(x['seg'] for r in out for d in r['visitors'] for x in d['target_reads'] if x['holds_packet']))
print('target reads of non-packet visitor rows by phase', Counter(x['seg'] for r in out for d in r['visitors'] for x in d['target_reads'] if not x['holds_packet']))
