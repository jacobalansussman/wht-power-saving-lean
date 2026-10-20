#!/usr/bin/env python3
"""Step 1c: exact frame-propagation sandwich for the merge 'helper h (a P1 visitor) takes over the role of helper f
(a P2 visitor, deleted)'. Equal-frame gates move BOTH registers to the gate frame, so lifting f's role by a frame H lifts
every register reached from f by a time-respecting gate path (A, B, tail; cut gates have no frame; cleanup is at the
full frame). Hence the merge is frame-legal only if H is inside S_f = intersection of u_t^perp over targets t reached
from f's role, intersected with the end frames of reached low-end helpers and the centre frames of centres reached
before B. Here H >= H1(h) = h's frame at the end of its own role (join of its A/B/tail gate frames), a lower bound
because h's residue still has to be cleared. Counts ordered pairs (h, f) with H1(h) <= S_f, by class and block relation,
with and without the timing condition (h's last A/B/tail gate before f reaches P2; f's earlier gates are dirt premixes).
Usage: python3 -I pairs.py FULL-WORD.json.gz survey.json pairs.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']; D = dict(segs)
gl = []  # (t, a, b, frame) for A, B, tail
for name in ('A', 'B', 'tail'):
    for i, e in enumerate(D[name]):
        gl.append((off[name] + i, registers(e), basis(e[-1])))
Bstart = off['B']
lowend = {int(q): basis(U) for q, U in w.get('helper_final', {}).items()}
ctr = {q: basis(U) for q, U in zip(w['centers'], w['center_frames'])}
def cls(q):
    if q < v: return 'x'
    if q < 2 * v: return 'y'
    hh = q - 2 * v
    return 'fac' if hh < 224 else ('face' if hh < 392 else 'M')
first = {}; last = {}; H1 = {}
for t, regs, U in gl:
    for q in regs:
        first.setdefault(q, t); last[q] = t; H1[q] = U   # monotone: frame after last gate = last gate frame
def downstream(f, t0):
    reach = {f: t0}
    for t, regs, U in gl:
        if t < t0: continue
        if any(q in reach for q in regs):
            for q in regs: reach.setdefault(q, t)
    return reach
vis = {}; tin = {}
for r in S['planes']:
    for d in r['visitors']: vis[d['q']] = (tuple(r['plane']), r['block'], d['cls'], d['out_spans_P']); tin[d['q']] = d['t_in']
helpers = [q for q in vis if q >= 2 * v]
Sf = {}; info = {}
for f in helpers:
    reach = downstream(f, tin[f])
    cons = []
    for q, t in reach.items():
        if v <= q < 2 * v: cons.append(perp(ports[q - v]))
        if q in lowend: cons.append(lowend[q])
        if q in ctr and t < Bstart: cons.append(ctr[q])
    # intersection of subspaces: vectors in all
    vecs = [x for x in range(1, 512) if all(contained((x,), C) for C in cons)]
    Sf[f] = basis(vecs)
    info[f] = dict(targets=sum(1 for q in reach if v <= q < 2 * v), lowends=sum(1 for q in reach if q in lowend), dimS=len(Sf[f]))
print('S_f dims over helper plane visitors (frame-propagation sandwich):', Counter(i['dimS'] for i in info.values()))
print('  of which with P_f inside S_f (sanity, must be all):', sum(1 for f in helpers if contained(vis[f][0], Sf[f])), '/', len(helpers))
print('H1 dims of helper plane visitors:', Counter(len(H1[q]) for q in helpers))
pairs = []; cnt = Counter()
for hq in helpers:
    for f in helpers:
        if hq == f or vis[hq][0] == vis[f][0]: continue
        if not contained(H1[hq], Sf[f]): continue
        timed = last[hq] < tin[f]
        rel = 'same-block' if vis[hq][1] == vis[f][1] else 'cross-block'
        cnt[(rel, timed)] += 1
        pairs.append(dict(h=hq, f=f, h_cls=vis[hq][2], f_cls=vis[f][2], rel=rel, timed=timed, dimH1=len(H1[hq]), dimS=len(Sf[f])))
print('ordered pairs with H1(h) <= S_f:', dict(cnt))
# looser: P1 + P2 <= S_f only (the relaxH2 condition, with the exact propagation sandwich)
c2 = Counter()
for hq in helpers:
    for f in helpers:
        if vis[hq][0] == vis[f][0]: continue
        if contained(join(vis[hq][0], vis[f][0]), Sf[f]): c2['cross' if vis[hq][1] != vis[f][1] else 'same'] += 1
print('looser test P1 + P2 <= S_f (relaxH2 condition, exact sandwich):', dict(c2))
json.dump(dict(info={str(k): x for k, x in info.items()}, pairs=pairs), open(sys.argv[3], 'w'))
