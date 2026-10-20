#!/usr/bin/env python3
"""Step 1f (the decisive necessary condition). Merge: h (a visitor of plane P1) takes over the role of f (a helper visitor
of plane P2, deleted). In the merged word h arrives at P2's packet gate with frame >= P1, so the arrival gate of f at
P2 (time t_in(f)) runs at a frame >= P1 + P2, and from that time on BOTH f's register (now h) and its arrival partner
carry frames >= P1 + P2. Every target that later holds content from either of them (exact marker propagation through the
rest of the program, cancellations included) reads it at a frame >= P1 + P2, so P1 + P2 must lie in u_t^perp.
S_arr(f) = intersection of u_t^perp over those targets. This is order-independent in the sense that any word with the
same dataflow from these two registers has the same content consumers.
Cut gates (centre scatter, frame None) are skipped; a centre reached before B adds its centre frame.
Also checks (I9) sources as the third register: a source s visits its own plane P_s; for a plane P' with u_s in P'
(another plane of the source's block or a cross-block coincidence), P_s + P' must lie in S_arr of the P' visitor.
Usage: python3 -I arrival_sandwich.py FULL-WORD.json.gz survey.json arrival.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w)
gl = []
for name, g in segs: gl += [(name, e) for e in g]
vis = {}; tin = {}; planes = {}
for r in S['planes']:
    P = tuple(r['plane']); planes[P] = r
    for d in r['visitors']: vis[(d['q'], P)] = (r['block'], d['cls']); tin[(d['q'], P)] = d['t_in']
def markers(regs, t0):
    z = {q: {q: F(1)} for q in regs}  # per-marker propagation
    hit = set(); hitc.clear()
    zz = defaultdict(dict)
    for q in regs: zz[q][q] = F(1)
    for t in range(t0 + 1, len(gl)):
        name, e = gl[t]
        if name == 'cut':
            continue   # centre scatter: frame-free interface (target level 0); handled by the centre-frame constraint
        if e[0] == 'add':
            a, b = e[1], e[2]
            if zz.get(b):
                if a in ctrf and t < Boff: hitc.add(a)
                da = zz[a]
                for k, x in zz[b].items():
                    y = da.get(k, F(0)) + F(e[3], e[4]) * x
                    if y: da[k] = y
                    else: da.pop(k, None)
                if v <= a < 2 * v and da: hit.add(a)
        elif zz.get(e[1]):
            zz[e[1]] = {k: x * F(e[2], e[3]) for k, x in zz[e[1]].items()}
    return hit
hitc = set(); ctrf = {q: basis(U) for q, U in zip(w['centers'], w['center_frames'])}; Boff = S['seg_offsets']['B']
def inter(us): return basis([x for x in range(1, 512) if all(dot(x, u) == 0 for u in us)])
Sarr = {}; rows = []
for (q, P), (blk, c) in vis.items():
    t0 = tin[(q, P)]
    if t0 < 0: continue
    name, e = gl[t0]
    partner = [r for r in registers(e) if r != q]
    T = markers([q] + partner, t0)   # inject right after the arrival gate: both registers' content from then on
    Sx = inter([ports[t - v] for t in T])
    for cc in hitc: Sx = basis([x for x in range(1, 512) if contained((x,), Sx) and contained((x,), ctrf[cc])])
    Sarr[(q, P)] = Sx
    rows.append(dict(q=q, plane=list(P), block=blk, cls=c, partner=partner, ntargets=len(T), dimS=len(Sarr[(q, P)])))
print('arrival sandwich dims (all plane visitors with an arrival gate):', dict(Counter(r['dimS'] for r in rows)))
print('  by visitor class:', dict(Counter((r['cls'], r['dimS']) for r in rows)))
print('  P inside S_arr (sanity):', sum(1 for (q, P), Sx in Sarr.items() if contained(P, Sx)), '/', len(Sarr))
# pairs: h visitor of P1 (any class), f helper visitor of P2 (to delete)
cnt = Counter(); ok = []
for (hq, P1), (b1, c1) in vis.items():
    for (f, P2), (b2, c2) in vis.items():
        if f < 2 * v or P1 == P2 or hq == f or (f, P2) not in Sarr: continue
        if contained(join(P1, P2), Sarr[(f, P2)]):
            cnt[('cross' if b1 != b2 else 'same', c1, c2)] += 1; ok.append((hq, list(P1), f, list(P2)))
print('ordered pairs (h at P1, f at P2) with P1 + P2 <= S_arr(f):', dict(cnt))
# I9: sources. For each source s with plane P_s (its visit), planes P' containing u_s
lab = {}
src_plane = {}
for (q, P), (b, c) in vis.items():
    if q < v: src_plane[q] = P
same = Counter(); cross = []
for s, Ps in src_plane.items():
    for P2, r in planes.items():
        if P2 == Ps or not contained((ports[s],), P2): continue
        rel = 'same' if r['block'] == planes[Ps]['block'] else 'cross'
        Q = join(Ps, P2)
        fits = [d['q'] for d in r['visitors'] if d['q'] >= 2 * v and (d['q'], P2) in Sarr and contained(Q, Sarr[(d['q'], P2)])]
        same[(rel, len(Q), bool(fits))] += 1
        if rel == 'cross': cross.append((s, list(Ps), list(P2), len(Q), fits))
print('I9: (source, other plane P2 containing u_s) by (block relation, rank of P_s + P2, some helper visitor of P2 fits):', dict(same))
# planes lying inside a dim-4 S_arr other than P itself
inside = Counter()
for (f, P2), Sx in Sarr.items():
    inside[(len(Sx), sum(1 for P in planes if P != P2 and contained(P, Sx)))] += 1
print('(dim S_arr, number of OTHER packet planes inside S_arr):', dict(inside))
# free carriers: any helper a (all 504, incl. the 112 non-plane M copies) whose whole own role (frame after its last
# A/B/tail gate, H1(a)) fits inside S_arr(f) of a deletable helper visitor f, with a's role over before f arrives at P2
Dg = dict(segs); off = S['seg_offsets']; H1 = {}; last = {}
for name in ('A', 'B', 'tail'):
    for i, e in enumerate(Dg[name]):
        for q in registers(e): H1[q] = basis(e[-1]); last[q] = off[name] + i
planevis = {q for (q, P) in vis}
car = Counter(); carl = []
for a in range(2 * v, n):
    for (f, P2), Sx in Sarr.items():
        if f < 2 * v or a == f: continue
        if contained(H1.get(a, ()), Sx):
            timed = last.get(a, -1) < tin[(f, P2)]
            car[('plane' if a in planevis else 'nonplane', len(H1.get(a, ())), timed)] += 1; carl.append((a, f, list(P2), timed))
print('carriers a with H1(a) <= S_arr(f) [(a plane visitor?, dim H1(a), a done before f arrives)]:', dict(car))
json.dump(dict(rows=rows, pairs=ok, i9_cross=cross, carriers=carl), open(sys.argv[3], 'w'))
