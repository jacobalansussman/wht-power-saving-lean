#!/usr/bin/env python3
"""Step 1g: the arrival sandwich with the copy time and copy source left free. For each third visitor f of P2 (the 168
face copies, which receive p2 by a copy gate) the merged helper h needs p2 at some time t_c before f's first reader
reads it, from some register g holding p2 (the other copy holder of P2; the proper-row visitor holds no p2). The copy gate runs at
K = P1 + P2 + frame(g, t_c); from then on h and g carry frames >= K. Necessary: K inside u_t^perp for every target that
later receives content from f's role (marker from f's arrival, f's own reach) or from g after t_c. Reach from a later
t_c is a subset, so every t_c in the window is tried. Counts (P1, f) cross-block with some (g, t_c) passing.
Usage: python3 -I latecopy.py FULL-WORD.json.gz survey.json latecopy.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']
gl = []
for name, g in segs: gl += [(name, e) for e in g]
paid, chains, off2, segs2 = frame_replay(w)
def frame_at(q, t):
    U = chains[q][0][2]
    for tt, seg, UU in chains[q]:
        if tt <= t: U = UU
    return U
ctrf = {q: basis(U) for q, U in zip(w['centers'], w['center_frames'])}; Boff = off['B']
MODE = os.environ.get('REACH', 'content')
def reach(q, t0):
    # 'content': exact marker coefficients (a register whose marker cancels to 0 stops carrying it) = weakest condition;
    # 'frame': every read propagates (frames climb at every read whatever the content)
    zz = {q: F(1)}; hit = set(); hc = set()
    for t in range(t0 + 1, len(gl)):
        name, e = gl[t]
        if name == 'cut': continue
        if e[0] != 'add':
            if e[1] in zz: zz[e[1]] *= F(e[2], e[3])
            continue
        a, b = e[1], e[2]
        if b in zz and (MODE == 'frame' or zz[b] != 0):
            zz[a] = zz.get(a, F(0)) + F(e[3], e[4]) * zz[b]
            if MODE == 'content' and zz[a] == 0: continue
            if v <= a < 2 * v: hit.add(a)
            if a in ctrf and t < Boff: hc.add(a)
    return hit, hc
def Sof(hit, hc):
    vecs = [x for x in range(1, 512) if all(dot(x, ports[t - v]) == 0 for t in hit) and all(contained((x,), ctrf[c]) for c in hc)]
    return basis(vecs)
planes = {tuple(r['plane']): r for r in S['planes']}
res = Counter(); best = []
for P2, r in planes.items():
    vs = r['visitors']
    for f in [d for d in vs if d['cls'] == 'face']:
        tf = f['t_in']
        # f's first read after arrival
        rd = [t for t in range(tf + 1, len(gl)) if gl[t][1][0] == 'add' and gl[t][1][2] == f['q'] and gl[t][0] != 'cut']
        t_end = rd[0] if rd else len(gl)
        Sf = Sof(*reach(f['q'], tf))
        cand = []
        for g in [d for d in vs if d['q'] != f['q'] and d['out_spans_P']]:   # g must hold p2 (the other copy holder)
            # times g holds p2: from its own arrival to t_end; try every change point of g's frame and each time g is touched
            ts = sorted({t for t in range(max(g['t_in'], 0), t_end) if gl[t][1][0] == 'add' and g['q'] in registers(gl[t][1])} | {max(g['t_in'], 0)})
            for tc in ts:
                Sg = Sof(*reach(g['q'], tc))
                Sx = basis([x for x in range(1, 512) if contained((x,), Sf) and contained((x,), Sg)])
                cand.append((g['q'], tc, frame_at(g['q'], tc), Sx))
        ok = 0
        for P1, r1 in planes.items():
            if r1['block'] == r['block']: continue
            for gq, tc, Fg, Sx in cand:
                if contained(join(join(P1, P2), Fg), Sx): ok += 1; best.append((f['q'], list(P1), gq, tc)); break
        res[ok > 0] += 1
        res['max_dim_S'] = max(res['max_dim_S'], max(len(c[3]) for c in cand) if cand else 0)
print('reach mode', MODE, ': face copies f with some cross-block P1 and some (g, t_c) passing:', res[True], 'of', res[True] + res[False], '; largest S over (g, t_c):', res['max_dim_S'])
json.dump(dict(passing=best), open(sys.argv[3], 'w'))
