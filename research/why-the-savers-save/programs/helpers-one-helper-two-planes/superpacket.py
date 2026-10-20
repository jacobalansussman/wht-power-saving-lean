#!/usr/bin/env python3
"""Step 2b (variant I2, super-packet). Merge P2's face copy f into P1's face copy h when both are read by one common
reader A: h receives p2 directly from P2's other holder g2 (copy gate at time t_c, frame >= P1 + P2 + frame(g2, t_c)),
A reads h once. Exact necessary conditions checked:
 (i) a common reader A of f and h (cross-block planes), outside cleanup;
 (ii) every other (non-A) read of h happens before t_c (otherwise those readers receive the extra p2);
 (iii) frame: P1 + P2 + frame(g2, t_c) inside u_t^perp for every target that receives g2's content after t_c (exact
      marker, cut skipped, centres before B add their frame) -- t_c ranges over every gate time of g2 in
      [g2 holds p2, A's read of f);
 (iv) coefficient ratio: A's read coefficients of h and f (reported; any nonzero ratio is realisable by scaling the copy).
Usage: python3 -I superpacket.py FULL-WORD.json.gz survey.json superpacket.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']
gl = []
for name, g in segs: gl += [(name, e) for e in g]
paid, chains, _, _ = frame_replay(w)
def frame_at(q, t):
    U = chains[q][0][2]
    for tt, seg, UU in chains[q]:
        if tt <= t: U = UU
    return U
ctrf = {q: basis(U) for q, U in zip(w['centers'], w['center_frames'])}; Boff = off['B']
memo = {}
def S_after(q, t0):
    if (q, t0) in memo: return memo[(q, t0)]
    zz = {q: F(1)}; hit = set(); hc = set()
    for t in range(t0 + 1, len(gl)):
        name, e = gl[t]
        if name == 'cut': continue
        if e[0] != 'add':
            if e[1] in zz: zz[e[1]] *= F(e[2], e[3])
            continue
        a, b = e[1], e[2]
        if zz.get(b):
            zz[a] = zz.get(a, F(0)) + F(e[3], e[4]) * zz[b]
            if zz[a]:
                if v <= a < 2 * v: hit.add(a)
                if a in ctrf and t < Boff: hc.add(a)
    Sx = basis([x for x in range(1, 512) if all(dot(x, ports[u - v]) == 0 for u in hit) and all(contained((x,), ctrf[c]) for c in hc)])
    memo[(q, t0)] = Sx; return Sx
faces = {}
for r in S['planes']:
    for d in r['visitors']:
        if d['cls'] == 'face': faces[d['q']] = (tuple(r['plane']), r['block'], d, [x for x in r['visitors'] if x['out_spans_P'] and x['q'] != d['q']])
readsof = defaultdict(list)
for t, (name, e) in enumerate(gl):
    if e[0] == 'add' and name not in ('cleanup', 'compensation', 'precondition', 'cut'): readsof[e[2]].append((t, e[1], F(e[3], e[4])))
cnt = Counter(); passing = []
for f, (P2, b2, df, g2s) in faces.items():
    rf = [(t, A, c) for t, A, c in readsof[f] if t > df['t_in']]
    for hq, (P1, b1, dh, _) in faces.items():
        if hq == f or b1 == b2: continue
        rh = [(t, A, c) for t, A, c in readsof[hq] if t > dh['t_in']]
        common = {A for _, A, _ in rf} & {A for _, A, _ in rh}
        if not common: continue
        cnt['common reader'] += 1
        for A in common:
            tAf = min(t for t, a, _ in rf if a == A); tAh = min(t for t, a, _ in rh if a == A)
            tA = min(tAf, tAh)          # the merged single read by A happens no later than the earlier original read
            others = [t for t, a, _ in rh if a != A]
            # A's gates strictly between its two original reads (they would see the moved packet early)
            lo, hi = sorted((tAf, tAh)); between = [t for t in range(lo + 1, hi) if A in registers(gl[t][1]) and gl[t][0] not in ('cut',)]
            ok = False
            for g in g2s:
                for tc in range(max(g['t_in'], 0), tA):
                    e = gl[tc][1]
                    if g['q'] not in registers(e) and tc != max(g['t_in'], 0): continue
                    if any(t >= tc for t in others): continue
                    K = join(join(frame_at(hq, tc), P2), frame_at(g['q'], tc))
                    if contained(K, S_after(g['q'], tc)) and contained(K, S_after(hq, tc)):
                        ok = True; passing.append(dict(f=f, h=hq, A=A, g2=g['q'], tc=tc, tAf=tAf, tAh=tAh, A_gates_between=len(between), dimK=len(K), P1=list(P1), P2=list(P2))); break
                if ok: break
            cnt['passes (ii)+(iii)' if ok else 'fails'] += 1
            if ok: cnt['passes, A has no gate between its two reads'] += (len(between) == 0)
print(dict(cnt))
json.dump(dict(counts=dict(cnt), passing=passing), open(sys.argv[3], 'w'), indent=1)
