#!/usr/bin/env python3
"""Step 1e: why the frame-propagation sandwich is tight. For each helper visitor f, S_frame(f) computed over gate windows
(A+B only; A+B+tail), and the number of frame-reached targets that are NOT content consumers. Plus, for the content-
compatible pairs (a) of content_sandwich.py, the earliest frame-reached target t with P1+P2 not inside u_t^perp, and the
path class (direct reader / via one intermediate / longer).
Usage: python3 -I frameprop_detail.py FULL-WORD.json.gz survey.json content.json frameprop.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); C = json.load(open(sys.argv[3])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']; D = dict(segs)
gl = []
for name in ('A', 'B', 'tail'):
    for i, e in enumerate(D[name]): gl.append((off[name] + i, name, registers(e), basis(e[-1])))
vis = {}; tin = {}
for r in S['planes']:
    for d in r['visitors']: vis[d['q']] = (tuple(r['plane']), r['block'], d['cls']); tin[d['q']] = d['t_in']
helpers = sorted(q for q in vis if q >= 2 * v)
def reach(f, t0, phases):
    rc = {f: (t0, None, 0)}  # reg -> (time, via, hops)
    for t, name, regs, U in gl:
        if t < t0 or name not in phases: continue
        a = [q for q in regs if q in rc]
        if a:
            src = a[0]
            for q in regs:
                if q not in rc: rc[q] = (t, src, rc[src][2] + 1)
    return rc
def inter(us): return basis([x for x in range(1, 512) if all(dot(x, u) == 0 for u in us)])
dims = defaultdict(Counter); nonc = Counter()
Sfr = {}
for f in helpers:
    for ph in (('A', 'B'), ('A', 'B', 'tail')):
        rc = reach(f, tin[f], ph)
        Ts = [q for q in rc if v <= q < 2 * v]
        Sx = inter([ports[q - v] for q in Ts]); dims[ph][len(Sx)] += 1
        if ph == ('A', 'B'): Sfr[f] = (Sx, rc)
for ph, c in dims.items(): print('S_frame dims, window', ph, dict(sorted(c.items())))
# pairs (a): first violating target in the A+B window
E = [tuple(x) for x in C['pairs_a']]
why = Counter(); ok = []
for hq, f in E:
    if vis[hq][1] == vis[f][1]: continue
    PP = join(vis[hq][0], vis[f][0]); Sx, rc = Sfr[f]
    bad = sorted((rc[q][0], q) for q in rc if v <= q < 2 * v and not inperp(PP, ports[q - v]))
    if not bad: ok.append((hq, f)); continue
    t, q = bad[0]; why[min(rc[q][2], 4)] += 1
print('cross-block content-compatible pairs:', sum(1 for a, b in E if vis[a][1] != vis[b][1]))
print('  surviving the A+B frame propagation:', len(ok))
print('  first violating target reached after hops (1 = direct partner of f, 4 = 4 or more):', dict(sorted(why.items())))
json.dump(dict(dims={str(k): dict(c) for k, c in dims.items()}, ok=ok, hops=dict(why)), open(sys.argv[4], 'w'))
# --- second test: the task's S_P = intersection over content consumers of BOTH packet copies (both holders' markers),
# and the identity of the first violating partner for the surviving content pairs
Sc = {int(k): tuple(x) for k, x in C['S_content'].items()}
holders = defaultdict(list)
for r in S['planes']:
    for d in r['visitors']:
        if d['out_spans_P']: holders[tuple(r['plane'])].append(d['q'])
gl2 = []
for name, g in segs: gl2 += [(name, e) for e in g]
def content_targets(f, t0):
    z = {f: F(1)}; hit = set()
    for t in range(t0 + 1, len(gl2)):
        name, e = gl2[t]
        if e[0] == 'add':
            a, b = e[1], e[2]
            if z.get(b):
                z[a] = z.get(a, F(0)) + F(e[3], e[4]) * z[b]
                if v <= a < 2 * v and z[a]: hit.add(a)
        elif z.get(e[1]): z[e[1]] = z[e[1]] * F(e[2], e[3])
    return hit
SP = {}
for P, hs in holders.items():
    T = set()
    for q in hs: T |= content_targets(q, tin[q])
    SP[P] = inter([ports[t - v] for t in T])
print('S_P over both copies: dims', dict(Counter(len(x) for x in SP.values())))
cnt = Counter(); E2 = []
for hq in helpers:
    for f in helpers:
        P1, P2 = vis[hq][0], vis[f][0]
        if vis[hq][1] == vis[f][1] or P1 == P2: continue
        if contained(join(P1, P2), SP[P2]): cnt[(vis[hq][2], vis[f][2], f in holders[P2])] += 1; E2.append((hq, f))
print('cross-block pairs with P1+P2 <= S_P2 (both copies), by (h cls, f cls, f is a copy holder):', dict(cnt))
why = Counter(); ok2 = []
for hq, f in E2:
    PP = join(vis[hq][0], vis[f][0]); Sx, rc = Sfr[f]
    bad = sorted((rc[q][0], q) for q in rc if v <= q < 2 * v and not inperp(PP, ports[q - v]))
    if not bad: ok2.append((hq, f)); continue
    t, q = bad[0]
    # walk back the path
    path = [q]; x = q
    while rc[x][1] is not None: x = rc[x][1]; path.append(x)
    path = path[::-1]
    def c(q):
        if q == f: return 'f'
        if q < v: return 'x'
        if q < 2 * v: return 'y'
        if q in holders[vis[f][0]]: return 'g2'
        return 'h'
    why[tuple(c(z) for z in path)] += 1
print('  surviving A+B frame propagation:', len(ok2))
print('  violating path shapes (f -> ... -> target):', why.most_common(8))
