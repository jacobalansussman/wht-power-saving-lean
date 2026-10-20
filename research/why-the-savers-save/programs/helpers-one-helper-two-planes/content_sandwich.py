#!/usr/bin/env python3
"""Step 1d: content-based sandwich (the hub's C.1 upper side) versus the frame-propagation sandwich of pairs.py.
For each helper plane visitor f: inject a marker z into f at its plane arrival t_in and propagate exactly through the rest
of the program (a += c*b, scalings); the content consumers are the targets that ever hold z. S_content(f) = intersection of
their u_t^perp. Also records which targets are reached by frames but never by content, and through which gate types.
Then: ordered pairs (h, f), h a helper visitor of P1, with (a) P1 + P2 <= S_content(f)  [task's relaxed-H2 condition]
and (b) H1(h) <= S_content(f) [h's actual role-1 frame], cross-block vs same-block, plus exact max b-matching on (a).
Usage: python3 -I content_sandwich.py FULL-WORD.json.gz survey.json content.json"""
import sys, os, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib import *
from collections import Counter, defaultdict
w = load(sys.argv[1]); S = json.load(open(sys.argv[2])); v, h, R, n, full, start, end = geometry(w); ports = w['ports']
segs, comp = program(w); off = S['seg_offsets']
gl = []
for name, g in segs: gl += [(name, e) for e in g]
vis = {}; tin = {}
for r in S['planes']:
    for d in r['visitors']: vis[d['q']] = (tuple(r['plane']), r['block'], d['cls'], d['out_spans_P']); tin[d['q']] = d['t_in']
helpers = sorted(q for q in vis if q >= 2 * v)
H1 = {}
for name in ('A', 'B', 'tail'):
    for e in dict(segs)[name]:
        for q in registers(e): H1[q] = basis(e[-1])
def content_targets(f, t0):
    z = defaultdict(F); z[f] = F(1); hit = set()
    for t in range(t0 + 1, len(gl)):
        name, e = gl[t]
        if e[0] == 'add':
            a, b = e[1], e[2]
            if z.get(b):
                z[a] = z.get(a, F(0)) + F(e[3], e[4]) * z[b]
                if v <= a < 2 * v and z[a]: hit.add(a)
        else:
            a = e[1]
            if z.get(a): z[a] = z[a] * F(e[2], e[3])
    return hit
def inter(us):
    return basis([x for x in range(1, 512) if all(dot(x, u) == 0 for u in us)])
Sc = {}; ntg = {}
for f in helpers:
    T = content_targets(f, tin[f]); ntg[f] = len(T)
    Sc[f] = inter([ports[t - v] for t in T])
print('content consumers per helper visitor:', sorted(Counter(ntg.values()).items()))
print('dim S_content:', Counter(len(Sc[f]) for f in helpers), ' P_f inside S_content:', sum(1 for f in helpers if contained(vis[f][0], Sc[f])))
ca = Counter(); cb = Counter(); E = []
for hq in helpers:
    for f in helpers:
        if vis[hq][0] == vis[f][0]: continue
        rel = 'cross' if vis[hq][1] != vis[f][1] else 'same'
        if contained(join(vis[hq][0], vis[f][0]), Sc[f]):
            ca[(rel, vis[hq][2], vis[f][2])] += 1; E.append((hq, f))
        if contained(H1[hq], Sc[f]): cb[(rel, len(H1[hq]))] += 1
print('(a) P1+P2 <= S_content(f), by (block relation, h class, f class):', dict(ca))
print('(b) H1(h) <= S_content(f):', dict(cb))
# exact max b-matching on (a), plane level: each plane P2 can lose at most 1 visitor (it must keep 2 at P2), each helper h
# can take at most one extra role, and a helper that takes a role cannot itself be deleted. Bipartite: h -> f, plus
# constraints; solve as max-flow with node capacities (h used once as giver, f deleted once, a register cannot be both,
# at most one deletion per plane).
import itertools
src, snk = 's', 't'
cap = defaultdict(int); adj = defaultdict(set)
def add(a, b, c):
    cap[(a, b)] += c; adj[a].add(b); adj[b].add(a)
for hq, f in E:
    add(('h', hq), ('f', f), 1)
for hq in {a for a, b in E}: add(src, ('h', hq), 1)
for f in {b for a, b in E}: add(('f', f), ('P', vis[f][0]), 1)
for P in {vis[b][0] for a, b in E}: add(('P', P), snk, 1)
def maxflow():
    flow = 0
    while True:
        par = {src: None}; q = [src]
        for x in q:
            for y in adj[x]:
                if y not in par and cap[(x, y)] > 0: par[y] = x; q.append(y)
        if snk not in par: return flow
        y = snk
        while par[y] is not None: x = par[y]; cap[(x, y)] -= 1; cap[(y, x)] += 1; y = x
        flow += 1
mf = maxflow()
print('exact max matching on (a) [upper bound: ignores that a giver cannot also be deleted]:', mf)
json.dump(dict(ntargets={str(f): ntg[f] for f in helpers}, S_content={str(f): list(Sc[f]) for f in helpers}, pairs_a=E, maxmatch_a=mf), open(sys.argv[3], 'w'))
