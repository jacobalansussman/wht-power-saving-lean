#!/usr/bin/env python3
"""an-build: exact frame-legality test of TWO-SIDED slack on our lowering of the U8 local word (k-lower/out/U8X-parsed.json.gz,
gcert/1, our own data file; semantics of our checker tools/gx/gxcore.py chk_replay: x = 0..v-1, y = v..2v-1, slots after;
phase A never names a y; at the A/B cut the retained totals descend onto the y's (scatter); one common frame per gate;
frames monotone except the totals' descent at the cut).

Two-sided endpoint contract (complex engine; an-synth 4.1): y_T STARTS on a line <s> and x_S ENDS in s^perp, with s
orthogonal to every port.  Since the 112 ports span F_2^9 (checked below), s cannot lie in the 9-dim window; and the
scatter forces ONE s for all y's of the call (a total descends through a NESTED chain of its readers' frames; distinct
lines are not nested; the reader graph is connected: checked below).  So s is a 10th local coordinate, shared by all y's.

The s-content of a frame is then a single bit.  Propagation at the FIXED gate order (pr717 F2: legal iff G <= U):
  G(r,i) = s is in r's frame after gate i  (forward: y's start with s; a gate joins frames; totals hold s at the cut);
  A(r,i) = r's frame must avoid s at gate i (backward: every x avoids s at all times; a gate joins frames).
  violation: a gate with some participant in G and some participant in A.
Classes of slots: 'x-only' (some A, never G), 'y-only' (some G, never A), 'crossing' (A at an earlier gate, G later),
'free' (neither).  Same-run hosting: an X-type idle host never contains s (X ends in s^perp), a Y-type host contains s
from its start; so a crossing slot has no data-array host and must go to a bank."""
import gzip, json, sys, os
from collections import Counter
P = os.path.expanduser(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)), '../../data/inputs/U8X-parsed.json.gz'))
d = json.load(gzip.open(P))
h, v, R = d['h'], d['v'], d['R']; n = 2 * v + R
fr = d['frames']; dim = [len(f) for f in fr]
A, B = d['A'], d['B']; G_all = A + B; nA = len(A)
ret = d['ret']; scat = d['scat']['table']
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
def cls(r): return 'x' if r < v else 'y' if r < 2 * v else 's'
out = {}
# --- facts that force a single shared s outside the window
def rank2(vs):
    basis = []
    for u in vs:
        for b in basis:
            if u ^ b < u: u ^= b
        if u: basis.append(u); basis.sort(reverse=True)
    return len(basis)
out['port_span_rank'] = rank2(d['ports'])
readers = {}
for T, row in enumerate(scat):
    for k, a, b in row: readers.setdefault(k, set()).add(T)
# connectivity of the bipartite reader graph (y's -- totals)
par = list(range(v + len(ret)))
def find(a):
    while par[a] != a: par[a] = par[par[a]]; a = par[a]
    return a
for k, Ts in readers.items():
    for T in Ts: par[find(T)] = find(v + k)
comps = len({find(T) for T in range(v)})
out['readers_per_total'] = sorted(len(s) for s in readers.values()); out['reader_graph_components'] = comps
print('ports span rank', out['port_span_rank'], 'of', h, '| readers per total', out['readers_per_total'], '| components', comps)
# --- which registers appear where
first = {}; last = {}
for i, g in enumerate(G_all):
    for r in regs(g):
        first.setdefault(r, i); last[r] = i
xB = sorted({r for g in B for r in regs(g) if r < v}); xlastA = sum(1 for r in range(v) if last.get(r, -1) < nA)
print('x registers touched in phase B:', len(xB), '| x whose last gate is in A:', xlastA)
out['x_in_B'] = len(xB)
# --- forward G
Gs = [False] * n
for T in range(v): Gs[v + T] = True
Gat = []                                       # per gate: s present after the join
for i, g in enumerate(G_all):
    if i == nA:                                # the cut: totals must contain s (their descent ends in a y frame <s>+...)
        for _, s_, f in ret: Gs[s_] = True
    rs = regs(g)
    j = any(Gs[r] for r in rs)
    for r in rs: Gs[r] = Gs[r] or j
    Gat.append(j)
# --- backward A
Av = [False] * n
for r in range(v): Av[r] = True                # x avoid s up to their end
Aat = [False] * len(G_all)
for i in range(len(G_all) - 1, -1, -1):
    g = G_all[i]; rs = regs(g)
    j = any(Av[r] for r in rs)
    for r in rs: Av[r] = Av[r] or j
    Aat[i] = j
    # (the cut needs no backward rule: the totals' frames before the cut only have to lie below f + s)
viol = [i for i in range(len(G_all)) if Gat[i] and Aat[i]]
print('violating gates (some participant holds s, some must avoid it):', len(viol), 'of', len(G_all),
      '| in A:', sum(1 for i in viol if i < nA), 'in B:', sum(1 for i in viol if i >= nA))
out['violating_gates'] = len(viol)
# --- per-slot classes
gA = {r: [] for r in range(n)}; gG = {r: [] for r in range(n)}
Gs2 = [False] * n
for T in range(v): Gs2[v + T] = True
for i, g in enumerate(G_all):
    if i == nA:
        for _, s_, f in ret: Gs2[s_] = True
    for r in regs(g):
        if Aat[i]: gA[r].append(i)
        if Gat[i]: gG[r].append(i)
classes = Counter(); span = Counter(); cls_of = {}
start, final = d['start'], d['final']
for r in range(2 * v, n):
    a_, g_ = bool(gA[r]), bool(gG[r]) or r in {s_ for _, s_, _ in ret}
    if a_ and g_:
        c = 'crossing' if (not gG[r] or max(gA[r]) < min(gG[r] or [10**9])) else 'violating'
    elif a_: c = 'x-only'
    elif g_: c = 'y-only'
    else: c = 'free'
    cls_of[r] = c; classes[c] += 1
    span[c] += dim[final[r]] - dim[start[r]]
tot = {s_ for _, s_, _ in ret}
print('slot classes:', dict(classes), '| of them totals:', Counter(cls_of[s_] for s_ in tot))
print('in-call span (final dim - start dim) by class:', dict(span))
out['slot_classes'] = dict(classes); out['span_by_class'] = dict(span)
# how do crossing slots cross: last x-avoid gate in A or B, first s gate
cross = [r for r in cls_of if cls_of[r] == 'crossing']
print('crossing slots: last must-avoid gate in A for', sum(1 for r in cross if max(gA[r]) < nA), 'of', len(cross),
      '; first s gate in B for', sum(1 for r in cross if gG[r] and min(gG[r]) >= nA))
# direct x-y gates (illegal under two-sided)
xy = sum(1 for g in B if any(r < v for r in regs(g)) and any(v <= r < 2 * v for r in regs(g)))
print('gates joining an x and a y directly:', xy)
out['direct_xy_gates'] = xy
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'twosided_local.json'), 'w'), indent=1)
