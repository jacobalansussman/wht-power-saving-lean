#!/usr/bin/env python3
"""an-build: information paths of the BACKWARD reading of a gcert/1 local word (own code).  The backward word applies, in
the SAME time order, the transposed inverse of every single add (Work/GCert/Chain/Phased.lean bwd_scalar: b1 = MA^-T,
then the transposed scatter y -> copies -> retained slots, then b2 = MB^-T): an add  t += c*s  becomes  s -= c*t, so
information moves from the target of each add to its source.  We follow which registers carry information of the y's.
usage: python3 -I bwdpaths.py FILE"""
import gzip, json, sys
d = json.load(gzip.open(sys.argv[1]))
v, R = d['v'], d['R']; n = 2 * v + R
A, B = d['A'], d['B']
def adds(g):
    return [(t, g[2]) for t, _, _ in g[3]] if g[0] == 'out' else [(g[2], s) for s, _, _ in g[3]]
tot = [s for _, s, _ in d['ret']]
# forward reading, for comparison: does x information reach the y's only through the cut?
def flow(direction):
    has = [False] * n
    if direction == 'bwd':
        # y's carry their own information from the start; they are untouched in A
        for r in range(v, 2 * v): has[r] = True
    else:
        for r in range(v): has[r] = True
    reach_x_or_y = set()
    seq = [('A', g) for g in A] + [('cut', None)] + [('B', g) for g in B]
    for ph, g in seq:
        if ph == 'cut':
            if direction == 'bwd':
                ycut = any(has[r] for r in range(v, 2 * v))
                for t in tot: has[t] = has[t] or ycut        # y -> copies -> retained slots
            else:
                c = any(has[t] for t in tot)
                for r in range(v, 2 * v): has[r] = has[r] or c   # retained slots -> copies -> y
            continue
        for t_, s_ in adds(g):
            if direction == 'bwd':
                if has[t_]: has[s_] = True                   # s -= c t
            else:
                if has[s_]: has[t_] = True                   # t += c s
    return has
for direction in ('fwd', 'bwd'):
    has = flow(direction)
    tgt = range(v, 2 * v) if direction == 'fwd' else range(v)
    print(direction, ': registers of the other bank reached:', sum(has[r] for r in tgt), 'of', v)
# in the backward reading, which x's are reached ONLY through the cut (i.e. not by a B chain from a y)?
has = [False] * n
for r in range(v, 2 * v): has[r] = True
for g in B:
    for t_, s_ in adds(g):
        if has[t_]: has[s_] = True
print('bwd, without the cut (y information only through B adds): x reached', sum(has[r] for r in range(v)), 'of', v)
has = [False] * n
for t in tot: has[t] = True
for g in B:
    for t_, s_ in adds(g):
        if has[t_]: has[s_] = True
print('bwd, from the retained slots after the cut, through B adds: x reached', sum(has[r] for r in range(v)), 'of', v)
