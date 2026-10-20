#!/usr/bin/env python3
"""an-build: the END step of the Lean word (Carrier/Live.lean g8 Mi = PX + PY + PS Mi: after phase B the slot rows are
set from the inverse, i.e. slot := Mti[s, x] x + Mti[s, s'] s'; backward word: h8(Mt)^T = PX + PY + Mt^T PS, x += Mt[s,x]^T s).
It joins the x's (two-sided: end in s^perp, no s) with every slot whose final row depends on x.  A slot that is forced to
hold s at its end (forward propagation of s from the y gates of phase B and from the dirt step) cannot share the x frame:
its x-dependence must cross s.  Each crossing carries one linear combination, so the number of crossings per call is at
least  rank_Q( Mt[F, x] ),  F = the x-dependent slots forced to hold s at the end.  Exact (Fractions).
usage: python3 -I endstep.py FILE DIRTJSON"""
import gzip, json, sys
from fractions import Fraction as F
d = json.load(gzip.open(sys.argv[1])); dj = json.load(open(sys.argv[2]))
v, R = d['v'], d['R']; n = 2 * v + R
A, B = d['A'], d['B']; nA = len(A)
def adds(g):
    return [(t, g[2], F(a, b)) for t, a, b in g[3]] if g[0] == 'out' else [(g[2], s, F(a, b)) for s, a, b in g[3]]
def regs(g): return [g[2]] + [t for t, _, _ in g[3]] + (g[4] if g[0] == 'out' else [])
row = [{r: F(1)} for r in range(n)]
def addto(t, s, c):
    for k, val in list(row[s].items()):
        nv = row[t].get(k, 0) + c * val
        if nv: row[t][k] = nv
        else: row[t].pop(k, None)
for g in A:
    for t, s, c in adds(g): addto(t, s, c)
copies = {k: dict(row[s]) for k, s, _ in d['ret']}
for T, ent in enumerate(d['scat']['table']):
    for k, a, b in ent:
        for key, val in copies[k].items():
            nv = row[v + T].get(key, 0) + F(a, b) * val
            if nv: row[v + T][key] = nv
            else: row[v + T].pop(key, None)
for g in B:
    for t, s, c in adds(g): addto(t, s, c)
xdep = [r for r in range(2 * v, n) if any(k < v for k in row[r])]
# forced s at the end: forward from the y gates (B) only (the most permissive: no dirt-step joins, no totals rule)
s = [False] * n
for r in range(v, 2 * v): s[r] = True
for g in A + B:
    rs = regs(g)
    if any(s[r] for r in rs):
        for r in rs: s[r] = True
forced = [r for r in xdep if s[r]]
def rank(rows):
    piv = {}; rk = 0
    for rr in rows:
        vec = {k: val for k, val in rr.items() if k < v}
        while vec:
            p = min(vec)
            if p in piv:
                pv = piv[p]; c = vec[p] / pv[p]
                for k, val in pv.items():
                    nv = vec.get(k, 0) - c * val
                    if nv: vec[k] = nv
                    else: vec.pop(k, None)
            else:
                piv[p] = vec; rk += 1; break
    return rk
print('slots whose final row depends on x (joined with the x''s by the end step):', len(xdep), 'of', R)
print('  of them forced to hold s at their end by the y gates of phase B alone:', len(forced))
print('  rank of their x-dependence (crossings needed per call, lower bound):', rank([row[r] for r in forced]))
print('  rank of the x-dependence of all x-dependent slots:', rank([row[r] for r in xdep]))
