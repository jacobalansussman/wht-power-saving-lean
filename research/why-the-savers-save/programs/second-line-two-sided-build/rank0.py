#!/usr/bin/env python3
"""an-build: exact ranks of the two s-crossing maps of a GB call on a gcert/1 word (own code, Fractions).
  time 0 (transpose of the dirt gate g1):  slot q += sum_T Mt[y_T, q] y_T  for the dirt slots q that meet an x in phase A;
  end    (transpose of g8, h8^T):          x_t  += sum_q Mt[q, x_t] q      for the x-dependent slots forced to hold s.
The crossings needed per GB call are the ranks of these two maps (a basis crosses, free additions distribute it).
usage: python3 -I rank0.py FILE DIRTJSON"""
import gzip, json, sys
from fractions import Fraction as F
d = json.load(gzip.open(sys.argv[1])); dj = json.load(open(sys.argv[2]))
v, R = d['v'], d['R']; n = 2 * v + R
A, B = d['A'], d['B']
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
def rank(vecs):
    piv = {}; rk = 0
    for vec in vecs:
        vec = dict(vec)
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
Kx = dj['K_meet_x_in_A']
# functional of y delivered to slot q at time 0:  column q of Mt[y, .]  ->  vector over T
vecs0 = [{T: row[v + T][q] for T in range(v) if q in row[v + T]} for q in Kx]
r0 = rank(vecs0)
s = [False] * n
for r in range(v, 2 * v): s[r] = True
for g in A + B:
    rs = regs(g)
    if any(s[r] for r in rs):
        for r in rs: s[r] = True
Fs = [r for r in range(2 * v, n) if s[r] and any(k < v for k in row[r])]
vecsE = [{k: val for k, val in row[q].items() if k < v} for q in Fs]
rE = rank(vecsE)
print('time-0 crossing map (dirt slots meeting an x in A):', len(Kx), 'slots, rank', r0)
print('end crossing map (x-dependent slots forced to hold s):', len(Fs), 'slots, rank', rE)
print('crossings per GB call at this gate order >=', r0 + rE, '; per replica (4 GB calls) >=', 4 * (r0 + rE), 'against the gain 896')
json.dump(dict(r0=r0, rE=rE, nKx=len(Kx), nF=len(Fs)), open('rank0.json', 'w'))
