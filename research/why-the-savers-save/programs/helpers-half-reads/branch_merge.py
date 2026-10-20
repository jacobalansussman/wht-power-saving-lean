"""d-halfread variant check: can a second copy's readers read the merge holder a itself (no half read), i.e. does
h3's branch fit into a's branch?  For each copy: replace every read of h3 in the packet window by a read of a at
the same gate frame; a's chain then visits the union of both branches' frames in time order -- it must stay nested,
and a must hold the same value at each moved read as h3's packet part (a unchanged since the copy).
Usage: python3 -I branch_merge.py WORD copies.json"""
import sys, json
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); C = json.load(open(sys.argv[2]))
gl = {}
for T, (ph, e) in enumerate(prog):
    if ph != 'compensation':
        for q in regs(e): gl.setdefault(q, []).append(T)
st = Counter(); ok_list = []
for r in C:
    Tm, _, a, b, k = r['merge']; Tc, _, h3, q = r['copy']
    rd = [(T, c) for T, ph, c, cl, qq, d in r['reads'] if Tc < T < r['Tend']]
    a_writes = [T for T in gl[a] if T > Tc and prog[T][1][1] == a]
    first_aw = min(a_writes) if a_writes else 10**9
    if any(T > first_aw for T, c in rd): st['a rewritten before a copy read'] += 1; continue
    ev = [(t, U) for t, U in chain[a] if t > Tc] + [(T, basis(prog[T][1][-1])) for T, c in rd]
    ev.sort(); cur = chain[a][0][1]
    for t, U in chain[a]:
        if t <= Tc: cur = U
    nested = True
    for t, U in ev:
        if not contained(cur, U): nested = False; break
        cur = U
    if not nested: st['branches not nested (a would have to descend)'] += 1; continue
    st['fits'] += 1; ok_list.append((r['P'], h3))
print(dict(st)); print(ok_list[:10])
