"""d-halfread: g-small-2 variant screen. For each second copy h3 of the U8 FULL-WORD, look for a register r that never
sits at P (M copy, target, other helper) that could carry the packet instead of h3 by half reads: r takes alpha at
L_a before the merge and k*beta from b after it, then delivers p to h3's consumers. Tests (1) frames: r's frame just
before the merge lies inside L_a, its next frame contains L_a (stop inserted if needed), it reaches a frame containing
P after the merge; (2) delivery: the net coefficients with which registers read r while r would hold p (from its first
frame containing P to the start of cleanup) are proportional to the net coefficients h3's consumers take from h3.
Usage: python3 -I carriers.py WORD copies.json"""
import sys, json
from fractions import Fraction as F
from collections import Counter, defaultdict
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); assert not bad
C = json.load(open(sys.argv[2])); start, end = ends(w)
plane_visitors = {q for q in range(n) for T, U in chain[q] if len(U) == 3}
Tclean = min(T for T, (ph, e) in enumerate(prog) if ph == 'cleanup')
reads_of = defaultdict(list)
for T, (ph, e) in enumerate(prog):
    if e[0] == 'add' and ph != 'compensation': reads_of[e[2]].append((T, e[1], F(e[3], e[4])))
def frame_at(q, T):
    U = chain[q][0][1]
    for t, X in chain[q]:
        if t <= T: U = X
    return U
def next_climb(q, T):
    for t, X in chain[q]:
        if t > T: return t, X
    return None, None
cls = lambda q: 'x' if q < v else ('y' if q < 2*v else ('M' if q not in plane_visitors else 'h-plane'))
stat = Counter(); hits = []
for r in C:
    Tm, phm, a, b, k = r['merge']; Tc, phc, h3, q = r['copy']; P = tuple(r['P']); La = tuple(r['La'])
    need = Counter()
    for T, c, qq in reads_of[h3]:
        if Tc < T < r['Tend']: need[c] += qq
    need = {c: x for c, x in need.items() if x}
    best = None
    for c in range(n):
        if c in (a, b, h3) or c in plane_visitors and c >= 2*v: pass
        if c in (a, b, h3): continue
        F0 = frame_at(c, Tm-1)
        if not contained(F0, La): stat['frame before merge not inside L_a'] += 0; continue
        T1, F1 = next_climb(c, Tm-1)
        if F1 is None or not contained(La, F1): continue
        # first frame containing P after the merge
        TP = next((t for t, X in chain[c] if t >= Tm and contained(P, X)), None)
        if TP is None: continue
        got = Counter()
        for T, d, qq in reads_of[c]:
            if TP <= T < Tclean: got[d] += qq
        got = {d: x for d, x in got.items() if x}
        ok_frames = True
        if not need: continue
        ratio = None; match = set(got) == set(need)
        if match:
            ratios = {need[d]/got[d] for d in need}; match = len(ratios) == 1
        stat[('frame-ok', cls(c), 'delivery matches' if match else 'delivery differs')] += 1
        if match: hits.append(dict(P=list(P), h3=h3, carrier=c, cls=cls(c), need={str(d): str(x) for d, x in need.items()}))
for s, x in sorted(stat.items(), key=str): print(x, s)
print('hits', len(hits)); print(json.dumps(hits[:10]))
json.dump(hits, open('carriers.json', 'w'), indent=0)
