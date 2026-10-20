"""d-halfread: can any non-plane helper (M copy) or target sit on ANY half line (L_a or L_b) of ANY of the 168 planes
before that plane's merge, with its next frame containing the line (a free rank-2 slot usable for a half read or for
relaying a target's 2-label term)?  Usage: python3 -I mcopy_lines.py WORD copies.json"""
import sys, json
from collections import Counter
sys.path.insert(0, __import__('os').path.dirname(__file__))
from wordlib import *
w = load(sys.argv[1]); v, R = w['report']['v'], w['report']['R']; n = 2*v+R
prog = program(w); chain, bad = frame_replay(w, prog); C = json.load(open(sys.argv[2]))
pv = {q for q in range(n) for T, U in chain[q] if len(U) == 3}
def frame_at(q, T):
    U = chain[q][0][1]
    for t, X in chain[q]:
        if t <= T: U = X
    return U
def next_climb(q, T):
    for t, X in chain[q]:
        if t > T: return t, X
    return None, None
st = Counter(); feas = []
for r in C:
    Tm = r['merge'][0]
    for side in ('La', 'Lb'):
        L = tuple(r[side])
        for c in range(v, n):
            if c in pv: continue
            kind = 'target' if c < 2*v else 'M copy'
            F0 = frame_at(c, Tm-1); T1, F1 = next_climb(c, Tm-1)
            if contained(F0, L) and F1 is not None and contained(L, F1):
                st[(kind, 'frame-feasible')] += 1; feas.append((r['P'], side, c))
            else: st[(kind, 'not')] += 1
print(dict(st)); print(feas[:10])
