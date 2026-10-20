"""d-halfread: why no M copy or target can carry a packet by half reads (first failing frame condition, over all
168 planes x candidate registers). Usage: python3 -I carriers_why.py WORD copies.json"""
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
stat = Counter(); detail = Counter()
for r in C:
    Tm = r['merge'][0]; P = tuple(r['P']); La = tuple(r['La']); Lb = tuple(r['Lb'])
    for c in range(v, n):
        if c in pv: continue
        kind = 'target' if c < 2*v else 'M copy'
        F0 = frame_at(c, Tm-1); T1, F1 = next_climb(c, Tm-1)
        capP = any(contained(P, X) for t, X in chain[c])
        if not capP: why = 'never reaches a frame containing P'
        elif not contained(F0, La):
            why = 'frame %d before the merge not inside L_a' % len(F0)
            if contained(F0, Lb): why += ' (inside L_b)'
        elif F1 is None or not contained(La, F1): why = 'next frame (rank %d) does not contain L_a' % (len(F1) if F1 else -1)
        else: why = 'frame-feasible'
        stat[(kind, why)] += 1
for s, x in sorted(stat.items()): print(x, s)
