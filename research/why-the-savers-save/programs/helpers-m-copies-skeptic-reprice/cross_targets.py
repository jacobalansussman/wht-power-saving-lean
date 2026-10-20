"""v-mcopy-code-reprice: can target-target reads at frames already in both chains carry the cut's packet content?
For each group G target t: other targets X sharing a chain frame with t (dim of shared frame), X's cut reads,
and whether t's post-cut packet vector lies in the span of the post-cut packet vectors of {X sharing a frame}."""
import sys, json
sys.path.insert(0, '.')
from wordlib import *
from fractions import Fraction as F
from collections import Counter, defaultdict
w = load(sys.argv[1]); gi = int(sys.argv[2]); v = w['report']['v']; n = 2*v + w['report']['R']
Gs = json.load(open('mc_groups.json')); allM = set(q for g in Gs for q in g['members'])
A = json.load(open('all_groups.json'))[gi]; T = A['targets']; mem = set(A['members'])
prog = program(w); start, end, chain, paid, runs = frames(w)
rows = [{q: F(1)} for q in range(n)]; P = []
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
        xs = sorted(s for s in rows[e[2]] if s < v)
        if len(xs) == 8: P += xs
    apply(e, rows)
    if ph == 'cut' and prog[t+1][0] != 'cut': postcut = [dict(r) for r in rows]
full = basis([1 << i for i in range(9)])
fr = {y: [tuple(U) for tag, U in chain[y] if U and tuple(U) != full] for y in range(v, 2*v)}
cutr = defaultdict(list)
for e in w['cut']: cutr[e[1]].append((e[2], F(e[3], e[4])))
def rank(M):
    M = [r[:] for r in M]; r = 0
    if not M: return 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        r += 1
    return r
vec = lambda y: [postcut[y].get(s, F(0)) for s in P]
for t in T:
    sh = defaultdict(list)
    for X in range(v, 2*v):
        if X == t: continue
        for G in set(fr[t]) & set(fr[X]): sh[X].append(len(G))
    outside = [X for X in sh if X not in T]
    dims = Counter(d for X in outside for d in sh[X])
    nz = [X for X in outside if any(postcut[X].get(s, 0) for s in P)]
    same = Counter(tuple(c for c, k in cutr[X]) for X in outside)
    base = [vec(X) for X in sh]
    r0 = rank(base); r1 = rank(base + [vec(t)])
    print('target', t, 'cut', [(c, str(k)) for c, k in cutr[t]], '| outside targets sharing a chain frame:', len(outside),
          'shared-frame dims', dict(dims), '| with nonzero group packet content post-cut:', len(nz),
          '| their cut-centre sets', dict(same), '| t in span of sharers (incl. group mates)?', r0 == r1, '(rank', r0, ')')
    base2 = [vec(X) for X in outside]
    print('     t in span of OUTSIDE sharers only?', rank(base2) == rank(base2 + [vec(t)]))
