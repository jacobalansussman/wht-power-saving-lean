"""v-mcopy-code-reprice: cut and target-target structure around group G.
Shows: centres read by each group target in the cut; post-cut packet-source vectors C_t and their rank;
B target-target gates touching group targets (partners, frames); frames shared by group-target chains with ANY other target."""
import sys, json
sys.path.insert(0, '.')
from wordlib import *
from fractions import Fraction as F
from collections import Counter, defaultdict
w = load(sys.argv[1]); gi = int(sys.argv[2]); v = w['report']['v']; n = 2*v + w['report']['R']
Gs = json.load(open('mc_groups.json')); allM = set(q for g in Gs for q in g['members'])
A = json.load(open('all_groups.json'))[gi]; T = A['targets']; mem = set(A['members'])
prog = program(w); start, end, chain, paid, runs = frames(w)
rows = [{q: F(1)} for q in range(n)]; pk = []
postcut = None
for t, (ph, i, e) in enumerate(prog):
    if e[0] == 'add' and ph == 'B' and e[1] in mem and e[2] not in allM:
        xs = sorted(s for s in rows[e[2]] if s < v)
        if len(xs) == 8: pk.append(xs)
    apply(e, rows)
    if ph == 'cut' and prog[t+1][0] != 'cut': postcut = [dict(r) for r in rows]
P = [s for p in pk for s in p]
print('group', gi, 'targets', T, 'centres', w['centers'])
for y in T:
    cg = [(e[2], str(F(e[3], e[4]))) for e in w['cut'] if e[1] == y]
    print(' target', y, 'cut reads (centre, coef):', cg)
# C_t restricted to packet sources, per packet
for y in T:
    print('  C_%d by packet:' % y, [''.join('+' if postcut[y].get(s, 0) > 0 else ('-' if postcut[y].get(s, 0) < 0 else '0') for s in p) for p in pk])
# rank of the 4 vectors over Q
import itertools
M = [[postcut[y].get(s, F(0)) for s in P] for y in T]
def rank(M):
    M = [r[:] for r in M]; r = 0
    for c in range(len(M[0])):
        piv = next((i for i in range(r, len(M)) if M[i][c] != 0), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        for i in range(len(M)):
            if i != r and M[i][c] != 0:
                f = M[i][c] / M[r][c]; M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        r += 1
    return r
print(' rank of post-cut packet content of the 4 targets:', rank(M))
# B y-y gates touching group targets
for ph, i, e in prog:
    if ph == 'B' and e[0] == 'add' and v <= e[1] < 2*v and v <= e[2] < 2*v and (e[1] in T or e[2] in T):
        print('  B y-y', e[1], '<-', e[2], 'k', F(e[3], e[4]), 'frame dim', len(basis(e[-1])), 'both in group' if (e[1] in T and e[2] in T) else 'OUTSIDE GROUP')
full = basis([1 << i for i in range(9)])
fr = {y: set(tuple(U) for tag, U in chain[y] if U and tuple(U) != full) for y in range(v, 2*v)}
for y in T:
    sh = {y2: len(fr[y] & fr[y2]) for y2 in range(v, 2*v) if y2 != y and fr[y] & fr[y2]}
    print(' target', y, 'chain frame dims', [len(U) for tag, U in chain[y]], 'targets sharing a chain frame:', sh)
